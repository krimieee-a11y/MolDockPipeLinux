#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module 4e (Controller) — OVERLAY-SAFE multi-GPU orchestrator for Vina-GPU (Linux)

What this script guarantees:
- Single-writer rule: only the controller writes state/manifest.csv, results/summary.csv, results/leaderboard.csv
- Workers write ONLY per-GPU shards: state/manifest_gpu<N>.csv (atomic .tmp then replace)
- Overlay-safe merge: preserves Modules 1–3 fields (ADMET/SDF/PDBQT, smiles, inchikey, tools_*)
  and overlays ONLY docking-owned fields from shards (vina_*, tools_vina, receptor_sha1, config_hash, updated_at)
- Deterministic conflict resolution across shards: prefer DONE with best (lowest) vina_score, else latest updated_at
- Pre-flight consolidation + backfill from existing *_out.pdbqt; final consolidation at end of run
- Atomic final writes + timestamped backups
- Idempotent rebuild of summary.csv & leaderboard.csv from merged manifest
"""

from __future__ import annotations
import argparse
import csv
import os
import re
import shlex
import signal
import subprocess
import sys
import time
from pathlib import Path
from datetime import datetime, timezone
from typing import Dict, List, Tuple

# ------------------------------ Paths & Files ------------------------------
BASE = Path(".").resolve()
DIR_PREP = BASE / "prepared_ligands"
DIR_RESULTS = BASE / "results"
DIR_STATE = BASE / "state"
DIR_LOGS = BASE / "logs"

DIR_STATE.mkdir(parents=True, exist_ok=True)
DIR_RESULTS.mkdir(parents=True, exist_ok=True)
DIR_LOGS.mkdir(parents=True, exist_ok=True)

FILE_MANIFEST = DIR_STATE / "manifest.csv"
FILE_SUMMARY = DIR_RESULTS / "summary.csv"
FILE_LEADERBOARD = DIR_RESULTS / "leaderboard.csv"

# ------------------------------ CSV Schema ------------------------------
MANIFEST_FIELDS = [
    "id","smiles","inchikey",
    "admet_status","admet_reason",
    "sdf_status","sdf_path","sdf_reason",
    "pdbqt_status","pdbqt_path","pdbqt_reason",
    "vina_status","vina_score","vina_pose","vina_reason",
    "config_hash","receptor_sha1","tools_rdkit","tools_meeko","tools_vina",
    "created_at","updated_at"
]

# Docking overlay policy
DOCKING_FIELDS = {
    "vina_status","vina_score","vina_pose","vina_reason",
    "receptor_sha1","tools_vina","config_hash","updated_at"
}

# ------------------------------ Utilities ------------------------------
def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00","Z")

def read_csv_dicts(path: Path) -> List[dict]:
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8") as f:
        return [dict(r) for r in csv.DictReader(f)]

def atomic_write_csv(path: Path, rows: List[dict], headers: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=headers)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in headers})
        f.flush()
        os.fsync(f.fileno())
    tmp.replace(path)

def load_manifest(path: Path) -> Dict[str, dict]:
    out: Dict[str, dict] = {}
    for r in read_csv_dicts(path):
        rid = r.get("id","")
        if not rid:
            continue
        row = {k: r.get(k, "") for k in MANIFEST_FIELDS}
        out[rid] = row
    return out

def save_manifest(path: Path, manifest: Dict[str, dict]) -> None:
    rows = [{k: v.get(k, "") for k in MANIFEST_FIELDS} for _, v in sorted(manifest.items())]
    atomic_write_csv(path, rows, MANIFEST_FIELDS)

def backup_manifest(path: Path) -> None:
    try:
        if path.exists():
            bdir = path.parent / "backups"
            bdir.mkdir(parents=True, exist_ok=True)
            ts = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
            (bdir / f"manifest.{ts}.csv").write_text(path.read_text(encoding="utf-8"), encoding="utf-8")
    except Exception:
        pass

def safe_float(x) -> float:
    try:
        return float(x)
    except Exception:
        return float("inf")

def ts_ord(x: str) -> float:
    try:
        s = (x or "").replace("Z","")
        return datetime.fromisoformat(s).timestamp()
    except Exception:
        return 0.0

def overlay_docking(base_row: dict, shard_row: dict) -> bool:
    changed = False
    for k in DOCKING_FIELDS:
        v = shard_row.get(k, "")
        if v != "" and base_row.get(k, "") != v:
            base_row[k] = v
            changed = True
    # helpful identity fill-ins if base is missing
    for k in ("smiles","inchikey","pdbqt_path"):
        if not base_row.get(k) and shard_row.get(k):
            base_row[k] = shard_row[k]
            changed = True
    if not base_row.get("created_at"):
        base_row["created_at"] = shard_row.get("created_at") or now_iso()
    # move updated_at forward
    if ts_ord(shard_row.get("updated_at","")) >= ts_ord(base_row.get("updated_at","")):
        base_row["updated_at"] = shard_row.get("updated_at") or now_iso()
        changed = True
    return changed

def pick_best_shard_rows(shard_rows: List[dict]) -> Dict[str, dict]:
    by_id: Dict[str, List[dict]] = {}
    for r in shard_rows:
        rid = (r.get("id") or "").strip()
        if rid:
            by_id.setdefault(rid, []).append(r)
    best: Dict[str, dict] = {}
    for rid, group in by_id.items():
        done = [g for g in group if str(g.get("vina_status","")).upper() == "DONE"]
        if done:
            best[rid] = min(done, key=lambda g: safe_float(g.get("vina_score","")))
        else:
            best[rid] = max(group, key=lambda g: ts_ord(g.get("updated_at","")))
    return best

# ------------------------------ Summary/Leaderboard ------------------------------
def rebuild_summaries(manifest: Dict[str, dict]) -> None:
    # summary
    s_rows = []
    for _, row in sorted(manifest.items()):
        sc = row.get("vina_score","")
        s_rows.append({
            "id": row.get("id",""),
            "inchikey": row.get("inchikey",""),
            "vina_score": sc,
            "pose_path": row.get("vina_pose",""),
            "created_at": row.get("updated_at","")
        })
    atomic_write_csv(FILE_SUMMARY, s_rows, ["id","inchikey","vina_score","pose_path","created_at"])

    # leaderboard: numeric sort ascending
    def _to_num(x):
        try: return float(x)
        except: return float("inf")
    ranked = sorted([r for r in s_rows if r.get("vina_score","") != ""], key=lambda r: _to_num(r["vina_score"]))
    leaders = []
    for i, r in enumerate(ranked, 1):
        leaders.append({
            "rank": i,
            "id": r["id"],
            "inchikey": r["inchikey"],
            "vina_score": r["vina_score"],
            "pose_path": r["pose_path"]
        })
    atomic_write_csv(FILE_LEADERBOARD, leaders, ["rank","id","inchikey","vina_score","pose_path"])

# ------------------------------ Preflight Consolidation ------------------------------
RES_RE = re.compile(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)", re.I)

def best_score_from_pose(pose_path: Path) -> str:
    try:
        txt = pose_path.read_text(errors="ignore")
        scores = [float(m.group(1)) for m in RES_RE.finditer(txt)]
        return f"{min(scores):.2f}" if scores else ""
    except Exception:
        return ""

def merge_into(dest: Dict[str, dict], src_rows: List[dict]) -> Tuple[int,int]:
    added, updated = 0, 0
    best_by_id = pick_best_shard_rows(src_rows)
    for rid, srow in best_by_id.items():
        if rid not in dest:
            row = {k: "" for k in MANIFEST_FIELDS}
            row["id"] = rid
            overlay_docking(row, srow)
            if not row.get("updated_at"):
                row["updated_at"] = now_iso()
            if not row.get("created_at"):
                row["created_at"] = now_iso()
            dest[rid] = row; added += 1
        else:
            base = dest[rid]
            if overlay_docking(base, srow):
                dest[rid] = base; updated += 1
    return added, updated

def backfill_from_results(manifest: Dict[str, dict]) -> Tuple[int,int]:
    added, updated = 0, 0
    for pose in DIR_RESULTS.glob("*_out.pdbqt"):
        lig_id = pose.name[:-len("_out.pdbqt")]
        pdbqt = DIR_PREP / f"{lig_id}.pdbqt"
        best = best_score_from_pose(pose)
        if lig_id in manifest:
            row = manifest[lig_id]
            changed = False
            if not row.get("vina_pose"):
                row["vina_pose"] = str(pose.resolve()); changed = True
            if not row.get("vina_score") and best:
                row["vina_score"] = best; changed = True
            if not row.get("vina_status"):
                row["vina_status"] = "DONE"; changed = True
            if changed:
                row["updated_at"] = now_iso(); updated += 1
        else:
            row = {k: "" for k in MANIFEST_FIELDS}
            row["id"] = lig_id
            if pdbqt.exists():
                row["pdbqt_path"] = str(pdbqt.resolve())
            row["vina_status"] = "DONE"
            row["vina_pose"] = str(pose.resolve())
            row["vina_score"] = best
            ts = now_iso()
            row["created_at"] = ts; row["updated_at"] = ts
            manifest[lig_id] = row; added += 1
    return added, updated

def preflight_consolidate() -> None:
    print("🧹 Pre-flight: consolidating leftover shards + results …", flush=True)
    merged = load_manifest(FILE_MANIFEST) if FILE_MANIFEST.exists() else {}
    # merge non-tmp shard files in state/
    shard_files = [p for p in sorted(DIR_STATE.glob("manifest_gpu*.csv")) if not str(p).endswith(".tmp")]
    total_added = total_updated = 0
    for shard in shard_files:
        rows = read_csv_dicts(shard)
        a,u = merge_into(merged, rows)
        total_added += a; total_updated += u
    # backfill from results
    b_add, b_upd = backfill_from_results(merged)

    backup_manifest(FILE_MANIFEST)
    save_manifest(FILE_MANIFEST, merged)
    rebuild_summaries(merged)
    print(f"   Shard overlay: +{total_added} added, {total_updated} updated", flush=True)
    print(f"   Backfill     : +{b_add} added, {b_upd} updated", flush=True)

# ------------------------------ GPU Scheduling ------------------------------
def detect_gpu_ids(limit: int | None = None) -> list[int]:
    """
    Auto-detect physical GPU IDs via `nvidia-smi -L`.
    If `limit` is None or <= 0, do NOT cap the number of GPUs.
    """
    try:
        out = subprocess.check_output(
            ["nvidia-smi", "-L"],
            text=True,
            stderr=subprocess.DEVNULL,
            timeout=2.0,
        )
        ids = sorted({int(m.group(1)) for m in re.finditer(r"GPU\s+(\d+):", out)})
        if not ids:
            ids = [0]
    except Exception:
        ids = [0]

    if limit and limit > 0 and len(ids) > limit:
        ids = ids[:limit]
    return ids

def parse_gpu_ids_arg(arg: str) -> List[int]:
    out = []
    for tok in arg.split(","):
        tok = tok.strip()
        if tok == "":
            continue
        try:
            out.append(int(tok))
        except Exception:
            pass
    return sorted(list(dict.fromkeys(out)))

def list_pending_ligands() -> List[str]:
    # any prepared_ligands/*.pdbqt that doesn't yet have results/<id>_out.pdbqt
    ids = []
    for p in sorted(DIR_PREP.glob("*.pdbqt")):
        lig_id = p.stem
        pose = DIR_RESULTS / f"{lig_id}_out.pdbqt"
        if not pose.exists():
            ids.append(lig_id)
    return ids

def round_robin_split(items: List[str], n: int) -> List[List[str]]:
    bins = [ [] for _ in range(max(1,n)) ]
    for i, it in enumerate(items):
        bins[i % len(bins)].append(it)
    return bins

def write_subset_file(gid: int, ids: List[str]) -> Path:
    subdir = DIR_STATE / "subsets"
    subdir.mkdir(parents=True, exist_ok=True)
    path = subdir / f"subset_gpu{gid}.list"
    path.write_text("\n".join(ids) + ("\n" if ids else ""), encoding="utf-8")
    return path

def spawn_worker(worker_path: Path, gid: int, ids: List[str]) -> subprocess.Popen:
    subset = write_subset_file(gid, ids)
    per_manifest = DIR_STATE / f"manifest_gpu{gid}.csv"
    logf = DIR_LOGS / f"gpu{gid}.log"
    cmd = [
        sys.executable, str(worker_path),
        "--gpu", str(gid),
        "--subset", str(subset),
        "--out-manifest", str(per_manifest),
        "--log", str(logf),
    ]
    print("Launching:", " ".join(shlex.quote(c) for c in cmd), flush=True)
    return subprocess.Popen(cmd)

def forward_sigint_to(procs: List[subprocess.Popen]) -> None:
    for p in procs:
        try:
            if p.poll() is None:
                p.send_signal(signal.SIGINT)
        except Exception:
            pass

# ------------------------------ Final Merge ------------------------------
def merge_per_gpu_manifests_all(base_manifest: Dict[str, dict]) -> Dict[str, dict]:
    merged = dict(base_manifest)
    shard_files = [p for p in sorted(DIR_STATE.glob("manifest_gpu*.csv")) if not str(p).endswith(".tmp")]
    for shard in shard_files:
        rows = read_csv_dicts(shard)
        merge_into(merged, rows)
    return merged

# ------------------------------ Main ------------------------------
def main():
    ap = argparse.ArgumentParser(description="Module 4e (Controller) — overlay-safe orchestrator")
    ap.add_argument("--worker", type=str, default="Module 4e (Worker) — SHARD SAFE.py",
                    help="Worker script path")
    ap.add_argument("--gpu-ids", type=str, default="",
                    help="Comma-separated GPU IDs (e.g., 0,1,2,3). If empty, auto-detect.")
    ap.add_argument(
    "--max-gpus",
    type=int,
    default=0,  # <= 0 means "use ALL detected GPUs" (new default)
    help="Cap number of GPUs when auto-detecting. "
         "Use 0 or negative to use ALL detected (default)."
)
)
    ap.add_argument("--dry-run", action="store_true", help="Plan only; do not launch workers")
    args = ap.parse_args()

    worker_path = Path(args.worker)
    if not worker_path.exists():
        raise SystemExit(f"❌ Worker script not found: {worker_path}")

    # 0) Pre-flight consolidation
    preflight_consolidate()

    # 1) GPUs
    if args.gpu_ids.strip():
        gpu_ids = parse_gpu_ids_arg(args.gpu_ids)
    else:
        # If CUDA_VISIBLE_DEVICES is set and no --gpu-ids given, prefer it
        env_gpus = os.environ.get("CUDA_VISIBLE_DEVICES", "").strip()
        if env_gpus:
            # Handle simple numeric CSV like "0,1,2"; ignore UUID entries
            try:
                gpu_ids = [int(x) for x in env_gpus.split(",") if x.strip().isdigit()]
            except Exception:
                gpu_ids = []  # fall through to detector
        if not env_gpus or not gpu_ids:
            cap = None if args.max_gpus <= 0 else args.max_gpus
            gpu_ids = detect_gpu_ids(limit=cap)

if not gpu_ids:
    raise SystemExit("❌ No GPUs available/detected. Supply --gpu-ids if needed.")

print("🎯 GPUs:", gpu_ids, flush=True)


    # 2) Pending ligands
    pending = list_pending_ligands()
    if not pending:
        print("✅ Nothing to dock — all ligands already have outputs. Finalizing summaries…", flush=True)
        merged = load_manifest(FILE_MANIFEST) if FILE_MANIFEST.exists() else {}
        rebuild_summaries(merged)
        return

    # 3) Split work
    splits = round_robin_split(pending, len(gpu_ids))
    for i, gid in enumerate(gpu_ids):
        print(f"  GPU {gid}: {len(splits[i])} ligands", flush=True)

    if args.dry_run:
        print("🧪 Dry run — not launching workers.", flush=True)
        return

    # 4) Launch workers
    procs: List[subprocess.Popen] = []
    try:
        for i, gid in enumerate(gpu_ids):
            chunk = splits[i]
            if not chunk:
                print(f"GPU {gid} — no ligands assigned.", flush=True)
                continue
            p = spawn_worker(worker_path, gid, chunk)
            procs.append(p)

        # 5) Monitor workers
        while procs:
            alive = [p for p in procs if p.poll() is None]
            if not alive:
                break
            time.sleep(0.5)

    except KeyboardInterrupt:
        print("\n⏹️  Ctrl+C — asking workers to stop…", flush=True)
        forward_sigint_to(procs)

    finally:
        # 6) Final merge + backfill + summaries
        base = load_manifest(FILE_MANIFEST) if FILE_MANIFEST.exists() else {}
        merged = merge_per_gpu_manifests_all(base)
        b_add, b_upd = backfill_from_results(merged)
        backup_manifest(FILE_MANIFEST)
        save_manifest(FILE_MANIFEST, merged)
        rebuild_summaries(merged)
        rcs = [p.poll() for p in procs]
        print("✅ Controller done. Worker return codes:", rcs, flush=True)
        print(f"   Manifest : {FILE_MANIFEST}", flush=True)
        print(f"   Summary  : {FILE_SUMMARY}", flush=True)
        print(f"   Leaders  : {FILE_LEADERBOARD}", flush=True)

if __name__ == "__main__":
    if sys.platform.startswith("linux"):
        os.environ.setdefault("PYTHONUNBUFFERED","1")
    main()
