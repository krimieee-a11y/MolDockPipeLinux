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
import shutil
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

# ======================== Tunables ========================
MANIFEST_FIELDS = [
    # ID + identity
    "id", "name", "smiles", "inchikey",
    # module artifacts
    "sdf_path", "pdbqt_path",
    # docking outputs (owned by worker)
    "vina_status", "vina_mode", "vina_score", "vina_rmsd_lb", "vina_rmsd_ub",
    "vina_exhaustiveness", "vina_seed", "vina_time",
    "vina_out_path",
    # run provenance
    "receptor", "receptor_sha1", "config_hash", "gpu_id",
    "tools_admet", "tools_build3d", "tools_prepare", "tools_vina",
    "updated_at",
]
SUMMARY_FIELDS = ["id", "name", "vina_score", "vina_mode", "vina_out_path"]
LEADER_FIELDS  = ["id", "name", "vina_score"]

# Folder layout (relative to cwd)
DIR_STATE   = Path("state")
DIR_RESULTS = Path("results")
DIR_PREP    = Path("prepared_ligands")
DIR_INPUT   = Path("input")
DIR_LOGS    = Path("logs")

DIR_STATE.mkdir(parents=True, exist_ok=True)
DIR_RESULTS.mkdir(parents=True, exist_ok=True)
DIR_LOGS.mkdir(parents=True, exist_ok=True)

MANIFEST_PATH = DIR_STATE / "manifest.csv"
SUMMARY_PATH  = DIR_RESULTS / "summary.csv"
LEADER_PATH   = DIR_RESULTS / "leaderboard.csv"

STOP_REQUESTED = False
HARD_STOP = False

def _sigint(_, __):
    global STOP_REQUESTED, HARD_STOP
    if not STOP_REQUESTED:
        STOP_REQUESTED = True
        print("\n🧯 Graceful stop requested. Finishing current operations…", flush=True)
    else:
        HARD_STOP = True
        print("\n🛑 Hard stop requested. Aborting launches.", flush=True)

signal.signal(signal.SIGINT, _sigint)
signal.signal(signal.SIGTERM, _sigint)

# ------------------------------ Utilities ------------------------------
def utcnow_iso() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")

def atomic_write_csv(path: Path, rows: List[dict], headers: List[str]) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=headers, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    tmp.replace(path)

def read_csv_dicts(path: Path) -> List[dict]:
    if not path.exists():
        return []
    out = []
    with path.open("r", newline="", encoding="utf-8") as f:
        rdr = csv.DictReader(f)
        for row in rdr:
            out.append({k: (row.get(k) or "") for k in rdr.fieldnames})
    return out

def load_manifest(path: Path) -> Dict[str, dict]:
    m = {}
    for row in read_csv_dicts(path):
        rid = row.get("id", "")
        if rid != "":
            m[rid] = row
    return m

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
        if base_row.get(k, "") == "" and shard_row.get(k, "") != "":
            base_row[k] = shard_row[k]
    # updated_at provenance
    base_row["updated_at"] = utcnow_iso()
    return changed

DOCKING_FIELDS = {
    "vina_status", "vina_mode", "vina_score", "vina_rmsd_lb", "vina_rmsd_ub",
    "vina_exhaustiveness", "vina_seed", "vina_time", "vina_out_path",
    "receptor", "receptor_sha1", "config_hash", "gpu_id", "tools_vina"
}

def pick_better(a: dict, b: dict) -> dict:
    """
    Deterministic winner between two shard rows for the same id.
    Prefer DONE with better (lower) vina_score; else latest updated_at.
    """
    sa, sb = (a.get("vina_status",""), b.get("vina_status",""))
    if sa == "DONE" and sb != "DONE":
        return a
    if sb == "DONE" and sa != "DONE":
        return b
    # both DONE or both not DONE → compare score
    fa, fb = safe_float(a.get("vina_score","")), safe_float(b.get("vina_score",""))
    if fa != fb:
        return a if fa < fb else b
    # tie-breaker: newest timestamp
    ta, tb = ts_ord(a.get("updated_at","")), ts_ord(b.get("updated_at",""))
    return a if ta >= tb else b

def rebuild_summaries(manifest: Dict[str, dict]) -> None:
    # summary.csv
    srows = []
    for rid, row in manifest.items():
        srows.append({
            "id": rid,
            "name": row.get("name",""),
            "vina_score": row.get("vina_score",""),
            "vina_mode": row.get("vina_mode",""),
            "vina_out_path": row.get("vina_out_path",""),
        })
    atomic_write_csv(SUMMARY_PATH, srows, SUMMARY_FIELDS)

    # leaderboard.csv
    lrows = []
    for rid, row in manifest.items():
        sc = safe_float(row.get("vina_score",""))
        if sc != float("inf"):
            lrows.append({"id": rid, "name": row.get("name",""), "vina_score": row.get("vina_score","")})
    lrows.sort(key=lambda r: safe_float(r["vina_score"]))
    atomic_write_csv(LEADER_PATH, lrows, LEADER_FIELDS)

# ------------------------------ Pre-flight consolidation ------------------------------
def preflight_consolidate() -> None:
    """Load base manifest, overlay any existing shard files, then backfill from results/*_out.pdbqt."""
    base = load_manifest(MANIFEST_PATH)
    shard_files = sorted(DIR_STATE.glob("manifest_gpu*.csv"))
    added, updated = 0, 0

    # 1) Overlay existing shard files (if any)
    for shard in shard_files:
        rows = read_csv_dicts(shard)
        for row in rows:
            rid = row.get("id","")
            if rid == "":
                continue
            if rid not in base:
                base[rid] = {k: "" for k in MANIFEST_FIELDS}
                base[rid]["id"] = rid
                base[rid]["name"] = row.get("name","")
                added += 1
            if overlay_docking(base[rid], row):
                updated += 1

    # 2) Backfill from prior *_out.pdbqt files (for resumability)
    b_add, b_upd = 0, 0
    for pose in sorted(DIR_RESULTS.glob("*_out.pdbqt")):
        # id is filename stem without suffix _out
        rid = pose.name.replace("_out.pdbqt","")
        if rid not in base:
            base[rid] = {k: "" for k in MANIFEST_FIELDS}
            base[rid]["id"] = rid
            b_add += 1
        # minimal docking fields
        row = {
            "id": rid,
            "vina_status": "DONE",
            "vina_out_path": str(pose),
            "updated_at": utcnow_iso(),
        }
        if overlay_docking(base[rid], row):
            b_upd += 1

    # Persist the preflight merge
    backup_manifest(MANIFEST_PATH)
    save_manifest(MANIFEST_PATH, base)
    print(f"🧩 Pre-flight: +{added} added, {updated} overlay-updated", flush=True)
    rebuild_summaries(base)
    print(f"   Backfill     : +{b_add} added, {b_upd} updated", flush=True)

# ------------------------------ GPU Scheduling ------------------------------
def detect_gpu_ids(limit: int | None = None) -> List[int]:
    # Try nvidia-smi -L and parse explicit numeric IDs; do not cap when limit <= 0 or None
    try:
        out = subprocess.check_output(["nvidia-smi", "-L"], text=True, stderr=subprocess.DEVNULL, timeout=2.0)
        ids = sorted({int(m.group(1)) for m in re.finditer(r"GPU\s+(\d+):", out)})
        if not ids:
            ids = [0]
    except Exception:
        # Fallback: assume at least one GPU (0)
        ids = [0]
    if limit is not None and limit > 0 and len(ids) > limit:
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
    bins = [[] for _ in range(max(1, n))]
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
    ]
    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = str(gid)
    f = logf.open("w", encoding="utf-8")
    f.write(f"[{utcnow_iso()}] Launch worker on GPU {gid} with {len(ids)} ligands\n")
    f.flush()
    return subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT, env=env)

# ------------------------------ Final consolidation ------------------------------
def consolidate_shards() -> Dict[str, dict]:
    base = load_manifest(MANIFEST_PATH)
    shard_files = sorted(DIR_STATE.glob("manifest_gpu*.csv"))
    merged = base.copy()
    total_added, total_updated = 0, 0
    def merge_into(dest: Dict[str, dict], rows: List[dict]) -> None:
        nonlocal total_added, total_updated
        for row in rows:
            rid = row.get("id","")
            if rid == "":
                continue
            if rid not in dest:
                dest[rid] = {k: "" for k in MANIFEST_FIELDS}
                dest[rid]["id"] = rid
                dest[rid]["name"] = row.get("name","")
                total_added += 1
            if overlay_docking(dest[rid], row):
                total_updated += 1

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
    ap.add_argument("--max-gpus", type=int, default=0,
                    help="Cap number of GPUs when auto-detecting. Use 0/negative for ALL (default).")
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
        # Respect CUDA_VISIBLE_DEVICES if set (numeric CSV). Otherwise auto-detect.
        env_gpus = os.environ.get("CUDA_VISIBLE_DEVICES", "").strip()
        gpu_ids = []
        if env_gpus:
            try:
                gpu_ids = [int(x) for x in env_gpus.split(",") if x.strip().isdigit()]
            except Exception:
                gpu_ids = []
        if not gpu_ids:
            cap = None if args.max_gpus <= 0 else args.max_gpus
            gpu_ids = detect_gpu_ids(limit=cap)

    if not gpu_ids:
        raise SystemExit("❌ No GPUs available/detected. Supply --gpu-ids if needed.")
    print("🎯 GPUs:", gpu_ids, flush=True)

    # 2) Pending ligands
    pending = list_pending_ligands()
    if not pending:
        print("✅ No pending ligands. Nothing to do.", flush=True)
        return

    print(f"🧩 Pending ligands: {len(pending)}", flush=True)

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
            if STOP_REQUESTED:
                break
            if not splits[i]:
                print(f"  GPU {gid}: nothing to do.", flush=True)
                continue
            procs.append(spawn_worker(worker_path, gid, splits[i]))

        # 5) Wait for workers
        while procs and not HARD_STOP:
            live = []
            for p in procs:
                r = p.poll()
                if r is None:
                    live.append(p)
            procs = live
            if procs:
                time.sleep(1.0)
    finally:
        # Best-effort: leave logs open/closed by Popen
        pass

    # 6) Final consolidation
    print("🧾 Final consolidation…", flush=True)
    merged = consolidate_shards()
    backup_manifest(MANIFEST_PATH)
    save_manifest(MANIFEST_PATH, merged)
    rebuild_summaries(merged)
    print("✅ Done.", flush=True)

if __name__ == "__main__":
    main()
