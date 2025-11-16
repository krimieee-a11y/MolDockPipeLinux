#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module 1 (BOILED-Egg YOLK Gate): ADMET Screening (Pipeline-Compatible)

- Reads input/input.csv (expects: id, smiles)
- Computes descriptors using RDKit (MW, HBD, HBA, ROTB, TPSA, WLOGP)
- Applies Lipinski, Veber, AND BOILED-Egg YOLK requirement
- Writes:
    output/admet.csv
    state/admet_pass.list
    state/admet_fail.list
- Fully compatible with the unified manifest schema used by:
    Module 2, Module 3, Module 4b GPU, Module 4e Multi-GPU

PASS RULE (strict):
    Lipinski OK
    Veber OK
    BOILED-Egg region == YOLK

Run:
    python "Module 1 (BOILED-Egg).py"
"""

from __future__ import annotations
import csv, time, hashlib
from pathlib import Path
from typing import Dict, Any, Tuple

# ------------------------ Paths ------------------------
ROOT = Path(__file__).resolve().parent
DIR_INPUT  = ROOT / "input"
DIR_OUTPUT = ROOT / "output"
DIR_STATE  = ROOT / "state"

FILE_INPUT    = DIR_INPUT / "input.csv"
FILE_ADMET    = DIR_OUTPUT / "admet.csv"
FILE_PASS     = DIR_STATE / "admet_pass.list"
FILE_FAIL     = DIR_STATE / "admet_fail.list"
FILE_MANIFEST = DIR_STATE / "manifest.csv"

# ------------------------ Optional deps ------------------------
try:
    from rdkit import Chem
    from rdkit.Chem import Crippen, rdMolDescriptors, Descriptors, Lipinski
    RDKit_OK = True
except Exception:
    RDKit_OK = False
    Chem = None

# ---------------------- Manifest Schema ------------------------
MANIFEST_FIELDS = [
    "id","smiles","inchikey",
    "admet_status","admet_reason",
    "sdf_status","sdf_path","sdf_reason",
    "pdbqt_status","pdbqt_path","pdbqt_reason",
    "vina_status","vina_score","vina_pose","vina_reason",
    "config_hash","receptor_sha1","tools_rdkit","tools_meeko","tools_vina",
    "created_at","updated_at"
]

# ------------------------ Utilities ------------------------
def ensure_dirs():
    for d in (DIR_INPUT, DIR_OUTPUT, DIR_STATE):
        d.mkdir(parents=True, exist_ok=True)

def now():
    return time.strftime("%Y-%m-%d %H:%M:%S")

def read_input_rows() -> list[dict]:
    if not FILE_INPUT.exists():
        raise FileNotFoundError("Missing input/input.csv (must contain id,smiles)")
    rows = []
    with FILE_INPUT.open("r", newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            if row.get("id") and row.get("smiles"):
                rows.append({"id": row["id"].strip(), "smiles": row["smiles"].strip()})
    return rows

# --------------------- Descriptor Logic ----------------------
def compute_descriptors(smiles: str) -> Tuple[Dict[str, Any], str]:
    if not RDKit_OK:
        return ({"valid": False, "error": "rdkit_not_available"}, "")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ({"valid": False, "error": "invalid_smiles"}, "")
    try:
        mw    = Descriptors.ExactMolWt(mol)
        hbd   = Lipinski.NumHDonors(mol)
        hba   = Lipinski.NumHAcceptors(mol)
        rotb  = Lipinski.NumRotatableBonds(mol)
        tpsa  = rdMolDescriptors.CalcTPSA(mol)
        wlogp = Crippen.MolLogP(mol)
        rings = Lipinski.RingCount(mol)
        inchikey = Chem.MolToInchiKey(mol)
        return ({
            "valid": True, "mw": mw, "hbd": hbd, "hba": hba,
            "rotb": rotb, "tpsa": tpsa, "wlogp": wlogp, "rings": rings
        }, inchikey)
    except Exception:
        return ({"valid": False, "error": "rdkit_error"}, "")

# --------------------- ADMET RULES -----------------------
def rule_lipinski(d):
    return (d["mw"] <= 500 and d["hbd"] <= 5 and d["hba"] <= 10 and d["wlogp"] <= 5)

def rule_veber(d):
    return (d["tpsa"] <= 140 and d["rotb"] <= 10)

def boiled_egg_region(d):
    if not d.get("valid"):
        return "GREY"
    tpsa = float(d["tpsa"])
    wlogp = float(d["wlogp"])
    if (tpsa < 79.0 and 0.4 <= wlogp <= 6.0):
        return "YOLK"
    if (tpsa <= 130.0 and -0.4 <= wlogp <= 5.6):
        return "WHITE"
    return "GREY"

def decide_pass(desc):
    if not desc.get("valid"):
        return ("FAIL", "invalid")

    if not rule_lipinski(desc):
        return ("FAIL", "lipinski")
    if not rule_veber(desc):
        return ("FAIL", "veber")

    egg = boiled_egg_region(desc)
    if egg != "YOLK":
        return ("FAIL", "boiled_egg:not_yolk")

    return ("PASS", "ok")

# --------------------- Manifest I/O -----------------------
def load_manifest() -> dict[str, dict]:
    if not FILE_MANIFEST.exists():
        return {}
    out = {}
    with FILE_MANIFEST.open("r", newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            fixed = {k: row.get(k, "") for k in MANIFEST_FIELDS}
            out[fixed["id"]] = fixed
    return out

def save_manifest(data: dict[str, dict]):
    with FILE_MANIFEST.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=MANIFEST_FIELDS)
        w.writeheader()
        for lid in sorted(data.keys()):
            row = data[lid]
            out = {k: row.get(k, "") for k in MANIFEST_FIELDS}
            w.writerow(out)

# ---------------------------- MAIN ----------------------------
def main():
    ensure_dirs()
    rows = read_input_rows()
    manifest = load_manifest()

    created_ts = now()
    admet_rows = []
    pass_ids, fail_ids = [], []

    for r in rows:
        lid = r["id"]
        smiles = r["smiles"]

        desc, inchikey = compute_descriptors(smiles)
        decision, reason = decide_pass(desc)

        egg = boiled_egg_region(desc) if desc.get("valid") else "GREY"

        # ------------ Write ADMET output row ------------
        admet_rows.append({
            "id": lid,
            "smiles": smiles,
            "inchikey": inchikey,
            "mw": desc.get("mw",""),
            "hbd": desc.get("hbd",""),
            "hba": desc.get("hba",""),
            "rotb": desc.get("rotb",""),
            "tpsa": desc.get("tpsa",""),
            "wlogp": desc.get("wlogp",""),
            "rings": desc.get("rings",""),
            "boiled_egg": egg,
            "admet_decision": decision,
            "admet_reason": reason
        })

        # ------------ PASS/FAIL list ------------
        if decision == "PASS":
            pass_ids.append(lid)
        else:
            fail_ids.append(lid)

        # ------------ Manifest update ------------
        m = manifest.get(lid, {k:"" for k in MANIFEST_FIELDS})
        m["id"] = lid
        m["smiles"] = smiles
        m["inchikey"] = inchikey
        m["admet_status"] = decision
        m["admet_reason"] = reason

        if not m.get("created_at"):
            m["created_at"] = created_ts
        m["updated_at"] = now()

        # Mark RDKit availability
        m["tools_rdkit"] = "RDKit OK" if RDKit_OK else "RDKit MISSING"

        manifest[lid] = m

    # ---------------- Write outputs -----------------
    # ADMET CSV
    with FILE_ADMET.open("w", newline="", encoding="utf-8") as f:
        fields = [
            "id","smiles","inchikey","mw","hbd","hba","rotb","tpsa","wlogp","rings",
            "boiled_egg","admet_decision","admet_reason"
        ]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in admet_rows:
            w.writerow(row)

    # PASS/FAIL lists
    with FILE_PASS.open("w", encoding="utf-8") as f:
        f.write("\n".join(pass_ids))
    with FILE_FAIL.open("w", encoding="utf-8") as f:
        f.write("\n".join(fail_ids))

    # Manifest
    save_manifest(manifest)

    print(f"✅ ADMET done. Rows={len(admet_rows)} PASS={len(pass_ids)} FAIL={len(fail_ids)}")
    print(f"   ADMET CSV : {FILE_ADMET}")
    print(f"   Manifest  : {FILE_MANIFEST}")

if __name__ == "__main__":
    main()
