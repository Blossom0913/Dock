#!/usr/bin/env python3
"""
Workflow: split SMILES into PDB, convert PDB to PDBQT, dock with AutoDock Vina.

Steps:
 1) SMI -> multi PDB (one per molecule) via Open Babel (--gen3d, -h)
 2) PDB -> PDBQT via Open Babel (-opdbqt)
 3) Dock each ligand PDBQT with receptor/box from config.txt (Vina)

Outputs:
 - <workdir>/pdb/*.pdb
 - <workdir>/pdbqt/*.pdbqt
 - <workdir>/out_vina/*_out.pdbqt and per-ligand logs
 - <workdir>/summary.txt with time consumption

Example run (vina env):
 conda run -n vina python workflow_from_smi.py \
   --smi \
   "/home/cxt/Dock/data/SMILES format-20250819T054226Z-1-002/SMILES format/Genesis.smi" \
   --workdir from_smi --limit 10 --cpu 8 --exhaustiveness 8 --num-modes 9
"""

from __future__ import annotations

import argparse
import glob
import os
import shlex
import subprocess
import sys
import time
from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple


CONFIG_FILE_NAME = "config.txt"


@dataclass
class DockBox:
    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float


def parse_simple_kv_config(path: str) -> Dict[str, str]:
    key_to_value: Dict[str, str] = {}
    if not os.path.exists(path):
        raise FileNotFoundError(f"Config file not found: {path}")
    with open(path, "r", encoding="utf-8") as file:
        for line in file:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if "=" not in stripped:
                continue
            key, value = stripped.split("=", 1)
            key_to_value[key.strip()] = value.strip()
    return key_to_value


def load_box_from_config(config_path: str) -> Tuple[str, DockBox]:
    kv = parse_simple_kv_config(config_path)
    required_keys = [
        "receptor",
        "center_x",
        "center_y",
        "center_z",
        "size_x",
        "size_y",
        "size_z",
    ]
    missing = [k for k in required_keys if k not in kv]
    if missing:
        raise ValueError(f"Missing keys in {config_path}: {', '.join(missing)}")

    receptor_rel = kv["receptor"]
    receptor_path = receptor_rel
    if not os.path.isabs(receptor_path):
        receptor_path = os.path.abspath(os.path.join(os.path.dirname(config_path), receptor_rel))
    if not os.path.exists(receptor_path):
        raise FileNotFoundError(f"Receptor not found: {receptor_path}")

    box = DockBox(
        center_x=float(kv["center_x"]),
        center_y=float(kv["center_y"]),
        center_z=float(kv["center_z"]),
        size_x=float(kv["size_x"]),
        size_y=float(kv["size_y"]),
        size_z=float(kv["size_z"]),
    )
    return receptor_path, box


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def run_obabel_split_smi_to_pdb(smi_path: str, out_pdb_dir: str) -> Tuple[int, str, str]:
    ensure_dir(out_pdb_dir)
    # Use numbered output files with padding; --gen3d to generate 3D conformers; -h to add hydrogens.
    out_pattern = os.path.join(out_pdb_dir, "mol_%06d.pdb")
    cmd = [
        "obabel",
        "-ismi",
        smi_path,
        "-opdb",
        "-O",
        out_pattern,
        "-m",
        "--gen3d",
        "-h",
    ]
    completed = subprocess.run(cmd, capture_output=True, text=True)
    return completed.returncode, completed.stdout, completed.stderr


def discover_files(patterns: Iterable[str]) -> List[str]:
    files: List[str] = []
    for pattern in patterns:
        files.extend(glob.glob(pattern))
    files = [f for f in files if os.path.isfile(f)]
    files.sort()
    return files


def obabel_pdb_to_pdbqt(pdb_path: str, pdbqt_path: str) -> Tuple[int, str, str]:
    ensure_dir(os.path.dirname(pdbqt_path))
    cmd = [
        "obabel",
        "-ipdb",
        pdb_path,
        "-opdbqt",
        "-O",
        pdbqt_path,
        "-h",
    ]
    completed = subprocess.run(cmd, capture_output=True, text=True)
    return completed.returncode, completed.stdout, completed.stderr


def vina_dock(
    receptor_path: str,
    ligand_pdbqt_path: str,
    box: DockBox,
    out_dir: str,
    exhaustiveness: int,
    num_modes: int,
    cpu: int,
) -> Tuple[int, str, str, str]:
    ensure_dir(out_dir)
    base = os.path.splitext(os.path.basename(ligand_pdbqt_path))[0]
    out_pdbqt = os.path.join(out_dir, f"{base}_out.pdbqt")
    cmd = [
        "vina",
        "--receptor",
        receptor_path,
        "--ligand",
        ligand_pdbqt_path,
        "--center_x",
        str(box.center_x),
        "--center_y",
        str(box.center_y),
        "--center_z",
        str(box.center_z),
        "--size_x",
        str(box.size_x),
        "--size_y",
        str(box.size_y),
        "--size_z",
        str(box.size_z),
        "--exhaustiveness",
        str(exhaustiveness),
        "--num_modes",
        str(num_modes),
        "--out",
        out_pdbqt,
    ]
    if cpu and cpu > 0:
        cmd.extend(["--cpu", str(cpu)])
    completed = subprocess.run(cmd, capture_output=True, text=True)
    return completed.returncode, out_pdbqt, completed.stdout, completed.stderr


def write_text(path: str, content: str) -> None:
    ensure_dir(os.path.dirname(path))
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)


def main() -> int:
    parser = argparse.ArgumentParser(description="SMI -> PDB -> PDBQT -> Vina workflow")
    parser.add_argument("--smi", required=True, help="Path to input .smi file")
    parser.add_argument("--config", default=CONFIG_FILE_NAME, help="Path to config.txt")
    parser.add_argument("--workdir", default="from_smi", help="Working directory under current dir")
    parser.add_argument("--limit", type=int, default=0, help="Max molecules to process (0=no limit)")
    parser.add_argument("--cpu", type=int, default=0, help="CPUs for Vina (0 autodetect)")
    parser.add_argument("--exhaustiveness", type=int, default=8)
    parser.add_argument("--num-modes", type=int, default=9)
    args = parser.parse_args()

    start_total = time.perf_counter()

    # Resolve paths
    smi_path = os.path.abspath(args.smi)
    if not os.path.exists(smi_path):
        raise FileNotFoundError(f"SMI file not found: {smi_path}")
    config_path = os.path.abspath(args.config)
    receptor_path, box = load_box_from_config(config_path)

    cwd = os.path.dirname(config_path)  # anchor under thousand-dock
    workdir = os.path.join(cwd, args.workdir)
    pdb_dir = os.path.join(workdir, "pdb")
    pdbqt_dir = os.path.join(workdir, "pdbqt")
    out_dir = os.path.join(workdir, "out_vina")
    logs_dir = os.path.join(workdir, "logs")
    ensure_dir(workdir)
    ensure_dir(logs_dir)

    # Step 1: split SMI -> PDB
    t0 = time.perf_counter()
    ret, out, err = run_obabel_split_smi_to_pdb(smi_path, pdb_dir)
    t1 = time.perf_counter()
    write_text(os.path.join(logs_dir, "obabel_split.log"), out + ("\n--- STDERR ---\n" + err if err else ""))
    if ret != 0:
        print("Open Babel split failed. See obabel_split.log")
        return 2

    pdb_files = discover_files([os.path.join(pdb_dir, "*.pdb")])
    if not pdb_files:
        print("No PDB files produced from SMI.")
        return 3

    # Apply limit early if requested
    if args.limit and len(pdb_files) > args.limit:
        pdb_files = pdb_files[: args.limit]

    # Step 2: PDB -> PDBQT
    t2 = time.perf_counter()
    pdbqt_files: List[str] = []
    conv_failures = 0
    for pdb_path in pdb_files:
        base = os.path.splitext(os.path.basename(pdb_path))[0]
        pdbqt_path = os.path.join(pdbqt_dir, f"{base}.pdbqt")
        r, so, se = obabel_pdb_to_pdbqt(pdb_path, pdbqt_path)
        write_text(os.path.join(logs_dir, f"obabel_{base}.log"), so + ("\n--- STDERR ---\n" + se if se else ""))
        if r == 0 and os.path.exists(pdbqt_path) and os.path.getsize(pdbqt_path) > 0:
            pdbqt_files.append(pdbqt_path)
        else:
            conv_failures += 1
    t3 = time.perf_counter()

    if not pdbqt_files:
        print("No PDBQT files produced.")
        return 4

    # Step 3: Dock with Vina
    t4 = time.perf_counter()
    dock_failures = 0
    for ligand_pdbqt in pdbqt_files:
        code, out_pdbqt, vout, verr = vina_dock(
            receptor_path=receptor_path,
            ligand_pdbqt_path=ligand_pdbqt,
            box=box,
            out_dir=out_dir,
            exhaustiveness=args.exhaustiveness,
            num_modes=args.num_modes,
            cpu=args.cpu,
        )
        base = os.path.splitext(os.path.basename(ligand_pdbqt))[0]
        write_text(os.path.join(logs_dir, f"vina_{base}.log"), vout + ("\n--- STDERR ---\n" + verr if verr else ""))
        if code != 0:
            dock_failures += 1
    t5 = time.perf_counter()

    # Summary
    total_s = time.perf_counter() - start_total
    summary = []
    summary.append(f"SMI: {smi_path}")
    summary.append(f"Receptor: {receptor_path}")
    summary.append(f"PDB generated: {len(pdb_files)}")
    summary.append(f"PDBQT generated: {len(pdbqt_files)} (convert failures: {conv_failures})")
    summary.append(f"Dock failures: {dock_failures}")
    summary.append("")
    summary.append(f"Time Step1 (SMI->PDB): {t1 - t0:.2f}s")
    summary.append(f"Time Step2 (PDB->PDBQT): {t3 - t2:.2f}s")
    summary.append(f"Time Step3 (Dock): {t5 - t4:.2f}s")
    summary.append(f"Total time: {total_s:.2f}s")
    write_text(os.path.join(workdir, "summary.txt"), "\n".join(summary) + "\n")

    print("\n".join(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


