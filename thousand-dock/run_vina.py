#!/usr/bin/env python3
"""
Batch runner for AutoDock Vina using parameters in config.txt.

- Reads receptor and grid box (center/size) from config.txt
- Docks all ligands matching provided glob patterns
- Writes output PDBQT and saves Vina stdout/stderr to per-ligand log files

Usage examples:
  conda run -n vina python run_vina.py --limit 1 --exhaustiveness 4 --num-modes 5
  conda run -n vina python run_vina.py --pattern "prepared/*.pdbqt" --pattern "MOL_*.pdbqt"
"""

from __future__ import annotations

import argparse
import glob
import os
import shlex
import subprocess
import sys
from dataclasses import dataclass
from typing import Dict, Iterable, List


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


def load_box_from_config(config_path: str) -> DockBox:
    kv = parse_simple_kv_config(config_path)
    required_keys = [
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
    return DockBox(
        center_x=float(kv["center_x"]),
        center_y=float(kv["center_y"]),
        center_z=float(kv["center_z"]),
        size_x=float(kv["size_x"]),
        size_y=float(kv["size_y"]),
        size_z=float(kv["size_z"]),
    )


def discover_ligands(patterns: Iterable[str]) -> List[str]:
    ligands: List[str] = []
    for pattern in patterns:
        ligands.extend(glob.glob(pattern))
    # De-duplicate while preserving order
    seen: set[str] = set()
    unique_ligands: List[str] = []
    for path in ligands:
        if path not in seen and os.path.isfile(path):
            seen.add(path)
            unique_ligands.append(path)
    return unique_ligands


def ensure_directory(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def run_vina(
    receptor_path: str,
    ligand_path: str,
    box: DockBox,
    out_dir: str,
    exhaustiveness: int,
    num_modes: int,
    cpu: int | None,
) -> int:
    base_name = os.path.splitext(os.path.basename(ligand_path))[0]
    out_pdbqt = os.path.join(out_dir, f"{base_name}_out.pdbqt")
    log_path = os.path.join(out_dir, f"{base_name}.log")

    cmd: List[str] = [
        "vina",
        "--receptor",
        receptor_path,
        "--ligand",
        ligand_path,
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
    if cpu is not None and cpu > 0:
        cmd.extend(["--cpu", str(cpu)])

    print(f"Running: {' '.join(shlex.quote(c) for c in cmd)}")
    sys.stdout.flush()

    completed = subprocess.run(cmd, capture_output=True, text=True)
    with open(log_path, "w", encoding="utf-8") as log_file:
        log_file.write(completed.stdout)
        if completed.stderr:
            log_file.write("\n--- STDERR ---\n")
            log_file.write(completed.stderr)

    if completed.returncode != 0:
        print(f"FAILED: {ligand_path} (exit={completed.returncode}) -> see {log_path}")
    else:
        print(f"DONE: {ligand_path} -> {out_pdbqt}")
    return completed.returncode


def main() -> int:
    parser = argparse.ArgumentParser(description="Batch run AutoDock Vina using config.txt")
    parser.add_argument(
        "--config",
        default=CONFIG_FILE_NAME,
        help="Path to config.txt containing receptor and box parameters",
    )
    parser.add_argument(
        "--pattern",
        action="append",
        default=["prepared/*.pdbqt", "MOL_*.pdbqt"],
        help="Glob pattern(s) for ligand PDBQT files (can be specified multiple times)",
    )
    parser.add_argument(
        "--outdir",
        default="out_vina",
        help="Directory to write outputs",
    )
    parser.add_argument("--exhaustiveness", type=int, default=8)
    parser.add_argument("--num-modes", type=int, default=9)
    parser.add_argument("--cpu", type=int, default=0, help="Number of CPUs to use (0 autodetect)")
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Limit number of ligands processed (0 means no limit)",
    )

    args = parser.parse_args()

    config_path = os.path.abspath(args.config)
    kv = parse_simple_kv_config(config_path)
    receptor = kv.get("receptor")
    if not receptor:
        raise ValueError(f"Missing 'receptor' in {config_path}")
    receptor_path = os.path.abspath(receptor)
    if not os.path.exists(receptor_path):
        # Try relative to config directory
        receptor_path = os.path.join(os.path.dirname(config_path), receptor)
        receptor_path = os.path.abspath(receptor_path)
    if not os.path.exists(receptor_path):
        raise FileNotFoundError(f"Receptor not found: {receptor_path}")

    dock_box = load_box_from_config(config_path)
    ligand_patterns = [os.path.join(os.path.dirname(config_path), p) if not os.path.isabs(p) else p for p in args.pattern]
    ligands = discover_ligands(ligand_patterns)

    if not ligands:
        print("No ligands matched the given patterns.")
        return 1

    ensure_directory(args.outdir)
    processed = 0
    failures = 0
    for ligand_path in ligands:
        exit_code = run_vina(
            receptor_path=receptor_path,
            ligand_path=ligand_path,
            box=dock_box,
            out_dir=args.outdir,
            exhaustiveness=args.exhaustiveness,
            num_modes=args.num_modes,
            cpu=(args.cpu if args.cpu > 0 else None),
        )
        processed += 1
        if exit_code != 0:
            failures += 1
        if args.limit and processed >= args.limit:
            break

    print(f"Processed ligands: {processed}, failures: {failures}")
    return 0 if failures == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())


