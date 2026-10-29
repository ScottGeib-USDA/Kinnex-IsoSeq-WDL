#!/usr/bin/env python3
"""Prepare preprocessing inputs and generate an Atlas Slurm submission script."""

from __future__ import annotations

import argparse
import shlex
import subprocess
import sys
from pathlib import Path

DEFAULT_REPO_URL = "https://github.com/ScottGeib-USDA/Kinnex-IsoSeq-WDL.git"
DEFAULT_CONDA_ENV = "/home/scott.geib/.conda/envs/isoseq_wdl"
DEFAULT_MINIWDL_CONFIG = "~/.config/miniwdl.cfg"
BUILDER_RELATIVE_PATH = Path("scripts/build_inputs_from_files.py")
WORKFLOW_RELATIVE_PATH = Path("workflows/preprocessing.wdl")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Find a HiFi BAM and biosample CSV in a run folder, clone or reuse "
            "Kinnex-IsoSeq-WDL, generate preprocessing inputs JSON, and write "
            "an Atlas Slurm script."
        )
    )
    parser.add_argument(
        "--run-dir",
        type=Path,
        default=Path.cwd(),
        help="Folder containing the HiFi BAM symlink(s) and biosample CSV (default: current directory).",
    )
    parser.add_argument(
        "--repo-url",
        default=DEFAULT_REPO_URL,
        help=f"Git repository URL to clone when needed (default: {DEFAULT_REPO_URL}).",
    )
    parser.add_argument(
        "--branch",
        default="main",
        help="Repository branch to clone (default: main).",
    )
    parser.add_argument(
        "--repo-dir",
        type=Path,
        help="Repository checkout path (default: <run-dir>/Kinnex-IsoSeq-WDL).",
    )
    parser.add_argument(
        "--biosample-csv",
        type=Path,
        help="Explicit biosample CSV path; otherwise detect one in --run-dir.",
    )
    parser.add_argument(
        "--inputs-json",
        type=Path,
        help="Generated inputs JSON path (default: <run-dir>/preprocessing.auto.inputs.json).",
    )
    parser.add_argument(
        "--slurm-script",
        type=Path,
        help="Generated Slurm script path (default: <run-dir>/run_preprocessing.slurm).",
    )
    parser.add_argument(
        "--miniwdl-config",
        default=DEFAULT_MINIWDL_CONFIG,
        help=f"miniwdl config path used by the batch job (default: {DEFAULT_MINIWDL_CONFIG}).",
    )
    parser.add_argument(
        "--conda-env",
        default=DEFAULT_CONDA_ENV,
        help=f"Conda environment activated by the batch job (default: {DEFAULT_CONDA_ENV}).",
    )
    parser.add_argument(
        "--segmentation-adapter-set",
        choices=["8-fold", "12-fold", "16-fold"],
        help="Optional packaged Kinnex segmentation adapter set.",
    )
    parser.add_argument(
        "--hifi-demux-barcodes",
        type=Path,
        help="Optional custom HiFi-demultiplexing barcode FASTA.",
    )
    parser.add_argument(
        "--indexed-primers",
        type=Path,
        help="Optional custom Iso-Seq indexed-primer FASTA.",
    )
    parser.add_argument(
        "--segmentation-adapters",
        type=Path,
        help="Optional custom Skera segmentation-adapter FASTA.",
    )
    return parser.parse_args()


def discover_biosample_csv(run_dir: Path, explicit_path: Path | None) -> Path:
    if explicit_path is not None:
        csv_path = explicit_path.resolve()
        if not csv_path.is_file():
            raise ValueError(f"biosample CSV does not exist or is not a file: {csv_path}")
        return csv_path

    csv_candidates = sorted(path.resolve() for path in run_dir.glob("*.csv") if path.is_file())
    named_candidates = [path for path in csv_candidates if "biosample" in path.name.lower()]

    if len(named_candidates) == 1:
        return named_candidates[0]
    if len(named_candidates) > 1:
        candidates = named_candidates
    elif len(csv_candidates) == 1:
        return csv_candidates[0]
    else:
        candidates = csv_candidates

    if not candidates:
        raise ValueError(f"no CSV files found directly in run directory: {run_dir}")

    listing = "\n".join(f"  - {path}" for path in candidates)
    raise ValueError(
        "could not identify one biosample CSV; pass --biosample-csv explicitly. "
        f"Candidates:\n{listing}"
    )


def find_hifi_bams(run_dir: Path) -> list[Path]:
    bams = sorted(
        path.resolve()
        for path in run_dir.glob("*.bam")
        if path.is_file() and ".hifi_reads." in path.name
    )
    if not bams:
        raise ValueError(
            f"no HiFi BAM files found directly in {run_dir}; expected filenames containing '.hifi_reads.'"
        )
    return bams


def ensure_repo(repo_dir: Path, repo_url: str, branch: str) -> Path:
    repo_dir = repo_dir.resolve()
    if repo_dir.is_dir() and (repo_dir / ".git").exists():
        checkout = repo_dir
    elif repo_dir.exists():
        raise ValueError(
            f"repository destination exists but is not a Git checkout: {repo_dir}; "
            "choose another path with --repo-dir or move that directory"
        )
    else:
        repo_dir.parent.mkdir(parents=True, exist_ok=True)
        print(f"Cloning {repo_url} ({branch}) into {repo_dir}")
        subprocess.run(
            ["git", "clone", "--branch", branch, "--depth", "1", repo_url, str(repo_dir)],
            check=True,
        )
        checkout = repo_dir

    builder = checkout / BUILDER_RELATIVE_PATH
    workflow = checkout / WORKFLOW_RELATIVE_PATH
    if not builder.is_file():
        raise ValueError(
            f"repository checkout does not contain {BUILDER_RELATIVE_PATH}: {checkout}; "
            "make sure the selected branch includes the inputs-builder script"
        )
    if not workflow.is_file():
        raise ValueError(f"preprocessing WDL not found: {workflow}")
    return checkout


def shell_quote(path: Path | str) -> str:
    return shlex.quote(str(path))


def write_slurm_script(
    script_path: Path,
    run_dir: Path,
    repo_dir: Path,
    inputs_json: Path,
    miniwdl_config: str,
    conda_env: str,
) -> None:
    output_path = script_path.resolve()
    config_path = Path(miniwdl_config).expanduser().resolve()
    workflow_path = repo_dir / WORKFLOW_RELATIVE_PATH
    content = f'''#!/bin/bash

#SBATCH -p atlas 
#SBATCH --time=40:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=48   # 20 processor core(s) per node X 2 threads per core
#SBATCH --partition=atlas    # standard node(s)
#SBATCH --job-name="isoseq"
#SBATCH --mail-user=scott.geib@usda.gov   # email address
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output="isoseq-%j-%N.out" # job standard output file (%j replaced by job id)
#SBATCH --error="isoseq-%j-%N.err" # job standard error file (%j replaced by job id)
#SBATCH --account=ag100pest

set -euo pipefail

cd {shell_quote(run_dir)}
module load apptainer
apptainer --version
source activate {shell_quote(conda_env)}
miniwdl run --cfg {shell_quote(config_path)} \\
  -i {shell_quote(inputs_json)} \\
  {shell_quote(workflow_path)}
'''
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(content, encoding="utf-8")
    output_path.chmod(0o755)


def main() -> int:
    args = parse_args()
    run_dir = args.run_dir.resolve()
    if not run_dir.is_dir():
        raise ValueError(f"run directory does not exist or is not a directory: {run_dir}")

    hifi_bams = find_hifi_bams(run_dir)
    biosample_csv = discover_biosample_csv(run_dir, args.biosample_csv)
    repo_dir = args.repo_dir or run_dir / "Kinnex-IsoSeq-WDL"
    repo_dir = ensure_repo(repo_dir, args.repo_url, args.branch)

    inputs_json = args.inputs_json or run_dir / "preprocessing.auto.inputs.json"
    inputs_json = inputs_json.resolve()
    builder_command = [
        sys.executable,
        str(repo_dir / BUILDER_RELATIVE_PATH),
        "preprocessing",
        "--repo-dir",
        str(repo_dir),
        "--run-dir",
        str(run_dir),
        "--biosample-csv",
        str(biosample_csv),
        "--output-json",
        str(inputs_json),
    ]
    if args.segmentation_adapter_set is not None:
        builder_command.extend(
            ["--segmentation-adapter-set", args.segmentation_adapter_set]
        )
    for option, resource_path in (
        ("--hifi-demux-barcodes", args.hifi_demux_barcodes),
        ("--indexed-primers", args.indexed_primers),
        ("--segmentation-adapters", args.segmentation_adapters),
    ):
        if resource_path is not None:
            builder_command.extend([option, str(resource_path.resolve())])

    print(f"HiFi BAM files detected: {len(hifi_bams)}")
    print(f"Biosample CSV: {biosample_csv}")
    subprocess.run(builder_command, check=True)

    slurm_script = args.slurm_script or run_dir / "run_preprocessing.slurm"
    slurm_script = slurm_script.resolve()
    write_slurm_script(
        script_path=slurm_script,
        run_dir=run_dir,
        repo_dir=repo_dir,
        inputs_json=inputs_json,
        miniwdl_config=args.miniwdl_config,
        conda_env=args.conda_env,
    )

    print(f"Slurm script: {slurm_script}")
    print("Submit it from the run directory with: sbatch " + shell_quote(slurm_script))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
