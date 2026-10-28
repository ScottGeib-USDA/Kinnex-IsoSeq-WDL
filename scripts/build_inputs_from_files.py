#!/usr/bin/env python3
"""Build workflow input JSON files from SMRT Link run outputs.

This helper fills the run-specific hifi_sources and biosample_csv inputs for the
public HPC templates by inspecting a SMRT Link run directory.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path


TEMPLATE_BY_WORKFLOW_AND_MODE = {
    ("kinnex_isoseq", "no_instrument_demux"): Path(
        "backends/hpc/kinnex_isoseq.no_instrument_demux.hpc.inputs.json"
    ),
    ("kinnex_isoseq", "from_instrument_demux"): Path(
        "backends/hpc/kinnex_isoseq.from_instrument_demux.hpc.inputs.json"
    ),
    ("preprocessing", "no_instrument_demux"): Path("backends/hpc/preprocessing.hpc.inputs.json"),
    ("preprocessing", "from_instrument_demux"): Path("backends/hpc/preprocessing.hpc.inputs.json"),
}

RAW_HIFI_SUFFIX = ".hifi_reads.bam"
HIFI_READS_MARKER = ".hifi_reads."
DEMUX_HIFI_PATTERN = re.compile(r"\.hifi_reads\.([^.]+)\.bam$")


@dataclass(frozen=True)
class DemuxEntry:
    bam: Path
    hifi_barcode: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a workflow inputs JSON from SMRT Link hifi_reads outputs. "
            "The tool auto-detects run mode and fills hifi_sources."
        )
    )
    parser.add_argument(
        "--repo-dir",
        type=Path,
        help=(
            "Path to the Kinnex-IsoSeq-WDL repository root that contains "
            "backends/hpc. When omitted, the tool checks the current working "
            "directory, then ./Kinnex-IsoSeq_WDL, then ./Kinnex-IsoSeq-WDL."
        ),
    )
    parser.add_argument(
        "--workflow",
        choices=["kinnex_isoseq", "preprocessing"],
        default="kinnex_isoseq",
        help="Workflow namespace to populate in the output JSON.",
    )
    parser.add_argument(
        "--run-dir",
        type=Path,
        default=Path.cwd(),
        help=(
            "SMRT Link run output directory that contains hifi_reads BAM files. "
            "Defaults to the current working directory."
        ),
    )
    parser.add_argument(
        "--biosample-csv",
        type=Path,
        help=(
            "Path to the 3-column biosample CSV used by the workflow. "
            "When omitted, the tool auto-detects one CSV under --run-dir."
        ),
    )
    parser.add_argument(
        "--mode",
        choices=["auto", "no_instrument_demux", "from_instrument_demux"],
        default="auto",
        help="Run mode. Auto infers mode from discovered hifi_reads BAM names.",
    )
    parser.add_argument(
        "--template-json",
        type=Path,
        help="Optional explicit template JSON path. Overrides built-in template selection.",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        help=(
            "Path to write the generated workflow inputs JSON. "
            "Defaults to <run-dir>/<workflow>.auto.inputs.json."
        ),
    )
    parser.add_argument(
        "--output-prefix",
        help=(
            "Optional output prefix for kinnex_isoseq.output_prefix. "
            "Defaults to the run directory name."
        ),
    )
    parser.add_argument(
        "--segmentation-adapter-set",
        choices=["8-fold", "12-fold", "16-fold"],
        help=(
            "Optional segmentation adapter set to write into the output JSON. "
            "If omitted, the template value is preserved."
        ),
    )
    parser.add_argument(
        "--path-prefix-rewrite",
        nargs=2,
        metavar=("OLD_PREFIX", "NEW_PREFIX"),
        help=(
            "Optional absolute-prefix rewrite applied to discovered BAM paths and "
            "the biosample CSV path before writing JSON."
        ),
    )
    return parser.parse_args()


def discover_hifi_bams(run_dir: Path) -> tuple[list[Path], list[DemuxEntry]]:
    if not run_dir.exists():
        raise ValueError(f"run directory does not exist: {run_dir}")
    if not run_dir.is_dir():
        raise ValueError(f"run directory is not a directory: {run_dir}")

    raw_candidates: list[Path] = []
    demux_candidates: list[DemuxEntry] = []
    for bam in sorted(run_dir.rglob("*.bam")):
        name = bam.name
        demux_match = DEMUX_HIFI_PATTERN.search(name)
        if demux_match:
            token = demux_match.group(1)
            if token.startswith("bc"):
                barcode = token if "--" in token else f"{token}--{token}"
                demux_candidates.append(DemuxEntry(bam=bam.resolve(), hifi_barcode=barcode))
                continue

        if name.endswith(RAW_HIFI_SUFFIX):
            raw_candidates.append(bam.resolve())
            continue

        # Ignore unrelated BAMs; only names with the hifi_reads marker are in scope.
        if HIFI_READS_MARKER not in name:
            continue

    return raw_candidates, demux_candidates


def discover_biosample_csv(run_dir: Path) -> Path:
    candidates = sorted(
        {
            path.resolve()
            for path in run_dir.glob("*.csv")
            if "biosample" in path.name.lower()
        }
    )
    if len(candidates) == 0:
        raise ValueError(
            "could not auto-detect biosample CSV in --run-dir; "
            "provide --biosample-csv"
        )
    if len(candidates) > 1:
        found = "\n".join(f"  - {path}" for path in candidates)
        raise ValueError(
            "found multiple biosample CSV candidates in --run-dir; "
            "provide --biosample-csv explicitly:\n"
            f"{found}"
        )
    return candidates[0]


def resolve_repo_root(args: argparse.Namespace) -> Path:
    if args.repo_dir is not None:
        candidates = [args.repo_dir.resolve()]
    else:
        cwd = Path.cwd().resolve()
        candidates = [
            cwd,
            cwd / "Kinnex-IsoSeq_WDL",
            cwd / "Kinnex-IsoSeq-WDL",
        ]

    for candidate in candidates:
        if (candidate / "backends" / "hpc").is_dir():
            return candidate

    searched = "\n".join(f"  - {candidate}" for candidate in candidates)
    raise ValueError(
        "could not find repository root with backends/hpc. "
        "Set --repo-dir explicitly. Checked:\n"
        f"{searched}"
    )


def infer_mode(
    requested_mode: str,
    raw_candidates: list[Path],
    demux_candidates: list[DemuxEntry],
) -> str:
    if requested_mode != "auto":
        return requested_mode

    has_raw = len(raw_candidates) > 0
    has_demux = len(demux_candidates) > 0

    if has_raw and not has_demux:
        return "no_instrument_demux"
    if has_demux and not has_raw:
        return "from_instrument_demux"
    if not has_raw and not has_demux:
        raise ValueError("no hifi_reads BAM files were found under --run-dir")

    raise ValueError(
        "both raw and demultiplexed hifi_reads BAM patterns were detected; "
        "set --mode explicitly"
    )


def rewrite_path(path: Path, rewrite: tuple[Path, Path] | None) -> str:
    resolved = path.resolve()
    if rewrite is None:
        return str(resolved)

    old_prefix, new_prefix = rewrite
    try:
        suffix = resolved.relative_to(old_prefix)
    except ValueError as exc:
        raise ValueError(
            f"path {resolved} is not under rewrite prefix {old_prefix}; "
            "adjust --path-prefix-rewrite or remove it"
        ) from exc

    return str((new_prefix / suffix).resolve())


def resolve_template_path(repo_root: Path, args: argparse.Namespace, mode: str) -> Path:
    if args.template_json is not None:
        template_path = args.template_json
    else:
        template_rel = TEMPLATE_BY_WORKFLOW_AND_MODE[(args.workflow, mode)]
        template_path = repo_root / template_rel

    if not template_path.exists():
        raise ValueError(f"template JSON not found: {template_path}")
    return template_path.resolve()


def build_hifi_sources(
    mode: str,
    raw_candidates: list[Path],
    demux_candidates: list[DemuxEntry],
    rewrite: tuple[Path, Path] | None,
) -> list[dict[str, str]]:
    if mode == "no_instrument_demux":
        if len(raw_candidates) != 1:
            raise ValueError(
                "no_instrument_demux mode requires exactly one raw hifi_reads BAM; "
                f"found {len(raw_candidates)}"
            )
        return [{"hifi_bam": rewrite_path(raw_candidates[0], rewrite)}]

    if mode == "from_instrument_demux":
        if len(demux_candidates) == 0:
            raise ValueError(
                "from_instrument_demux mode requires demultiplexed hifi_reads BAMs; none were found"
            )
        return [
            {
                "hifi_bam": rewrite_path(entry.bam, rewrite),
                "hifi_barcode": entry.hifi_barcode,
            }
            for entry in demux_candidates
        ]

    raise ValueError(f"unsupported mode: {mode}")


def main() -> int:
    args = parse_args()
    repo_root = resolve_repo_root(args)
    run_dir = args.run_dir.resolve()

    raw_candidates, demux_candidates = discover_hifi_bams(run_dir)
    mode = infer_mode(args.mode, raw_candidates, demux_candidates)

    biosample_csv = (
        args.biosample_csv.resolve()
        if args.biosample_csv is not None
        else discover_biosample_csv(run_dir)
    )

    rewrite: tuple[Path, Path] | None = None
    if args.path_prefix_rewrite is not None:
        old_prefix = Path(args.path_prefix_rewrite[0]).resolve()
        new_prefix = Path(args.path_prefix_rewrite[1]).resolve()
        rewrite = (old_prefix, new_prefix)

    template_path = resolve_template_path(repo_root, args, mode)
    with template_path.open("rt", encoding="utf-8") as handle:
        output_json = json.load(handle)

    ns = args.workflow
    output_json[f"{ns}.hifi_sources"] = build_hifi_sources(
        mode=mode,
        raw_candidates=raw_candidates,
        demux_candidates=demux_candidates,
        rewrite=rewrite,
    )
    output_json[f"{ns}.biosample_csv"] = rewrite_path(biosample_csv, rewrite)

    if args.segmentation_adapter_set is not None:
        output_json[f"{ns}.segmentation_adapter_set"] = args.segmentation_adapter_set

    if ns == "kinnex_isoseq":
        output_json[f"{ns}.output_prefix"] = args.output_prefix or run_dir.name

    output_json_path = args.output_json or (run_dir / f"{ns}.auto.inputs.json")

    output_json_path.parent.mkdir(parents=True, exist_ok=True)
    with output_json_path.open("wt", encoding="utf-8") as handle:
        json.dump(output_json, handle, indent=2)
        handle.write("\n")

    print(f"Mode: {mode}")
    print(f"Template: {template_path}")
    print(f"Run directory: {run_dir}")
    print(f"Biosample CSV: {biosample_csv}")
    print(f"Wrote: {output_json_path.resolve()}")
    print(f"HiFi sources: {len(output_json[f'{ns}.hifi_sources'])}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)