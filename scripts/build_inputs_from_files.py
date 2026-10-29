#!/usr/bin/env python3
"""Build workflow input JSON files from local run files.

For preprocessing and kinnex_isoseq, this helper fills run-specific hifi_sources
and biosample_csv inputs by inspecting a run directory. For secondary_analysis,
it discovers FLNC BAMs.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import subprocess
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
    ("secondary_analysis", "flnc"): Path("backends/hpc/secondary_analysis.hpc.inputs.json"),
}

RAW_HIFI_SUFFIX = ".hifi_reads.bam"
HIFI_READS_MARKER = ".hifi_reads."
DEMUX_HIFI_PATTERN = re.compile(r"\.hifi_reads\.([^.]+)\.bam$")
BIOSAMPLE_ALLOWED_SAMPLE_NAME = re.compile(r"^[A-Za-z0-9_-]+$")


@dataclass(frozen=True)
class DemuxEntry:
    bam: Path
    hifi_barcode: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a workflow inputs JSON from local BAM/CSV files. "
            "For preprocessing/kinnex modes it auto-detects hifi_sources; "
            "for secondary_analysis it discovers FLNC BAMs."
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
        "workflow",
        choices=["kinnex_isoseq", "preprocessing", "secondary_analysis"],
        help="Workflow namespace to populate in the output JSON.",
    )
    parser.add_argument(
        "--run-dir",
        type=Path,
        default=Path.cwd(),
        help=(
            "Run directory to scan for inputs. For preprocessing/kinnex it scans "
            "for hifi_reads BAMs; for secondary_analysis it scans for *.flnc.bam. "
            "Defaults to the current working directory."
        ),
    )
    parser.add_argument(
        "--biosample-csv",
        type=Path,
        help=(
            "Path to the 3-column biosample CSV used by the workflow. "
            "When omitted, the tool auto-detects one CSV in --run-dir. "
            "Not used for secondary_analysis."
        ),
    )
    parser.add_argument(
        "--updated-biosample-csv",
        type=Path,
        help=(
            "Output path for a biosample CSV updated from a single detected HiFi "
            "barcode. Defaults to <biosample-name>.updated.csv beside the input."
        ),
    )
    parser.add_argument(
        "--mode",
        choices=["auto", "no_instrument_demux", "from_instrument_demux"],
        default="auto",
        help=(
            "Run mode for preprocessing/kinnex workflows. "
            "Auto infers mode from discovered hifi_reads BAM names. "
            "Not used for secondary_analysis."
        ),
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
    parser.add_argument(
        "--genome-fasta",
        type=Path,
        help="Optional custom uncompressed reference genome FASTA.",
    )
    parser.add_argument(
        "--annotation-gtf-gz",
        type=Path,
        help="Optional matching compressed annotation GTF for a custom genome.",
    )
    parser.add_argument(
        "--path-prefix-rewrite",
        nargs=2,
        metavar=("OLD_PREFIX", "NEW_PREFIX"),
        help=(
            "Optional absolute-prefix rewrite applied to discovered BAM paths and "
            "the biosample CSV and reference override paths before writing JSON."
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


def discover_flnc_bams(run_dir: Path) -> list[Path]:
    flnc_candidates = sorted(
        {
            bam.resolve()
            for bam in run_dir.rglob("*.bam")
            if bam.name.endswith(".flnc.bam")
        }
    )
    if len(flnc_candidates) == 0:
        raise ValueError(
            "no FLNC BAM files were found under --run-dir; expected one or more *.flnc.bam files"
        )
    return flnc_candidates


def validate_flnc_sm_tags_with_samtools(flnc_bams: list[Path]) -> list[str]:
    warnings: list[str] = []
    samtools = shutil.which("samtools")
    if samtools is None:
        warnings.append(
            "samtools was not found in PATH; skipping optional FLNC SM-tag validation"
        )
        return warnings

    for bam in flnc_bams:
        result = subprocess.run(
            [samtools, "view", "-H", str(bam)],
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode != 0:
            message = result.stderr.strip() or result.stdout.strip() or "unknown samtools error"
            warnings.append(f"{bam}: could not read BAM header with samtools ({message})")
            continue

        rg_lines = [line for line in result.stdout.splitlines() if line.startswith("@RG\t")]
        if not rg_lines:
            warnings.append(f"{bam}: no @RG records found; secondary_analysis expects SM tags")
            continue

        sm_values: list[str] = []
        missing_sm = 0
        for rg_line in rg_lines:
            sm_value = None
            for field in rg_line.split("\t")[1:]:
                if field.startswith("SM:"):
                    sm_value = field[3:]
                    break
            if sm_value is None or sm_value.strip() == "":
                missing_sm += 1
                continue
            sm_values.append(sm_value)

        if missing_sm > 0:
            warnings.append(f"{bam}: {missing_sm} @RG record(s) missing non-empty SM tag")

        unique_sm_values = set(sm_values)
        if len(unique_sm_values) > 1:
            preview = ", ".join(sorted(unique_sm_values)[:3])
            warnings.append(
                f"{bam}: multiple SM values detected in @RG ({preview}); "
                "secondary_analysis expects one SM per BAM"
            )

    return warnings


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


def is_valid_barcode_pair(value: str) -> bool:
    return value.count("--") == 1 and all(part != "" for part in value.split("--"))


def validate_biosample_csv(path: Path) -> None:
    expected_header = ["HiFi Barcode", "cDNA Barcode", "Bio Sample"]
    errors: list[str] = []
    seen_cdna_per_hifi: dict[str, set[str]] = {}

    with path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader, None)
        if header is None:
            raise ValueError(f"biosample CSV is empty: {path}")
        if [col.strip() for col in header] != expected_header:
            raise ValueError(
                "biosample CSV header must be exactly: "
                "HiFi Barcode,cDNA Barcode,Bio Sample"
            )

        for row_number, row in enumerate(reader, start=2):
            # Empty rows are allowed and ignored.
            if len(row) == 0 or all(cell.strip() == "" for cell in row):
                continue

            if len(row) != 3:
                errors.append(f"line {row_number}: expected 3 columns, found {len(row)}")
                continue

            hifi_barcode, cdna_barcode, bio_sample = (cell.strip() for cell in row)

            if hifi_barcode == "" or cdna_barcode == "" or bio_sample == "":
                errors.append(f"line {row_number}: each row must have three non-empty fields")
                continue

            if not is_valid_barcode_pair(hifi_barcode):
                errors.append(
                    f"line {row_number}: HiFi Barcode must be an exact two-part pair separated by '--'"
                )

            if not is_valid_barcode_pair(cdna_barcode):
                errors.append(
                    f"line {row_number}: cDNA Barcode must be an exact two-part pair separated by '--'"
                )

            if len(bio_sample) > 40:
                errors.append(f"line {row_number}: Bio Sample must be 40 characters or fewer")

            if not BIOSAMPLE_ALLOWED_SAMPLE_NAME.fullmatch(bio_sample):
                errors.append(
                    f"line {row_number}: Bio Sample may contain only letters, numbers, '_' and '-'"
                )

            seen_cdna = seen_cdna_per_hifi.setdefault(hifi_barcode, set())
            if cdna_barcode in seen_cdna:
                errors.append(
                    f"line {row_number}: duplicate cDNA Barcode '{cdna_barcode}' within HiFi Barcode group '{hifi_barcode}'"
                )
            else:
                seen_cdna.add(cdna_barcode)

    if errors:
        raise ValueError("biosample CSV validation failed:\n  - " + "\n  - ".join(errors))


def synchronize_biosample_hifi_barcode(
    csv_path: Path,
    demux_candidates: list[DemuxEntry],
    output_path: Path | None,
) -> Path:
    detected_barcodes = {entry.hifi_barcode for entry in demux_candidates}
    if not detected_barcodes:
        return csv_path

    with csv_path.open("rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle))

    if not rows:
        raise ValueError(f"biosample CSV is empty: {csv_path}")

    existing_barcodes = {
        row[0].strip()
        for row in rows[1:]
        if row and any(cell.strip() for cell in row) and row[0].strip()
    }
    missing_barcodes = detected_barcodes - existing_barcodes
    if not missing_barcodes:
        return csv_path

    if len(detected_barcodes) != 1 or len(existing_barcodes) != 1:
        detected = ", ".join(sorted(detected_barcodes))
        existing = ", ".join(sorted(existing_barcodes)) or "<none>"
        raise ValueError(
            "cannot safely update biosample CSV HiFi Barcode values: BAMs contain "
            f"[{detected}], while the CSV contains [{existing}]. Update the CSV "
            "barcode groups to match the BAM filenames or provide an unambiguous "
            "single-barcode-group CSV."
        )

    detected_barcode = next(iter(detected_barcodes))
    if not is_valid_barcode_pair(detected_barcode):
        raise ValueError(
            f"HiFi barcode parsed from BAM filename is not a barcode pair: {detected_barcode}"
        )

    updated_path = output_path or csv_path.with_name(
        f"{csv_path.stem}.updated{csv_path.suffix}"
    )
    updated_path = updated_path.resolve()
    if updated_path == csv_path.resolve():
        raise ValueError(
            "--updated-biosample-csv must be different from --biosample-csv; "
            "the original sample sheet is preserved"
        )

    updated_rows = [rows[0][:]]
    for row in rows[1:]:
        if not row or all(cell.strip() == "" for cell in row):
            updated_rows.append(row)
            continue
        updated_row = row[:]
        updated_row[0] = detected_barcode
        updated_rows.append(updated_row)

    updated_path.parent.mkdir(parents=True, exist_ok=True)
    with updated_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerows(updated_rows)

    validate_biosample_csv(updated_path)
    print(
        f"Updated biosample CSV HiFi Barcode column from "
        f"{next(iter(existing_barcodes))} to {detected_barcode}: {updated_path}"
    )
    return updated_path


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


def build_reference_overrides(
    args: argparse.Namespace,
    workflow: str,
    rewrite: tuple[Path, Path] | None,
) -> dict[str, str]:
    requested = {
        "hifi_demux_barcodes": args.hifi_demux_barcodes,
        "indexed_primers": args.indexed_primers,
        "segmentation_adapters": args.segmentation_adapters,
        "genome_fasta": args.genome_fasta,
        "annotation_gtf_gz": args.annotation_gtf_gz,
    }
    supplied = {name: path for name, path in requested.items() if path is not None}
    supported = {
        "kinnex_isoseq": set(requested),
        "preprocessing": {
            "hifi_demux_barcodes",
            "indexed_primers",
            "segmentation_adapters",
        },
        "secondary_analysis": {"genome_fasta", "annotation_gtf_gz"},
    }[workflow]
    unsupported = sorted(set(supplied) - supported)
    if unsupported:
        flags = ", ".join(f"--{name.replace('_', '-')}" for name in unsupported)
        raise ValueError(f"{flags} cannot be used with --workflow {workflow}")
    if "genome_fasta" in supplied and "annotation_gtf_gz" not in supplied:
        raise ValueError("--genome-fasta requires --annotation-gtf-gz with the matching annotation")

    overrides: dict[str, str] = {}
    for name, path in supplied.items():
        if not path.is_file():
            raise ValueError(f"--{name.replace('_', '-')} is not a file: {path}")
        overrides[name] = rewrite_path(path, rewrite)
    return overrides


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

    mode = "flnc" if args.workflow == "secondary_analysis" else "auto"
    raw_candidates: list[Path] = []
    demux_candidates: list[DemuxEntry] = []
    flnc_candidates: list[Path] = []
    biosample_csv: Path | None = None

    if args.workflow == "secondary_analysis":
        flnc_candidates = discover_flnc_bams(run_dir)
        if args.mode != "auto":
            raise ValueError("--mode is not used with --workflow secondary_analysis")
        if args.biosample_csv is not None:
            raise ValueError("--biosample-csv is not used with --workflow secondary_analysis")
        for warning in validate_flnc_sm_tags_with_samtools(flnc_candidates):
            print(f"WARNING: {warning}", file=sys.stderr)
    else:
        raw_candidates, demux_candidates = discover_hifi_bams(run_dir)
        mode = infer_mode(args.mode, raw_candidates, demux_candidates)
        biosample_csv = (
            args.biosample_csv.resolve()
            if args.biosample_csv is not None
            else discover_biosample_csv(run_dir)
        )
        biosample_csv = synchronize_biosample_hifi_barcode(
            biosample_csv,
            demux_candidates if mode == "from_instrument_demux" else [],
            args.updated_biosample_csv,
        )
        validate_biosample_csv(biosample_csv)

    rewrite: tuple[Path, Path] | None = None
    if args.path_prefix_rewrite is not None:
        old_prefix = Path(args.path_prefix_rewrite[0]).resolve()
        new_prefix = Path(args.path_prefix_rewrite[1]).resolve()
        rewrite = (old_prefix, new_prefix)

    reference_overrides = build_reference_overrides(args, args.workflow, rewrite)

    template_path = resolve_template_path(repo_root, args, mode)
    with template_path.open("rt", encoding="utf-8") as handle:
        output_json = json.load(handle)

    ns = args.workflow

    if reference_overrides:
        override_key = f"{ns}.reference_overrides"
        template_overrides = output_json.get(override_key, {})
        if not isinstance(template_overrides, dict):
            raise ValueError(f"{override_key} in the template must be a JSON object")
        template_overrides.update(reference_overrides)
        output_json[override_key] = template_overrides
        if "indexed_primers" in reference_overrides:
            output_json.pop(f"{ns}.isoseq_primers_set", None)
        if "segmentation_adapters" in reference_overrides:
            output_json.pop(f"{ns}.segmentation_adapter_set", None)

    if ns == "secondary_analysis":
        output_json[f"{ns}.flnc_bams"] = [
            rewrite_path(flnc_bam, rewrite)
            for flnc_bam in flnc_candidates
        ]
    else:
        output_json[f"{ns}.hifi_sources"] = build_hifi_sources(
            mode=mode,
            raw_candidates=raw_candidates,
            demux_candidates=demux_candidates,
            rewrite=rewrite,
        )
        output_json[f"{ns}.biosample_csv"] = rewrite_path(biosample_csv, rewrite)

    if args.segmentation_adapter_set is not None:
        if ns == "secondary_analysis":
            raise ValueError(
                "--segmentation-adapter-set is only valid for preprocessing or kinnex_isoseq workflows"
            )
        output_json[f"{ns}.segmentation_adapter_set"] = args.segmentation_adapter_set

    if ns in {"kinnex_isoseq", "secondary_analysis"}:
        output_json[f"{ns}.output_prefix"] = args.output_prefix or run_dir.name

    output_json_path = args.output_json or (run_dir / f"{ns}.auto.inputs.json")

    output_json_path.parent.mkdir(parents=True, exist_ok=True)
    with output_json_path.open("wt", encoding="utf-8") as handle:
        json.dump(output_json, handle, indent=2)
        handle.write("\n")

    print(f"Mode: {mode}")
    print(f"Template: {template_path}")
    print(f"Run directory: {run_dir}")
    if biosample_csv is not None:
        print(f"Biosample CSV: {biosample_csv}")
    print(f"Wrote: {output_json_path.resolve()}")
    if ns == "secondary_analysis":
        print(f"FLNC BAMs: {len(output_json[f'{ns}.flnc_bams'])}")
    else:
        print(f"HiFi sources: {len(output_json[f'{ns}.hifi_sources'])}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)