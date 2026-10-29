#!/usr/bin/env python3
"""Copy final FLNC BAMs and indexes into a flat, sample-named results directory."""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path


def resolve_manifest(source: Path) -> Path:
    if source.is_dir():
        if (source / "_LAST").exists():
            source = source / "_LAST" / "outputs.json"
        else:
            source = source / "outputs.json"
    if not source.is_file():
        raise ValueError(f"miniwdl output manifest not found: {source}")
    return source.resolve()


def collect_results(manifest: Path, results_dir: Path) -> None:
    manifest = resolve_manifest(manifest)
    samtools = shutil.which("samtools")
    if samtools is None:
        raise ValueError("samtools must be available in PATH to read FLNC sample names")
    with manifest.open(encoding="utf-8") as handle:
        document = json.load(handle)
    outputs = document.get("outputs", document)
    bams = outputs.get("preprocessing.flnc_bams")
    indexes = outputs.get("preprocessing.flnc_bam_pbis")
    if not isinstance(bams, list) or not bams:
        raise ValueError("manifest must contain non-empty preprocessing.flnc_bams")
    if not isinstance(indexes, list) or len(indexes) != len(bams):
        raise ValueError("manifest must contain one preprocessing.flnc_bam_pbis entry per BAM")

    planned: list[tuple[Path, Path]] = []
    destinations: set[Path] = set()
    for bam_value, index_value in zip(bams, indexes):
        bam = Path(bam_value)
        index = Path(index_value)
        bam = bam if bam.is_absolute() else manifest.parent / bam
        index = index if index.is_absolute() else manifest.parent / index
        if not bam.is_file() or not index.is_file():
            raise ValueError(f"missing output BAM or PBI: {bam}, {index}")
        header = subprocess.run(
            [samtools, "view", "-H", str(bam)],
            check=True, capture_output=True, text=True,
        ).stdout
        read_groups = [line.split("\t")[1:] for line in header.splitlines() if line.startswith("@RG\t")]
        samples: set[str] = set()
        if not read_groups:
            raise ValueError(f"{bam}: no @RG records")
        for fields in read_groups:
            values = [field[3:] for field in fields if field.startswith("SM:")]
            if len(values) != 1 or not values[0]:
                raise ValueError(f"{bam}: each @RG must have exactly one non-empty SM tag")
            samples.add(values[0])
        if len(samples) != 1:
            raise ValueError(f"{bam}: expected one distinct SM value, found {sorted(samples)}")
        sample = next(iter(samples))
        if not re.fullmatch(r"[A-Za-z0-9_-]{1,40}", sample):
            raise ValueError(f"{bam}: unsafe preprocessing sample name: {sample!r}")
        destination_bam = results_dir / f"{sample}.{bam.name}"
        for source, destination in (
            (bam, destination_bam),
            (index, Path(str(destination_bam) + ".pbi")),
        ):
            if destination in destinations or destination.exists() or destination.is_symlink():
                raise ValueError(f"refusing to overwrite or duplicate result: {destination}")
            destinations.add(destination)
            planned.append((source, destination))

    results_dir.mkdir(parents=True, exist_ok=True)
    for source, destination in planned:
        with source.open("rb") as reader, destination.open("xb") as writer:
            try:
                shutil.copyfileobj(reader, writer)
            except BaseException:
                destination.unlink()
                raise
        print(f"Collected: {destination}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "manifest", type=Path,
        help="Output JSON, execution directory, or run folder containing _LAST/outputs.json",
    )
    parser.add_argument("--results-dir", type=Path, required=True)
    args = parser.parse_args()
    collect_results(args.manifest.resolve(), args.results_dir.resolve())


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)