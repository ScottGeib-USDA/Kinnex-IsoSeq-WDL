<h1 align="center">PacBio Kinnex Iso-Seq Pipeline</h1>

> [!WARNING]
> **BETA:** This workflow is under active development. Interfaces, defaults,
> runtime behavior, and output bundles may change substantially. Use
> tagged releases for reproducible testing.

This repository contains workflows for analyzing human PacBio Kinnex Iso-Seq
sequencing data with [Workflow Description Language
(WDL)](https://openwdl.org/).

- [The wdl-dockerfiles repository](https://github.com/PacificBiosciences/wdl-dockerfiles)
  defines the Docker images, which PacBio hosts on
  [quay.io](https://quay.io/organization/pacbio).
- WDL runtime blocks select task container images. The
  [tools and containers guide](docs/tools_containers.md) lists them.
- Supply reusable GRCh38 and Kinnex resources through an immutable
  [reference container or typed file overrides](docs/reference_container.md).
- User-facing entrypoints live in the [workflows directory](workflows).

## Workflow

The repository has three user-facing workflows.

`kinnex_isoseq` runs end-to-end analysis of Kinnex Iso-Seq HiFi data. It takes
HiFi BAMs from a single acquisition and produces classified, filtered
isoforms.

`preprocessing` generates transcript BAMs, also known as FLNC BAMs, from Kinnex
HiFi data. It takes HiFi BAMs from a single acquisition and writes
sample-annotated transcript BAMs containing poly-A-trimmed, 5' -> 3'
transcript sequences.

`secondary_analysis` takes sample-annotated transcript BAMs from any number of
acquisitions and produces classified, filtered isoforms for each sample named
in the BAM headers. It merges multiple BAMs for the same sample.

| <img src="docs/figures/workflow_diagram.png" width="800" alt="PacBio Kinnex Iso-Seq Pipeline"> |
| :---: |

- [End-to-end HiFi BAM workflow](workflows/kinnex_isoseq.wdl)
- [Preprocessing workflow](workflows/preprocessing.wdl)
- [FLNC-level secondary analysis workflow](workflows/secondary_analysis.wdl)

## Quick start guide

This example uses:

- `miniwdl` as the workflow engine,
- Slurm HPC as the backend,
- Apptainer as the container runtime,
- the PacBio human reference data bundle,
- an instrument acquisition that was not demultiplexed on instrument.

### Setup

- Clone this repository

```bash
git clone git@github.com:PacificBiosciences/Kinnex-IsoSeq-WDL.git
cd Kinnex-IsoSeq-WDL
```

- Install [miniwdl](docs/backend-hpc.md#miniwdl)
- Install [apptainer](https://apptainer.org/docs/admin/main/installation.html#install-from-pre-built-packages)
- Obtain the immutable reference-container URI published for this workflow
  release and ensure compute nodes can authenticate to its registry.

### Fill out your inputs.json

Copy the `no_instrument_demux` kinnex_isoseq template:
[inputs.json](backends/hpc/kinnex_isoseq.no_instrument_demux.hpc.inputs.json)
and the adjacent example
[biosamples.example.csv](backends/hpc/biosamples.example.csv).

Place the copied sample sheet at
`<local_path_prefix>/samples/biosamples.csv`, or update the `biosample_csv`
input to its location. Replace every `<local_path_prefix>` placeholder with
paths visible to the workflow jobs, and edit the sample sheet to describe the
barcode and sample assignments for your run.

You can also generate most run-specific JSON fields directly from SMRT Link
run output by scanning the run directory for `hifi_reads` BAMs:

```bash
python3 scripts/build_inputs_from_files.py \
  kinnex_isoseq \
  --run-dir <local_path_prefix>/r84039_241003_193617 \
  --biosample-csv <local_path_prefix>/samples/biosamples.csv \
  --output-json backends/hpc/kinnex_isoseq.auto.inputs.json
```

If you run the helper from a run directory, it can auto-detect a single
`biosample*.csv` and writes `<workflow>.auto.inputs.json` in the current
directory:

```bash
cd <local_path_prefix>/r84039_241003_193617
python3 /path/to/Kinnex-IsoSeq-WDL/scripts/build_inputs_from_files.py \
  preprocessing
```

The helper fills `hifi_sources`, `biosample_csv`, and `output_prefix` using the
mode-specific HPC template. Use `--mode no_instrument_demux` or
`--mode from_instrument_demux` when you need to force a mode.

For instrument-demultiplexed BAMs, the helper reads the HiFi barcode from the
BAM filename. If the BAMs have one unique barcode and the biosample CSV has one
different HiFi Barcode group, it writes a corrected copy named
`<biosample-name>.updated.csv` and uses that copy in the inputs JSON. The
original sample sheet is preserved. Use `--updated-biosample-csv` to choose a
different output path. When there are multiple unmatched barcode groups, the
helper stops rather than guessing how sample assignments should map.

Optional reference-file flags add or replace fields under the workflow's
`reference_overrides` object: `--hifi-demux-barcodes`, `--indexed-primers`,
`--segmentation-adapters`, `--genome-fasta`, and `--annotation-gtf-gz`. For
example, use `--indexed-primers /refs/custom_primers.fasta` to select custom
primers. `--genome-fasta` requires a matching `--annotation-gtf-gz`. Barcode and
primer flags apply to `preprocessing` and `kinnex_isoseq`; genome and annotation
flags apply to `kinnex_isoseq` and `secondary_analysis`. When custom primers or
adapters are supplied, the helper removes the corresponding packaged-set
selector from the generated JSON.

To create a `secondary_analysis` JSON from existing FLNC BAMs, run from a
directory that contains one or more `*.flnc.bam` files (or set `--run-dir`):

```bash
python3 scripts/build_inputs_from_files.py \
  secondary_analysis
```

### Prepare a preprocessing Slurm job

After each successful run, the generated Slurm script collects final FLNC BAMs
and their `.bam.pbi` indexes into that run's flat `results/` folder. It uses
`samtools view -H` to read the single sample name from the BAM's `@RG SM` tags
and names the copies `<sample>.<original-bam-name>`. Original workflow outputs
are left intact so miniwdl links and call caching remain valid. Existing result
files are not overwritten. Samtools must be available in the job's `PATH`.

For an already completed run, run the collector from its sample folder:

```bash
python3 /path/to/Kinnex-IsoSeq-WDL/scripts/collect_flnc_results.py . \
  --results-dir results
```

The collector uses `_LAST/outputs.json` to select the latest execution. You can
instead pass a specific timestamped execution folder or its `outputs.json`.
Only final outputs listed in that manifest are exported, not intermediate BAMs.

From a parent folder containing one subfolder per run, use the Atlas setup
helper to clone or reuse the workflow repository, generate a separate
preprocessing inputs JSON for each run, and write one Slurm batch script that
processes them in sequence:

```bash
cd /path/to/folder-of-runs
python3 /path/to/Kinnex-IsoSeq-WDL/scripts/create_preprocessing_slurm.py
sbatch run_preprocessing_all.slurm
```

The helper defaults to cloning
`https://github.com/ScottGeib-USDA/Kinnex-IsoSeq-WDL.git` into
`Kinnex-IsoSeq-WDL` under the parent folder. Set `--repo-url` and `--branch` if
your changes are published elsewhere. The selected branch must contain
`scripts/build_inputs_from_files.py`. Each immediate child folder containing a
HiFi BAM is treated as a run and must contain one biosample CSV, unless
`--biosample-csv` supplies a shared CSV for all runs. Each run gets its own
`preprocessing.auto.inputs.json`. The helper writes one Slurm script with the
Atlas account, partition, email, Apptainer setup, Conda environment, and one
miniwdl invocation per run. Failed runs are reported at the end after the
remaining runs have been attempted.

The helper infers the packaged segmentation adapter set from each run folder
name: `8Plex` selects `8-fold`, `12Plex` selects `12-fold`, and `16Plex`
selects `16-fold`. An explicit `--segmentation-adapter-set` applies that value
to every run and overrides folder-name inference. Other supplied resource
options also apply to every run.

The generated job sources Atlas Conda's shell hook from
`/project/ag100pest/scott.geib/miniconda3/etc/profile.d/conda.sh` before running
`conda activate /home/scott.geib/.conda/envs/isoseq_wdl`. Override these paths
with `--conda-sh` and `--conda-env` if your cluster uses a different Conda
installation or environment.

For non-standard barcode or adapter resources, pass the corresponding FASTA
paths to the helper:

```bash
python3 /path/to/Kinnex-IsoSeq-WDL/scripts/create_preprocessing_slurm.py \
  --hifi-demux-barcodes /shared/refs/hifi_barcodes.fasta \
  --indexed-primers /shared/refs/custom_primers.fasta \
  --segmentation-adapters /shared/refs/custom_adapters.fasta
```

These are written to `preprocessing.reference_overrides` in the generated JSON.
The files must be readable from the Slurm compute nodes. A custom indexed-primer
or segmentation-adapter FASTA replaces its corresponding packaged-set selector.

### Run the workflow

```bash
miniwdl run --cfg ~/.config/miniwdl.cfg \
  -i kinnex_isoseq.no_instrument_demux.hpc.inputs.json \
  workflows/kinnex_isoseq.wdl
```

## Beyond quick start

1. [Select a backend environment](#selecting-a-backend)
2. [Configure a workflow execution engine and container runtime](#configuring-a-workflow-engine-and-container-runtime)
3. [Set up your reference resources](#setting-up-your-reference-resources)
4. [Fill out the inputs JSON file for your run](#filling-out-the-inputs-json)
5. [Run the workflow](#running-the-workflow)

### Selecting a backend

Slurm HPC is the supported public backend target. We maintain example
configurations for miniwdl with miniwdl-slurm and for Sprocket with its
Slurm/Apptainer backend. Cromwell can use the same WDL entrypoints and HPC
input templates, but this repository does not provide a Cromwell backend
configuration.

| Backend | Status | Documentation |
| :- | :- | :- |
| HPC | Supported on Slurm with miniwdl, Sprocket, or Cromwell | [HPC backend guide](docs/backend-hpc.md) |
| AWS HealthOmics | No maintained template | Not currently provided |
| Azure | No maintained template | Not currently provided |
| GCP | No maintained template | Not currently provided |

### Configuring a workflow engine and container runtime

You need a workflow engine to run WDL workflows. For the supported public HPC
path, install [`miniwdl`](https://miniwdl.readthedocs.io/en/latest/) with
[`miniwdl-slurm`](https://github.com/miniwdl-ext/miniwdl-slurm),
[`sprocket`](https://github.com/stjude-rust-labs/sprocket), or
[`Cromwell`](https://cromwell.readthedocs.io/en/stable/backends/Backends/).
Make sure Slurm commands and a Singularity-compatible container runtime are
available on the submit host and cluster nodes. See the
[HPC backend guide](docs/backend-hpc.md) for setup details.

### Setting up your reference resources

By default, `ref_name` is set to `"GRCh38_gencode49"`. The workflow resolves
that name to the repository-owned, digest-pinned reference container. The
container supplies the genome, annotation, Pigeon support files, and Kinnex
barcode, adapter, and primer resources without requiring shared-filesystem
resource setup. You can instead supply an explicit immutable
`reference_container` URI, which takes precedence over `ref_name`.

Use `reference_overrides` to replace any user-selectable resource. A run can
use container defaults, combine a container with overrides, or supply every
resource required by its entrypoint.
When all effective resources are overridden, the default named container is
not unpacked. An overridden genome must be an uncompressed FASTA. The workflow
generates its matching FASTA index and requires an explicit matching annotation
override. Packaged Pigeon support files do not fall back for a custom genome:
polyA, CAGE, and junction support are each optional and require their own
override to be used. This rule is especially important for mouse and other
non-GRCh38 genomes. See the
[reference-resource specification](docs/reference_container.md) for validation
rules and cross-resource compatibility responsibilities.

### Filling out the inputs JSON

Each workflow run reads its inputs from JSON. Public input templates live next
to the workflow entrypoints:

- [workflows/kinnex_isoseq.inputs.json](workflows/kinnex_isoseq.inputs.json)
- [workflows/preprocessing.inputs.json](workflows/preprocessing.inputs.json)
- [workflows/secondary_analysis.inputs.json](workflows/secondary_analysis.inputs.json)

HPC input templates are available in the [HPC backend directory](backends/hpc):

- [backends/hpc/kinnex_isoseq.no_instrument_demux.hpc.inputs.json](backends/hpc/kinnex_isoseq.no_instrument_demux.hpc.inputs.json)
- [backends/hpc/kinnex_isoseq.from_instrument_demux.hpc.inputs.json](backends/hpc/kinnex_isoseq.from_instrument_demux.hpc.inputs.json)
- [backends/hpc/preprocessing.hpc.inputs.json](backends/hpc/preprocessing.hpc.inputs.json)
- [backends/hpc/secondary_analysis.hpc.inputs.json](backends/hpc/secondary_analysis.hpc.inputs.json)
- [backends/hpc/biosamples.example.csv](backends/hpc/biosamples.example.csv), a 4 × 12 barcode-combination sample-sheet template for preprocessing and end-to-end runs

For the end-to-end workflow, use the `no_instrument_demux` template for one raw
HiFi BAM, or the `from_instrument_demux` template for BAMs that have been
demultiplexed using the HiFi barcode.

Copy the template for the entrypoint you want to run and replace every
`<local_path_prefix>` placeholder with paths visible to the workflow jobs. The
templates select the digest-pinned `GRCh38_gencode49` reference published for
this workflow release. For preprocessing or end-to-end runs,
`segmentation_adapter_set` selects the packaged `8-fold`, `12-fold`, or `16-fold`
Kinnex segmentation-adapter FASTA; omission selects `8-fold`. The matching
`isoseq_primers_set` choices are `IsoSeq-v2` and `IsoSeq96`, with `IsoSeq-v2`
selected on omission. Do not combine either selector with its matching custom
override (`segmentation_adapters` or `indexed_primers`).

If you do not use an optional input, leave it out of the input JSON. Do not set
it to `null`, because some WDL runners do not accept null values. This also
applies to optional fields inside inputs such as `hifi_sources` and
`reference_overrides`. An empty string or list still counts as a value.

The generic `workflows/*.inputs.json` templates list optional fields so you can
see what is available. Replace an optional placeholder with a real value when
you need it. Otherwise, delete the field. The HPC templates show this pattern.

### Running the workflow

Run with `miniwdl` from the repository root using the input JSON you copied and
edited. See the [HPC backend guide](docs/backend-hpc.md) for additional
engine-specific execution notes.

End-to-end HiFi BAM entrypoint, using the `no_instrument_demux` HPC template:

```bash
miniwdl run --cfg ~/.config/miniwdl.cfg \
  -i backends/hpc/kinnex_isoseq.no_instrument_demux.hpc.inputs.json \
  workflows/kinnex_isoseq.wdl
```

Preprocessing entrypoint:

```bash
miniwdl run --cfg ~/.config/miniwdl.cfg \
  -i backends/hpc/preprocessing.hpc.inputs.json \
  workflows/preprocessing.wdl
```

FLNC-level secondary analysis entrypoint:

```bash
miniwdl run --cfg ~/.config/miniwdl.cfg \
  -i backends/hpc/secondary_analysis.hpc.inputs.json \
  workflows/secondary_analysis.wdl
```

## Workflow inputs

The workflows use three input categories:

- *reference resources* are resolved from a named or explicit immutable
  [`reference_container`](docs/reference_container.md), typed per-file
  `reference_overrides`, or both.
- *sample sheets* are CSV files that describe run-specific sample and barcode
  assignments. `kinnex_isoseq` and `preprocessing` use
  [`biosample_csv`](docs/biosample_csv.md).
- *inputs.json* files define the datasets, sample files, run mode, tuning
  parameters, and backend settings for a workflow run.

The end-to-end `kinnex_isoseq` entrypoint starts from `hifi_sources`.
Omit `hifi_barcode` for HiFi demux mode, or provide one `hifi_barcode` per
source BAM for cDNA demux-only mode.

The standalone `preprocessing` entrypoint uses the same `hifi_sources`,
reference-resource, adapter-selection, and `biosample_csv` inputs as the
end-to-end workflow, and stops after producing FLNC BAMs.

The FLNC-level `secondary_analysis` entrypoint starts from `flnc_bams`. Each
FLNC BAM must contain `@RG` records with exactly one distinct non-empty `SM`
value. The workflow groups BAMs by the exact original `SM` and uses sanitized
sample prefixes for the final aligned BAM filenames. It aligns each input BAM
independently, then merges same-sample aligned BAMs before isocall.

Detailed workflow input and output interfaces are described in:

- [End-to-end HiFi BAM workflow interface](docs/kinnex_isoseq.md)
- [Preprocessing workflow interface](docs/preprocessing.md)
- [FLNC-level secondary analysis workflow interface](docs/secondary_analysis.md)
- [Reference-container contract](docs/reference_container.md)
- [Resource bundle layout](docs/resource_bundle.md)
- [Biosample CSV specification](docs/biosample_csv.md)

## Tool versions and Docker images

Task runtime blocks define the container images used by each workflow task.
Registry-relative PacBio image references, including named reference
containers, use the workflow `container_registry` input when provided and
otherwise default to `"quay.io/pacbio"`. An explicit complete
`reference_container` URI is never rewritten by that setting.

Tool and container details are documented in
[tools and containers](docs/tools_containers.md).

## Resource requirements

Resource requirements depend on input size, sample count, and the entrypoint
being run. Task-level CPU, base memory, retry, and image settings are defined in
WDL runtime blocks.

## Version information

Current development version: **0.4.0**.

For a complete changelog, see the [changelog](CHANGELOG.md) or the git history.

## DISCLAIMER

TO THE GREATEST EXTENT PERMITTED BY APPLICABLE LAW, THIS WEBSITE AND ITS
CONTENT, INCLUDING ALL SOFTWARE, SOFTWARE CODE, SITE-RELATED SERVICES, AND DATA,
ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF
ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY
WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS
FOR A PARTICULAR PURPOSE. ALL WARRANTIES ARE REJECTED AND DISCLAIMED. YOU ASSUME
TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THE FOREGOING. PACBIO IS NOT
OBLIGATED TO PROVIDE ANY SUPPORT FOR ANY OF THE FOREGOING, AND ANY SUPPORT
PACBIO DOES PROVIDE IS SIMILARLY PROVIDED WITHOUT REPRESENTATION OR WARRANTY OF
ANY KIND. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A REPRESENTATION
OR WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE
WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACBIO.
