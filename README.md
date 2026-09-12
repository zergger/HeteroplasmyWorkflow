# icHET 1.0

icHET detects candidate organellar heteroplasmy from aligned sequencing reads and produces per-site evidence tables. This repository extends the original [HeteroplasmyWorkflow](https://github.com/vtphan/HeteroplasmyWorkflow) with a read/alignment screen, multiple-testing-adjusted scores, an adaptive mitochondrial Percentage threshold, coverage-sufficiency screens, and explicit indel/linkage output.

The current workflow is in [`code/workflow`](code/workflow), including the SAM-header repair in the XA-retaining filter.

## Installation

The tested calling runtime is **Python 3.8.13**. **Python 3.11 and later are not supported** by this release. Use the supplied environment and dependencies, including Bokeh 1.4.0 for the retained HTML visualization.

```bash
git clone --branch main https://github.com/zergger/HeteroplasmyWorkflow.git
cd HeteroplasmyWorkflow
conda env create --file environment.yml
conda activate ichet-1.0
```

[`requirements.txt`](requirements.txt) records the tested Python package versions; [`environment.yml`](environment.yml) also declares bwa-mem2, samtools, sambamba and minimap2. `bash install_packages.sh` is an equivalent environment-creation entry. On Linux, make the `en_US.UTF-8` locale available: the existing annotation parser requests it explicitly. The release validation uses an existing compatible environment; a fresh Conda solve is not itself a tested release result.

## Run from FASTQ

Prepare the following inputs before invoking the retained orchestration entry:

- A combined nuclear and organellar reference, with an exact organelle sequence-ID match between the combined FASTA, single-organellar FASTA and configuration. Reference composition affects ambiguous alignment and threshold interpretation.
- A single-record mitochondrial FASTA and annotation CSV. Annotation columns are `Name,Type,Minimum,Maximum,Length,# Intervals,Direction`; coordinates are 1-based inclusive. Recognized annotation types are `CDS`, `rRNA` and `tRNA`.
- Paired short reads named `sample01_1.fastq.gz` and `sample01_2.fastq.gz`, and one sample ID per line in `read_ids.txt`. Single-end reads use `sample01.fastq.gz` with `PE = 0`.

Edit [`config.txt`](config.txt), keeping sample names and paths free of whitespace and shell metacharacters because the retained wrappers build shell commands. Run from the repository root:

```bash
mkdir -p output
bwa-mem2 index data/reference/combined.fasta
python code/workflow/run_hpc.py config.txt read_ids.txt
```

Use a fresh output directory for a new input/configuration. Existing wrappers reuse intermediate files by existence and can suppress worker exceptions; inspect `LOG_FILE`, `code/workflow/log_align_analyze_sort.txt`, every expected sample output and its contents instead of relying only on the process exit status. See [workflow details and diagnostics](code/workflow/README.md).

`read_type = short` uses bwa-mem2. A retained experimental `read_type = long` route uses minimap2 with `mm_preset` and single-end reads. `threads` controls threads per external command and `max_processes` controls concurrent sample workers.

## Filtering and interpretation

| Setting or output | Meaning in this implementation |
| --- | --- |
| `rm_alter_align = 1` | Default organelle filter rejects soft-clipped records and records containing `XA:`. XA is an alignment-ambiguity signal, not molecular proof of NUMT origin. |
| `rm_alter_align = 0` | Retains XA-bearing records while still rejecting soft clipping; both current filters preserve SAM headers. Treat this as a sensitivity branch. |
| `alignment_quality = 30` | Minimum MAPQ for the upstream samtools screen; paired short reads also require proper pairing and exclude secondary/supplementary records. |
| `score_threshold = 0.05` | Select `Score_q_value < 0.05`; the q value derives from the implemented candidate family. |
| `percentage_threshold = 0.01` | Fallback threshold. A per-sample threshold file overrides it. For mitochondrial runs the wrapper derives `1 / (R + 1)` from its samtools-depth summaries, where `R` is mean mitochondrial/mean nuclear depth. |
| `Percentage` | Minor-category fraction in the six-category A/C/G/T/D/I profile; it is not always the frequency of one specified ALT allele. |
| `count_threshold = 95` | Select `Total > N_req_95`. Values `99` and `999` choose the corresponding columns; this parameter is not a minimum ALT-read count. |
| `d_threshold = 0` | Disables nearest-neighbor-distance exclusion. |

Zero mean depth in either depth summary gives a wrapper threshold of zero and requires inspection. A sample-wide depth ratio is not a locus-specific contamination estimate. Retained calls and indel class counts do not establish molecular ancestry, event-specific allele frequency, biological truth, or a universal detection limit. Linear-reference boundaries, amplification bias and read filtering can affect recovery.

## Outputs and provenance

Per-organellar directories contain filtered SAMs, threshold files, candidate `csv/` tables (including category counts, q values and `N_req_*` columns), optional linkage tables, and `Result/` call summaries. Coordinates in exported candidate and final tables are 1-based. Record the Git commit, configuration, reference identities, tool versions, sample list and logs for every run. Historical wrapper output names include dates and random run IDs, so whole-run filename/HTML byte identity is not guaranteed.

The older carrot examples and visualization assets remain in the repository for their historical use; they are not the current 1.0 acceptance example. See [`examples/README.md`](examples/README.md), [`CHANGELOG.md`](CHANGELOG.md) and the existing [`LICENSE`](LICENSE).

## Citation

In preparation.

## References

Phan V, Pham DT, Melton C, Ramsey AJ, Daigle BJ Jr, Mandel JR. (2019). icHET: interactive visualization of cytoplasmic heteroplasmy. *Bioinformatics*, 35(21), 4411–4412.
