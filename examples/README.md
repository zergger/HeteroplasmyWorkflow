# Examples

`example_config.txt` and `example_read_ids.txt` mirror the root configuration template. Paired short-read files use `sample01_1.fastq.gz` and `sample01_2.fastq.gz`; single-end files use `sample01.fastq.gz`. Supply real combined/single-organelle references and annotations before running the full workflow.

The older `workflow/` carrot inputs and `visualization/` CSVs are retained historical assets. Their sample paths and old workflow conventions need adaptation. Annotation CSV columns are `Name,Type,Minimum,Maximum,Length,# Intervals,Direction`, with 1-based inclusive coordinates and recognized feature types CDS, rRNA and tRNA.
