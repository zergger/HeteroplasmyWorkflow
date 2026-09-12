# Workflow entry points

From the repository root, use `python code/workflow/run_hpc.py config.txt read_ids.txt`. The scripts import sibling modules, read `defaults.ini` and label outputs with `VERSION`. `run_local.py` is retained as a legacy alternative, not the acceptance entry for this release.

Full input conventions, dependencies and parameter meanings are in the [root README](../../README.md).

## Stages

1. `s02_hpc_align.py`: alignment, duplicate removal, coverage-derived mitochondrial threshold and SAM selection.
2. `s03_heteroplasmy_likelihood.py`: base/indel counts, likelihood scores, BH-adjusted q values, coverage requirements and optional within-read linkage.
3. `s04_sort_candidates.py`: sort candidate evidence by score.
4. `s05_select_sites.py`: apply q, Percentage and coverage gates; `d_threshold = 0` disables distance exclusion.
5. `s06_location_conservation.py` and `s07_plot_heteroplasmy.py`: retained conservation/HTML presentation stages.

The four-argument filter CLI is unchanged:

```bash
python code/workflow/filter_samfiles_cp_mt.py input.sam output None MT
python code/workflow/filter_samfiles_rm_numt.py input.sam output None MT
```

The first keeps XA records, the second rejects them; both exclude soft clipping and preserve every incoming SAM header in order in enabled outputs. `None` disables an organelle. These scripts select organelles from an already quality/flag-filtered SAM; they do not perform the upstream MAPQ or proper-pair screen themselves.

## Diagnostics and existing interface limitations

- Missing `en_US.UTF-8`: make that locale available before annotation parsing.
- `invalid mode: rU`: use the declared Python 3.8 environment.
- Bokeh import/API errors: use Bokeh 1.4.0 and compatible Jinja2 from the supplied requirements.
- S03 emits no candidate CSV when no processable candidate exists; handle that state explicitly before sorting/selection.
- Missing or empty per-sample files: inspect both configured logs and the wrapper worker output. Wrapper exit zero is insufficient for acceptance, and file-existence reuse does not validate a previous run's inputs.
- Check input/output paths and create the configured log parent before a run. Avoid whitespace or shell syntax in paths/sample IDs.
- The legacy direct S05 CLI parses its score argument as an integer; use the main orchestration or the native `process` interface for `0.05`.
- `count_threshold` selects a required-coverage column (`95`, `99`, `999`), not an ALT count. Keep `Total > N_req_*`, `Percentage >= threshold` and `q < threshold` distinct.
