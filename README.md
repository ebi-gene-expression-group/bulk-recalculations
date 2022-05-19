# Atlas Bulk (re)processing and recalculations

This set of Snakemake workflows aims to replace the Atlas new experiments processing, reprocessing and recalculations operations done from Atlas-Prod codebase which had direct involvement of the LSF CLI and could only run on the original cluster. Recalculations are the operations that need to happen on load or after an E! Update.

It contains data analysis rules for:
- RNA-Seq baseline analysis
- Microarray differential analysis
- RNA-Seq differential analysis

## Prerequisites

 * [Snakemake](https://snakemake.readthedocs.io/) installed
 * LSF batch scheduler
 * Set up configuration variables at `run_sorting_hat_test_data.sample.sh` for goal 'reprocess' or 'recalculations'.

## Run pipeline

```
./run_sorting_hat_test_data.sh $experiments_dir
```

The experiments path contains directories with Atlas accession names E-* (e.g. E-MTAB-5577), containing at least congifuration files in xml format after curation process.

A recalculations run requires that reprocess has been performed a priori. Recalaculations out produces a subset of the outputs produced during (re)processing.
