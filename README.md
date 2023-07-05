# Expression Atlas Bulk (re)processing and recalculations

This set of Snakemake workflows replaces the Atlas new experiments processing, reprocessing and recalculations operations done from Atlas-Prod codebase which had direct involvement of the LSF CLI and could only run on the original cluster. 

It contains data analysis rules for:
- RNA-Seq baseline analysis
- Microarray differential analysis
- RNA-Seq differential analysis
- Proteomics baseline analysis
- Proteomics differential analysis
- Cell type proportion predictions for RNA-Seq baseline and differential analysis

A recalculations run requires that reprocess has been performed a priori, and it is currently not available for proteomics experiments. Recalculations are the operations that need to happen on load or after an E! Update, and generate a subset of the outputs produced during (re)processing.

## Prerequisites

 * Snakemake (tested with version 6.6.1)
 * LSF batch scheduler
 * Set up configuration variables at `run_sorting_hat_test_data.sample.sh` for goal 'reprocess' or 'recalculations'.

## Run pipeline

```
./run_sorting_hat_test_data.sh EXPS_DIR
```

The experiments path contains one or more directories with Atlas accession names E-* (e.g. E-MTAB-5577), having at least configuration files in xml format after curation process. Completed processing by iRAP Single Lib (ISL) is necessary before new experiment processing.

Optionally, worflow execution can be tailored to specific accessions or species by defining these variables in the sorting-hat script.
