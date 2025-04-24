# Expression Atlas Bulk (re)processing and recalculations

This set of Snakemake workflows replaces the Atlas new experiments processing, reprocessing and recalculations operations done from Atlas-Prod codebase which had direct involvement of the LSF CLI and could only run on the original cluster.

It contains data analysis rules for:
- RNA-Seq baseline analysis
- Microarray differential analysis
- RNA-Seq differential analysis
- Proteomics baseline analysis
- Proteomics differential analysis

A recalculations run requires that reprocess has been performed a priori, and it is currently not available for proteomics experiments. Recalculations are the operations that need to happen on load or after an E! Update, and generate a subset of the outputs produced during (re)processing.

## Prerequisites

 * Snakemake (tested with version 7.32.4)
 * LSF or SLURM batch schedulers
 * Set up configuration variables at `run_sorting_hat_test_data.sample.sh` for goal 'reprocess' or 'recalculations'.

## Run pipeline

```
./run_sorting_hat_test_data.sh EXPS_DIR
```

The experiments path contains one or more directories with Atlas accession names E-* (e.g. `E-MTAB-5577`), having at least configuration files in xml format after curation process.

Optionally, worflow execution can be tailored to specific accessions or species by defining these variables in the sorting-hat script.

### New experiment processing and re-processing

Completed processing by iRAP Single Lib (ISL) is necessary before new experiment processing. It will run all rules avilable for the experiment type.

### Recalculations

This is necessary for the Ensembl Update (E! Update) part of a Data Release. E! Update brings annotations from Biomart, E! Mysql databases and ftp sites for all the relevant organisms in Expressiona Atlas, and leaves them in a format that can be consumed for the decoration process and the web applications. Ensembl Update validators (for Biomart attributes and GTF URL validations) are performed before running recalculations.


For differial experiments, the following outputs are generated (which correspond to rules):

1. **Percentile ranks rule**
- Output: `{accession}-percentile-ranks.tsv`

2. **Differential tracks rule**
- Outputs: `{accession}.{contrast_id}.genes.pval.bedGraph`
- Outputs: `{accession}.{contrast_id}.genes.log2foldchange.bedGraph`

3. **Differential GSEA rule**
- Outputs: `{accession}.{contrast_id}.{ext_db}.gsea.tsv`
- Outputs: `{accession}.{contrast_id}.{ext_db}.gsea_list.tsv`

4. **atlas_experiment_summary rule**
- Output: `{accession}-atlasExperimentSummary.Rdata`


For baseline RNA-seq experiments, the following outputs are generated:

1. **Baseline tracks rule**
- Outputs: `{accession}.{assay_id}.genes.expressions_{metric}.bedGraph`
  (where metric could be fpkm/tpm)

2. **Baseline heatmap rule**
- Outputs: `{accession}-heatmap-{metric}.pdf`
- Output: `{accession}-heatmap.pdf`

3. **Baseline coexpression rule**
- Outputs: `{accession}-{metric}-coexpressions.tsv.gz`
- Output: `{accession}-coexpressions.tsv.gz`

4. **atlas_experiment_summary rule**
- Output: `{accession}-atlasExperimentSummary.Rdata`

For proteomics experiments (`proteomics_baseline`, `proteomics_baseline_dia`, `proteomics_differential`), recalculations are not implemented.


| Rule                        |  Baseline  | Differential RNA-seq | Differential microarray |
|-----------------------------|------------|----------------------|-------------------------|
| atlas_experiment_summary    |  &#x2713;  |       &#x2713;       |         &#x2713;        |
| check_differential_gsea     |            |       &#x2713;       |         &#x2713;        |
| differential_gsea           |            |       &#x2713;       |         &#x2713;        |
| differential_tracks         |            |       &#x2713;       |         &#x2713;        |
| percentile_ranks            |            |       &#x2713;       |         &#x2713;        |
| baseline_coexpression       |  &#x2713;  |                      |                         |
| baseline_heatmap            |  &#x2713;  |                      |                         |
| baseline_tracks             |  &#x2713;  |                      |                         |
| link_baseline_coexpression  |  &#x2713;  |                      |                         |
| link_baseline_heatmap       |  &#x2713;  |                      |                         |
| touch_inputs_baseline       |  &#x2713;  |                      |                         |