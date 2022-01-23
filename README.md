# Atlas Bulk Recalculations

This set of Snakemake workflows aims to replace the Atlas recalculations operations done from Atlas-Prod codebase which had direct involvement of the LSF CLI and could only run on the original cluster. Recalculations are the operations that need to happen on load or after an E! Update.

It deals with:
- Microarray recalculations
- RNA-Seq differential recalculations
- RNA-Seq baseline recalculations

Ongoing features to handle new experiments processing / reprocessing. Reprocess RNA-seq experiments (baseline or differential), to be used when the experiment has been reprocessed in iRAP Single Lib (ISL), for instance due to a change in assembly or annotation.


## Temporary comments

Run like:

```
atlas-bulk-recalculations/test-data/E-MTAB-5577 on  master [?] via 🅒 snakemake on ☁️
❯ snakemake --config accession="E-MTAB-5577" tool="all" --use-conda --conda-frontend mamba -j 2 -s ../../Snakefile
```

```
snakemake --dry-run --config accession="E-MTAB-5577" tool="all" contrast_ids="g1_g2:g1_g3" contrast_labels="'PATL1 siRNA' vs 'HBB siRNA (control)':'C1' vs 'C2'" gff_file="../gff/organism.gtf" -j 2 -s ../../Snakefile
```

Example microarray differential run, standing on analysis/E-MTAB-5577 (test-data/E-MTAB-5577 in this case)

```
snakemake --dry-run --config accession="E-MTAB-5577" tool="all-diff" contrast_ids="g1_g2:g1_g3" contrast_labels="'PATL1 siRNA' vs 'HBB siRNA (control)':'C1' vs 'C2'" gff_file="../gff/organism.gtf" -j 2 -s ../../Snakefile
```

```
snakemake --dry-run --config accession="E-MTAB-5577" tool="all-diff" contrast_ids="g1_g2:g1_g3" contrast_labels="'PATL1 siRNA' vs 'HBB siRNA (control)':'C1' vs 'C2'" gff_file="../gff/organism.gtf" organism='arabidopsis_thaliana' BIOENTITIES_PROPERTIES_PATH=../bioentities -j 2 -s ../../Snakefile
```

```
snakemake --dry-run --config accession="E-MTAB-5577" tool="all-baseline" assay_ids="g1:g1" assay_labels="Assay1:Assay2" gff_file="../gff/organism.gtf" metric="tpm:fpkm" organism='arabidopsis_thaliana' -j 2 -s ../../Snakefile
```

To run sorting hat, which will produce running calls for all studies in the working directory:
```
snakemake --use-conda --conda-frontend mamba --config gtf_dir=gff sm_options="--use-conda --conda-frontend mamba -j 2" -j 1 -s ../Snakefile-sorting-hat
```
