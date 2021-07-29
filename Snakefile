from sys import exit
import yaml

# atom: set grammar=python:

metadata_summary = {}

def read_metadata_summary():
    global metadata_summary
    if not metadata_summary:
        with open(config['metadata_summary'], 'r') as fh:
            metadata_summary = yaml.safe_load(fh)

def get_contrast_labels():
    if 'contrast_labels' in config:
        return f"{config['contrast_labels']}".split("&&")
    else:
        read_metadata_summary()
        return [metadata_summary['contrasts'][x] for x in metadata_summary['contrasts']]

def get_contrast_ids():
    if 'contrast_ids' in config:
        return f"{config['contrast_ids']}".split("::")
    else:
        read_metadata_summary()
        return [x for x in metadata_summary['contrasts']]

def get_assay_labels():
    if 'assay_labels' in config:
        return f"{config['assay_labels']}".split("&&")
    else:
        read_metadata_summary()
        return [metadata_summary['assays'][x] for x in metadata_summary['assays']]

def get_assay_ids():
    if 'assay_ids' in config:
        return f"{config['assay_ids']}".split("::")
    else:
        return [x for x in metadata_summary['assays']]

def get_ext_db():
    if 'ext' in config:
        return f"{config['ext']}".split(":")
    else:
        return ["go", "reactome", "interpro"]

def get_metrics():
    if 'metric' in config:
        return config['metric'].split(":")
    else:
        return ['tpms', "fpkms"]

plot_labels = {"go": "GO terms", "reactome": "Reactome Pathways", "interpro": "Interpro domains"}

def get_ext_db_labels():
    output = []
    for ext_db in get_ext_db():
        output.append(plot_labels[ext_db])
    return output

def check_config_required(fields, method=""):
    exit=False
    for f in fields:
        if f not in config:
            print(f"{f} required to be set in config")
            ex=True
            if method:
                print(f" for method {method}")
    if exit:
        exit(2)

def get_outputs():
    """
    First method to be executed since it is run by rule all.
    """
    import os.path

    tool_outputs = {}
    tool_outputs['percentile-ranks'] = f"{config['accession']}-percentile-ranks.tsv"
    outputs = []

    # Read this now so that it is available for all other needs
    read_metadata_summary()
    required_config=['tool']
    check_config_required(fields=required_config)
    if 'percentile-ranks' in config['tool'] or config['tool']=="all-diff":
        outputs.append(f"{config['accession']}-percentile-ranks.tsv")
    if 'differential-tracks' in config['tool'] or config['tool']=="all-diff":
        check_config_required(fields=['contrast_ids', 'metadata_summary'], method='differential-tracks')
        # fake elements to mix contrasts labels and ids
        outputs.extend(expand(config['accession']+".{id}.{type}", id=get_contrast_ids(), type=["genes.pval.bedGraph", "genes.log2foldchange.bedGraph"]))
    if 'baseline-tracks' in config['tool'] or config['tool']=="all-baseline":
        check_config_required(fields=['metadata_summary'], method='differential-tracks')
        # combine metric (fpkm / tpm) with assay_id/assay_label (zip based)
        # in a product manner
        outputs.extend(expand(config['accession']+".{a_id}.genes.expressions_{metric}.bedGraph",
                            a_id=get_assay_ids(),
                            metric=get_metrics()))
    if 'differential-gsea' in config['tool'] or config['tool']=="all-diff":
        check_config_required(fields=['contrast_ids', 'organism', 'bioentities_properties'], method='differential-gsea')
        outputs.extend(
                expand(config['accession']+".{c_id}.{ext_db}.{type}",
                        c_id=get_contrast_ids(),
                        ext_db=get_ext_db(),
                        type=["gsea.tsv", "gsea_list.tsv"]))
    if 'atlas-experiment-summary' in config['tool'] or 'all' in config['tool']:
        outputs.append(f"{config['accession']}-atlasExperimentSummary.Rdata")
    if 'baseline-heatmap' in config['tool'] or 'all-baseline' in config['tool']:
        outputs.extend(expand(f"{config['accession']}"+"-{metric}.tsv", metric=get_metrics() ))
    print(outputs)
    return outputs

def get_contrast_label(wildcards):
    """
    Meant to be used within a rule, where a specific contrast_id is set in
    the wildcards. It retrieves the contrast id to label from the
    metadata_summary.
    """
    global metadata_summary
    if wildcards['contrast_id'] in metadata_summary['contrasts']:
        return metadata_summary['contrasts'][wildcards['contrast_id']]

def get_ext_db_label(wildcards):
    """
    Meant to be used within a rule, where a specific ext_db is set in
    the wildcards.
    """
    return plot_labels[wildcards['ext_db']]

def get_assay_label(wildcards):
    """
    Meant to be used within a rule, where a specific assay_id is set in
    the wildcards. It retrieves the assay id to label from the
    metadata_summary.
    """
    global metadata_summary
    return metadata_summary['assays'][wildcards['assay_id']]

wildcard_constraints:
    accession="E-\D+-\d+"

rule all:
    input:
        required_outputs=get_outputs()

rule percentile_ranks:
    conda: "envs/atlas-internal.yaml"
    log: "log/{accession}-percentile_ranks.log"
    output:
        percentile_ranks_merged="{accession}-percentile-ranks.tsv"
    shell:
        """
        mkdir -p log
        exec &> "{log}"
        rm -f {wildcards.accession}*-percentile-ranks.tsv
        for analytics in $(ls {wildcards.accession}*-analytics.tsv.unrounded); do
            {workflow.basedir}/bin/calculate_percentile_ranks.R $analytics
        done
        # multiple ranks file will be generated for the microarray case,
        # in that case these need to be merged by gene id (first column),
        # each contrast on a different columns, all tab separated.
        # NA should be added for Gene x contrast for which there is no value
        percentile_ranks=( $(ls {wildcards.accession}*-percentile-ranks.tsv) )
        if [ ${{#percentile_ranks[@]}} -gt 1 ]; then
            # more than one file, requires merging by Gene.ID
            {workflow.basedir}/bin/merge_by_gene_id.R ../{output.percentile_ranks_merged} {wildcards.accession}_*-percentile-ranks.tsv
            # remove only microarray derived multiple percentile ranks
            # (per array design <accession>_<arraydesign>-percentile-ranks.tsv)
            rm -f {wildcards.accession}_*-percentile-ranks.tsv
        else
            mv ${{percentile_ranks[0]}} {output.percentile_ranks_merged}
        fi
        """

rule differential_tracks:
    conda: "envs/irap.yaml"
    log: "log/{accession}.{contrast_id}-differential_tracks.log"
    input:
        gff=config['gff']
        # analytics will be derived below since it could be either {accession}-{arraydesign}-analytics.tsv
        # or just {accession}-analytics.tsv for RNA-Seq
    params:
        contrast_label=get_contrast_label
    output:
        # based on https://stackoverflow.com/questions/58187715/using-the-expand-function-in-snakemake-to-perform-a-shell-command-multiple-tim
        # fake=temp("fake_diff_tracks.{accession}.{contrast_id}")
        pval="{accession}.{contrast_id}.genes.pval.bedGraph",
        log2fold="{accession}.{contrast_id}.genes.log2foldchange.bedGraph"
    shell:
        """
        mkdir -p log
        exec &> "{log}"
        source {workflow.basedir}/bin/tracks_functions.sh
        set +e
        analyticsFile=$(grep -l -P "\\t{wildcards.contrast_id}\." {wildcards.accession}_A-*-analytics.tsv)
        if [ $? -ne 0 ]; then
            # rnaseq case
            analyticsFile={wildcards.accession}-analytics.tsv
        fi
        set -e
        generate_differential_tracks {wildcards.accession} {wildcards.contrast_id} $analyticsFile {input.gff} {params.contrast_label:q} ./
        """

rule differential_gsea:
    conda: "envs/irap.yaml"
    log: "log/{accession}.{contrast_id}.{ext_db}-differential_gsea.log"
    params:
        organism=config['organism'],
        BIOENTITIES_PROPERTIES_PATH=config['bioentities_properties'],
        contrast_label=get_contrast_label,
        ext_db_label=get_ext_db_label
    output:
        gsea="{accession}.{contrast_id}.{ext_db}.gsea.tsv",
        gsea_list="{accession}.{contrast_id}.{ext_db}.gsea_list.tsv"
    shell:
        """
        mkdir -p log
        exec &> "{log}"
        export BIOENTITIES_PROPERTIES_PATH={params.BIOENTITIES_PROPERTIES_PATH}
        source {workflow.basedir}/bin/gsea_functions.sh
        set +e
        analyticsFile=$(grep -l -P "\\t{wildcards.contrast_id}\." {wildcards.accession}_A-*-analytics.tsv)
        if [ $? -ne 0 ]; then
            # rnaseq case
            analyticsFile={wildcards.accession}-analytics.tsv
        fi
        set -e
        pvalColNum=$(get_contrast_colnum $analyticsFile {wildcards.contrast_id} "p-value")
        log2foldchangeColNum=$(get_contrast_colnum $analyticsFile {wildcards.contrast_id} "log2foldchange")
        plotTitle="
        Top 10 {params.ext_db_label} enriched in
        {params.contrast_label}
        (Fisher-exact, FDR < 0.1)"
        annotationFile=$(find_properties_file {params.organism} {wildcards.ext_db})
        {workflow.basedir}/bin/gxa_calculate_gsea.sh {wildcards.accession} $annotationFile $analyticsFile $pvalColNum $log2foldchangeColNum ./ {wildcards.contrast_id} "$plotTitle" {params.organism} {wildcards.ext_db}
        """

rule baseline_tracks:
    conda: "envs/irap.yaml"
    log: "log/{accession}-{assay_id}-{metric}-baseline_tracks.log"
    params:
        assay_label=get_assay_label
    input:
        gff=config['gff'],
        analytics="{accession}-{metric}.tsv"
    output:
        bedGraph="{accession}.{assay_id}.genes.expressions_{metric}.bedGraph"
    shell:
        """
        mkdir -p log
        exec &> {log:q}
        source {workflow.basedir}/bin/tracks_functions.sh
        echo "Past sourcing"
        generate_baseline_tracks {wildcards.accession} {wildcards.assay_id} {input.analytics} {input.gff} ./ {params.assay_label:q}
        """

rule baseline_coexpression:
    conda: "envs/clusterseq.yaml"
    log: "log/{accession}-{metric}-baseline_coexpression.log"
    input:
        expression="{accession}-{metric}.tsv.undecorated.aggregated"
    output:
        coexpression_comp="{accession}-{metric}-coexpressions.tsv.gz"
    shell:
        """
        mkdir -p log
        exec &> "{log}"
        {workflow.basedir}/run_coexpression_for_experiment.R {input.expression} {output.coexpression_comp}
        """

rule link_baseline_coexpression:
    """
    There is a case where coexpression might not be calculated, when the dataset
    has less than 3 columns. In that case it might be that the input files for this
    never appear, not sure whether this will timeout without errors or not.
    """
    log: "{accession}-{metric}-link_baseline_coexpression.log"
    input:
        expand("{accession}-{metric}-coexpressions.tsv.gz", metric=get_metrics(), accession=["{accession}"])
    shell:
        """
        mkdir -p log
        exec &> "{log}"
        if [ -s {wildcards.accession}-tpm-coexpressions.tsv.gz ]; then
            ln -s {wildcards.accession}-tpm-coexpressions.tsv.gz {wildcards.accession}-coexpressions.tsv.gz
        elif [ -s {wildcards.accession}-fpkm-coexpressions.tsv.gz ]; then
            ln -s {wildcards.accession}-fpkm-coexpressions.tsv.gz {wildcards.accession}-coexpressions.tsv.gz
        else
            echo "Error: neither TPM nor FPKM coexpressions.tsv.gz file found"
            exit 1
        fi
        """


rule baseline_heatmap:
    conda: "envs/atlas-internal.yaml"
    log: "log/{accession}-{metric}-baseline_heatmap.log"
    input:
        expression="{accession}-{metric}.tsv"
    output:
        heatmap="{accession}-heatmap-{metric}.pdf"
    shell:
        """
        mkdir -p log
        exec &> "{log}"
        {workflow.basedir}/bin/generateBaselineHeatmap.R --configuration {wildcards.accession}-configuration.xml \
		--input  input.expression \
		--output output.heatmap
        """

rule atlas_experiment_summary:
    conda: "envs/atlas-internal.yaml"
    log: "log/{accession}-atlas_experiment_summary.log"
    output:
        rsummary="{accession}-atlasExperimentSummary.Rdata"
    shell:
        """
        mkdir -p log
        exec &> "{log}"
        {workflow.basedir}/bin/createAtlasExperimentSummary.R \
	          --source ./ \
	          --accession {wildcards.accession} \
	          --output {output.rsummary}
        """
