from sys import exit
import yaml
import os

os.makedirs('logs', exist_ok='True')

# atom: set grammar=python:

metadata_summary = {}

def read_metadata_summary():
    global metadata_summary
    if not metadata_summary:
        with open(config['metadata_summary'], 'r') as fh:
            metadata_summary = yaml.safe_load(fh)

read_metadata_summary()


def get_methods_template_baseline():
    if 'methods_base' in config:
        return config['methods_base']
    else:
        return None

def get_methods_template_differential():
    if 'methods_dif' in config:
        return config['methods_dif']
    else:
        return None

def get_isl_dir():
    if 'isl_dir' in config:
        return config['isl_dir']
    else:
        return None

def get_tmp_dir():
    if 'tmp_dir' in config:
        return config['tmp_dir']
    else:
        return None

def read_skip_steps_file():
    if 'skip_steps_file' in config:
        global skip_steps
        with open(config['skip_steps_file'], 'r') as stream:
            try:
                skip_steps=yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                print(exc)
        return skip_steps


def get_from_config_or_metadata_summary(label):
    if label in config:
        return config[label]
    else:
        return metadata_summary[label]

def get_sdrf():
    return get_from_config_or_metadata_summary('sdrf')

def get_gff():
    return get_from_config_or_metadata_summary('gff')

def get_organism():
    return get_from_config_or_metadata_summary('organism')

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

def get_metrics_recalculations():
    import glob
    if 'metric' in config:
        return config['metric'].split(":")
    else:
        metric_grabbed = []
        for metric in ['tpms', 'fpkms']:
            # glob returns a list of matching paths, so if the metric is available a non-empty list is returned.
            if glob.glob(f"{config['accession']}-{metric}.tsv"):
                metric_grabbed.append(metric)
        if (len(metric_grabbed )>0):
            return metric_grabbed     # ideally: ['tpms', "fpkms"]
        else:
            sys.exit("No metric available for baseline analyses.")

def get_metrics_reprocess():
    """
    This could be based on files retrieved from iRAP/ISL
    """
    #import glob
    if 'metric' in config:
        return config['metric'].split(":")
    else:
        return ['tpms', 'fpkms']



#metrics = get_metrics()
plot_labels = {"go": "GO terms", "reactome": "Reactome Pathways", "interpro": "Interpro domains"}

def get_ext_db_labels():
    output = []
    for ext_db in get_ext_db():
        output.append(plot_labels[ext_db])
    return output

def check_config_required(fields, method=""):
    ex=False
    for f in fields:
        if f not in config:
            print(f"{f} required to be set in config")
            ex=True
            if method:
                print(f" for method {method}")
    if ex:
        exit(2)

def skip(acc, tool):
    if skip_accession != None and acc in skip_accession['skips'][tool]:
        return False
    else:
        return True

def get_number_columns(tsvfile):
    import csv
    with open(tsvfile) as file:
        tsv_file = csv.reader(file, delimiter="\t")
        num_cols = len(next(tsv_file))
    # return number of columns in a tsv file
    return num_cols

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


def get_mem_mb(wildcards, attempt):
    """
    To adjust resources in the rules 
    attemps = reiterations + 1
    Max number attemps = 8
    """
    mem_avail = [ 2, 2, 4, 8, 16, 64, 128, 256 ]  
    return mem_avail[attempt-1] * 1000



localrules: check_differential_gsea, link_baseline_coexpression, link_baseline_heatmap, copy_raw_gene_counts_from_isl, copy_normalised_counts_from_isl, copy_transcript_files_from_isl, copy_transcript_relative_isoforms, create_tracks_symlinks, check_mvaPlot_rnaseq


wildcard_constraints:
    accession="E-\D+-\d+",
    metric="tpms|fpkms"

rule percentile_ranks:
    conda: "envs/atlas-internal.yaml"
    log: "logs/{accession}-percentile_ranks.log"
    resources: mem_mb=get_mem_mb
    output:
        percentile_ranks_merged="{accession}-percentile-ranks.tsv"
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
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
            {workflow.basedir}/bin/merge_by_gene_id.R {output.percentile_ranks_merged} {wildcards.accession}_*-percentile-ranks.tsv
            # remove only microarray derived multiple percentile ranks
            # (per array design <accession>_<arraydesign>-percentile-ranks.tsv)
            rm -f {wildcards.accession}_*-percentile-ranks.tsv
        else
            if [[ "${{percentile_ranks[0]}}" != "{output.percentile_ranks_merged}" ]]; then
                echo "microarray experiment"
                mv ${{percentile_ranks[0]}} {output.percentile_ranks_merged}
            fi
        fi
        """

rule differential_tracks:
    conda: "envs/irap.yaml"
    log: "logs/{accession}.{contrast_id}-differential_tracks.log"
    resources: mem_mb=get_mem_mb
    input:
        gff=get_gff()
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
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
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
    log: "logs/{accession}.{contrast_id}.{ext_db}-differential_gsea.log"
    resources: mem_mb=get_mem_mb
    threads: 8
    params:
        organism=get_organism(),
        BIOENTITIES_PROPERTIES_PATH=config['bioentities_properties'],
        contrast_label=get_contrast_label,
        ext_db_label=get_ext_db_label
    output:
        gsea="{accession}.{contrast_id}.{ext_db}.gsea.tsv",
        gsea_list="{accession}.{contrast_id}.{ext_db}.gsea_list.tsv"
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
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
        annotationFile=$(find_properties_file {params.organism} {wildcards.ext_db})
        if [ -s "$annotationFile" ]; then
            pvalColNum=$(get_contrast_colnum $analyticsFile {wildcards.contrast_id} "p-value")
            log2foldchangeColNum=$(get_contrast_colnum $analyticsFile {wildcards.contrast_id} "log2foldchange")
            plotTitle="
            Top 10 {params.ext_db_label} enriched in
            {params.contrast_label}
            (Fisher-exact, FDR < 0.1)"
            {workflow.basedir}/bin/gxa_calculate_gsea.sh {wildcards.accession} $annotationFile $analyticsFile $pvalColNum $log2foldchangeColNum ./ {wildcards.contrast_id} "$plotTitle" {params.organism} {wildcards.ext_db} {threads}
        else
            touch {wildcards.accession}.{wildcards.contrast_id}.{wildcards.ext_db}.gsea.tsv
            touch {wildcards.accession}.{wildcards.contrast_id}.{wildcards.ext_db}.gsea_list.tsv
        fi
        """


rule check_differential_gsea:
    """
    Whether there is a lack of annotation for 'ext_db' in rule differential_gsea, or the
    analysis produces no results passing significance cutoff, gsea output files could be empty.
    Here we aim to delete those empty files.
    """
    log: "logs/{accession}.{contrast_id}.{ext_db}-check-differential_gsea.log"
    input:
        gsea="{accession}.{contrast_id}.{ext_db}.gsea.tsv",
        gsea_list="{accession}.{contrast_id}.{ext_db}.gsea_list.tsv"
    output:
        temp_gsea=temp("logs/{accession}.{contrast_id}.{ext_db}.check_differential_gsea.done"),
        temp_gsea_list=temp("logs/{accession}.{contrast_id}.{ext_db}.check_differential_gsea_list.done")
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        if [ -e {input.gsea} ]; then
            if ! [ -s {input.gsea} ] || [ $(sed -n '$=' {input.gsea} ) -lt 2 ] ; then
                echo "Deleting empty file "{input.gsea}  >> {log}
                rm {input.gsea}
            else
                echo "File "{input.gsea}" is ok"  >> {log}
            fi
        else
            echo "File "{input.gsea}" does not exist"  >> {log}
        fi

        if [ -e {input.gsea_list} ]; then
            if ! [ -s {input.gsea_list} ] || [ $(sed -n '$=' {input.gsea_list} ) -lt 2 ] ; then
                echo "Deleting empty file "{input.gsea_list}  >> {log}
                rm {input.gsea_list}
            else
                echo "File "{input.gsea_list}" is ok"  >> {log}
            fi
        else
            echo "File "{input.gsea_list}" does not exist"  >> {log}
        fi
        touch {output.temp_gsea}
        touch {output.temp_gsea_list}
        """


rule baseline_tracks:
    conda: "envs/irap.yaml"
    log: "logs/{accession}-{assay_id}-{metric}-baseline_tracks.log"
    resources: mem_mb=get_mem_mb
    params:
        assay_label=get_assay_label
    input:
        gff=get_gff(),
        analytics="{accession}-{metric}.tsv"
    output:
        bedGraph="{accession}.{assay_id}.genes.expressions_{metric}.bedGraph"
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> {log:q}
        source {workflow.basedir}/bin/tracks_functions.sh
        echo "Past sourcing"
        generate_baseline_tracks {wildcards.accession} {wildcards.assay_id} {input.analytics} {input.gff} ./ {params.assay_label:q}
        """

rule baseline_coexpression:
    conda: "envs/clusterseq.yaml"
    log: "logs/{accession}-{metric}-baseline_coexpression.log"
    resources: mem_mb=get_mem_mb
    params: num_retries=5
    input:
        expression="{accession}-{metric}.tsv.undecorated.aggregated"
    threads: 16
    output:
        coexpression_comp="{accession}-{metric}-coexpressions.tsv.gz"
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        {workflow.basedir}/bin/run_coexpression_for_experiment.R {input.expression} {output.coexpression_comp} {workflow.basedir} {threads} {params.num_retries}
        """

rule link_baseline_coexpression:
    """
    There is a case where coexpression might not be calculated, when the dataset
    has less than 3 columns. In that case it might be that the input files for this
    never appear, not sure whether this will timeout without errors or not.
    """
    log: "logs/{accession}-link_baseline_coexpression.log"
    input: lambda wildcards: f"{wildcards.accession}-tpms-coexpressions.tsv.gz" if 'tpms' in metrics else f"{wildcards.accession}-fpkms-coexpressions.tsv.gz"
    output: "{accession}-coexpressions.tsv.gz"
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        ln -s {input} {output}
        """

rule baseline_heatmap:
    conda: "envs/atlas-internal.yaml"
    log: "logs/{accession}-{metric}-baseline_heatmap.log"
    resources: mem_mb=get_mem_mb
    input:
        expression="{accession}-{metric}.tsv"
    output:
        heatmap="{accession}-heatmap-{metric}.pdf"
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        {workflow.basedir}/bin/generateBaselineHeatmap.R --configuration {wildcards.accession}-configuration.xml \
		--input  {input.expression} \
		--output {output.heatmap}
        """

rule link_baseline_heatmap:
    log: "logs/{accession}-link_baseline_heatmap.log"
    input: lambda wildcards: f"{wildcards.accession}-heatmap-tpms.pdf" if 'tpms' in metrics else f"{wildcards.accession}-heatmap-fpkms.pdf"
    output: "{accession}-heatmap.pdf"
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        ln -s {input} {output}
        """

rule atlas_experiment_summary:
    conda: "envs/atlas-internal.yaml"
    log: "logs/{accession}-atlas_experiment_summary.log"
    resources: mem_mb=get_mem_mb
    input:
        sdrf=get_sdrf()
    output:
        rsummary="{accession}-atlasExperimentSummary.Rdata"
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        export SDRF_PATH={input.sdrf}
        {workflow.basedir}/bin/createAtlasExperimentSummary.R \
	          --source ./ \
	          --accession {wildcards.accession} \
	          --output {output.rsummary}
        """



# rules below are specific for reprocessing

# baseline_rnaseq_experiment

rule copy_raw_gene_counts_from_isl:
    """
    Copy raw gene counts file.
    """
    log: "logs/{accession}-copy_raw_gene_counts_from_isl.log"
    output:
        raw_counts_undecorated="{accession}-raw-counts.tsv.undecorated"
    params:
        organism=get_organism(),
        isl_dir=get_isl_dir()
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        expIslDir={params.isl_dir}/{wildcards.accession}/{params.organism}
        echo "ISL dir: $expIslDir"

        [ ! -z $expIslDir+x ] || (echo "Env var $expIslDir needs to defined" && exit 1)
        if [ -s "$expIslDir/genes.raw.htseq2.tsv" ]; then
            rsync -avz $expIslDir/genes.raw.htseq2.tsv {output.raw_counts_undecorated}
        elif [ -s "$expIslDir/genes.raw.featurecounts.tsv" ]; then
            rsync -avz $expIslDir/genes.raw.featurecounts.tsv {output.raw_counts_undecorated}
        else
            echo "Neither genes.raw.htseq2.tsv nor genes.raw.featurecounts.tsv found on $expIslDir"
            exit 1
        fi
        """


rule copy_normalised_counts_from_isl:
    """
    Copy fpkm and tpm gene expression files.
    Replaces copy_unit_matrices_from_isl in experiment_loading_routines.sh
    """
    log: "logs/{accession}-{metric}-copy_normalised_counts_from_isl.log"
    params:
        organism=get_organism(),
        isl_dir=get_isl_dir()  
    output:
        normalised_counts_undecorated="{accession}-{metric}.tsv.undecorated"
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        expIslDir={params.isl_dir}/{wildcards.accession}/{params.organism}
        echo "ISL dir: $expIslDir"

        [ ! -z $expIslDir+x ] || (echo "snakemake param exp_isl_dir needs to defined in rule" && exit 1)

        if [[ "{wildcards.metric}" == "tpms" ]]; then
            metrictype="tpm"
        else
            metrictype="fpkm"
        fi

        if [ -s "$expIslDir/genes.$metrictype.htseq2.tsv" ]; then
            rsync -avz $expIslDir/genes.$metrictype.htseq2.tsv {output.normalised_counts_undecorated}
        elif [ -s "$expIslDir/genes.$metrictype.featurecounts.tsv" ]; then
            rsync -avz $expIslDir/genes.$metrictype.featurecounts.tsv {output.normalised_counts_undecorated}
        else
            echo "$expIslDir/genes.$metrictype.htseqORfeaturecounts.tsv not found"
            exit 1
        fi
        """


rule copy_transcript_files_from_isl:
    """
    This rule attemps to copy Kallisto TPM transcripts if metrics 'tpms' exists.
    If file does not exist for an accession, this rule can be skipped.
    """
    log: "logs/{accession}-copy_transcript_files_{metric}_from_isl.log"
    params:
        organism=get_organism(),
        isl_dir=get_isl_dir()  
    output:
        transcripts="{accession}-transcripts-{metric}.tsv.undecorated"
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        expIslDir={params.isl_dir}/{wildcards.accession}/{params.organism}
        echo "ISL dir: $expIslDir"

        [ ! -z $expIslDir+x ] || (echo "snakemake param exp_isl_dir needs to defined in rule" && exit 1)

        if [ -s "$expIslDir/transcripts.tpm.kallisto.tsv" ] ; then
            rsync -avz $expIslDir/transcripts.tpm.kallisto.tsv {output.transcripts}
        else
            echo "$expIslDir/transcripts.tpm.kallisto.tsv not found - skipping"
        fi
        """

rule copy_transcript_relative_isoforms:
    """
    Copy transcripts relative isoform usage files.
    If file does not exist for an accession, this rule can be skipped.
    """
    log: "logs/{accession}-copy_transcript_relative_isoforms.log"
    output:
        temp("logs/{accession}-copy_transcript_relative_isoforms.done")
    params:
        transcripts_relative_isoforms="{accession}-transcripts.riu.tsv",
        organism=get_organism(),
        isl_dir=get_isl_dir()  
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        expIslDir={params.isl_dir}/{wildcards.accession}/{params.organism}
        echo "ISL dir: $expIslDir"

        [ ! -z $expIslDir+x ] || (echo "snakemake param exp_isl_dir needs to defined in rule" && exit 1)

        if [ -s "$expIslDir/transcripts.riu.kallisto.tsv" ] ; then
            rsync -avz $expIslDir/transcripts.riu.kallisto.tsv {params.transcripts_relative_isoforms}
        else
            echo "$expIslDir/transcripts.riu.kallisto.tsv not found for {wildcards.accession} - skipping"
        fi
        touch {output}
        """

rule rnaseq_qc:
    """
    QC step for baseline rnaseq experiments.
    (WIP)
    """
    conda: "envs/perl-atlas-modules.yaml"
    log: "logs/{accession}-rnaseq_qc.log"
    output: "qc/{accession}-irap-single-lib-report.tsv"
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"

        {workflow.basedir}/bin/rnaseqQC.sh {wildcards.accession} {workflow.basedir}
        qcExitCode=$?

        if [ "$qcExitCode" -eq 2 ]; then
            echo "Experiment {wildcards.accession} has been disqualified due to insufficient quality, exiting"
            exit 1
        elif [ "$qcExitCode" -ne 0 ]; then
            echo "ERROR: QC for {wildcards.accession} failed" >&2
            exit "$qcExitCode"
        fi
        """

# move_qc_folder
# mv ${expTargetDir}/.qc ${expTargetDir}/qc

rule quantile_normalise_expression:
    """
    Quantile normalize and summarize expression in tpms and fpkms.
    To maintain previous logic, the output file is marked as temporal
    """
    conda: "envs/quantile.yaml"
    log: "logs/{accession}-quantile_normalise_expression_{metric}.log"
    input:
        config_xml="{accession}-configuration.xml",
        expression="{accession}-{metric}.tsv.undecorated"
    output:
        qn_expression=temp("{accession}-{metric}.tsv.undecorated.quantile_normalized")
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        {workflow.basedir}/bin/quantile_normalize.sh  -c {input.config_xml} -s {input.expression} -d {output.qn_expression} -b {workflow.basedir}/bin
        if [ $? -ne 0 ]; then
            echo "ERROR: Failed to quantile normalize {wildcards.metric} for {wildcards.accession}  " >&2
            exit 1
        fi
        """

rule summarize_expression:
    """
    Summarize expression, either into median per biological replicate or into quartile per assay group.
    """
    conda: "envs/perl-atlas-modules.yaml"
    log: "logs/{accession}-summarize_expression_{metric}.log"
    input:
        config_xml="{accession}-configuration.xml",
        qn_expression="{accession}-{metric}.tsv.undecorated.quantile_normalized"
    output:
        sum_expression="{accession}-{metric}.tsv.undecorated.aggregated"
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"

        perl {workflow.basedir}/bin/gxa_summarize_expression.pl  \
            --aggregate-quartiles \
            --configuration {input.config_xml}  \
            < {input.qn_expression}  \
            > {output.sum_expression}
        if [ $? -ne 0 ]; then
	        echo "ERROR: Failed to summarize gene expression in {wildcards.metric} for {wildcards.accession} " >&2
	        exit 1
        fi
        """

rule transcripts_na_check:
    """
    Replace NAs with 0 in Kallisto TPM transcripts, if the file exists
    (input and output file here is the same)
    """
    conda: "envs/atlas-internal.yaml"
    log: "logs/{accession}-rule-transcripts_na_check_{metric}.log"
    params:
        organism=get_organism(),
        isl_dir=get_isl_dir()  
    input:
        transcripts="{accession}-transcripts-{metric}.tsv.undecorated"
    output:
        temp("logs/{accession}-rule-transcripts_na_check_{metric}.done")
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        expIslDir={params.isl_dir}/{wildcards.accession}/{params.organism}
        echo "ISL dir: $expIslDir"

        [ ! -z $expIslDir+x ] || (echo "snakemake param exp_isl_dir needs to defined in rule" && exit 1)

        if [ -s "$expIslDir/transcripts.raw.kallisto.tsv" ] ; then
            {workflow.basedir}/bin/transcripts_expr_values_check.R {input.transcripts} $expIslDir/transcripts.raw.kallisto.tsv
            echo "transcripts NA check -  executed for {input.transcripts} "
        else
            echo "$expIslDir/transcripts.raw.kallisto.tsv not found for {wildcards.accession} - skipping rule_transcripts_na_check for {input.transcripts}"
        fi
        touch {output}
        """


rule quantile_normalise_transcripts:
    """
    Quantile normalize transcripts in tpms, if the file exists.
    """
    conda: "envs/quantile.yaml"
    log: "logs/{accession}-quantile_normalise_transcripts_{metric}.log"
    input:
        xml="{accession}-configuration.xml",
        transcripts_na_check=rules.transcripts_na_check.output
    output:
        temp("logs/{accession}-quantile_normalise_transcripts_{metric}.done")
    params:
        transcripts="{accession}-transcripts-{metric}.tsv.undecorated" ,
        qntranscripts="{accession}-transcripts-{metric}.tsv.undecorated.quantile_normalized"
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"

        # transcripts_na_check rule done
        echo {input.transcripts_na_check}

        if [ -s {params.transcripts} ] ; then

            {workflow.basedir}/bin/quantile_normalize.sh -c {input.xml} -s {params.transcripts} -d {params.qntranscripts} -b {workflow.basedir}/bin

            if [ $? -ne 0 ]; then
                echo "ERROR: Failed to quantile normalize for {params.transcripts}  " >&2
                exit 1
            fi
        else
            echo "File {params.transcripts} not found "
        fi
        touch {output}
        """

rule summarize_transcripts:
    """
    Summarize transcript expression in tpms, if the file exists,
    either into median per biological replicate or into quartile per assay group.
    """
    conda: "envs/perl-atlas-modules.yaml"
    log: "logs/{accession}-summarize_transcripts_{metric}.log"
    input:
        xml="{accession}-configuration.xml",
        getqn=rules.quantile_normalise_transcripts.output
    output:
        temp("logs/{accession}-summarize_transcripts_{metric}.done")
    params:
        qn_transcripts="{accession}-transcripts-{metric}.tsv.undecorated.quantile_normalized",
        agg_transcripts="{accession}-transcripts-{metric}.tsv.undecorated.aggregated"
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"

        # transcripts qn rule done
        echo {input.getqn}

        if [ -s {params.qn_transcripts} ] ; then

            perl {workflow.basedir}/bin/gxa_summarize_expression.pl  \
                --configuration {input.xml}  \
                < {params.qn_transcripts}  \
                > {params.agg_transcripts}
            if [ $? -ne 0 ]; then
                echo "ERROR: Failed to summarize transcript expressions TPMs for {params.qn_transcripts}  " >&2
                exit 1
            fi
            # maintain previous behaviour - quantile_normalized files are temporary
            rm {params.qn_transcripts}
        else
            echo "File {params.qn_transcripts} not found. Transcript summary not performed "
        fi
        touch {output}
        """


rule generate_methods_baseline_rnaseq:
    """
    Fetches metadata about the analysis methods used in ISL/iRap to preprocess the experiment,
    to generate analysis methods.
    """
    conda: "envs/perl-atlas-modules.yaml"
    log: "logs/{accession}-generate_methods_baseline_rnaseq.log"
    params:
        organism=get_organism(),
        template=get_methods_template_baseline(),
        isl_dir=get_isl_dir()  
    output:
        methods=temp("{accession}-analysis-methods.tsv_baseline_rnaseq")
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"

        if [ ! -s {params.template} ] ; then
            echo "Methods template not found "
            exit 1
        fi

        source {workflow.basedir}/bin/reprocessing_routines.sh
        expIslDir={params.isl_dir}/{wildcards.accession}/{params.organism}
        echo "ISL dir: $expIslDir"

        [ ! -z $expIslDir+x ] || (echo "snakemake param exp_isl_dir needs to defined in rule" && exit 1)

        # IRAP methods version file
        if [ ! -s "$expIslDir/irap.versions.tsv" ] ; then
            echo "$expIslDir/irap.versions.tsv not found for {wildcards.accession} "
            exit 1
        fi

        # set env variables for mapper and quantification methods from used in irap.
        get_methods_from_irap "$expIslDir/irap.versions.tsv"
        [ ! -z ${{baseline_mapper+x}} ] || (echo "Env var baseline_mapper not defined." && exit 1)
        [ ! -z ${{baseline_quantMethod+x}} ] || (echo "Env var baseline_quantMethod not defined." && exit 1)
        [ ! -z ${{de_mapper+x}} ] || (echo "Env var de_mapper not defined." && exit 1)
        [ ! -z ${{de_quantMethod+x}} ] || (echo "Env var de_mapper not defined." && exit 1)
        [ ! -z ${{de_deMethod+x}} ] || (echo "Env var de_mapper not defined." && exit 1)

        echo $baseline_mapper
        echo $baseline_quantMethod
        echo $de_mapper
        echo $de_quantMethod
        echo $de_deMethod

        perl {workflow.basedir}/bin/gxa_generate_methods.pl "$expIslDir/irap.versions.tsv" {wildcards.accession} {params.organism} {params.template} "${{baseline_mapper:?}}" "${{baseline_quantMethod:?}}" "${{de_mapper:?}}" "${{de_quantMethod:?}}" "${{de_deMethod:?}}" > {output.methods}

        if [ $? -ne 0 ]; then
            echo "ERROR: Failed to generate analysis methods for {wildcards.accession}" >&2
            exit 1
        fi
        cp {output.methods} {wildcards.accession}-analysis-methods.tsv
        """

rule decorate_expression_baseline:
    """
    Decorate rna-seq baseline experiment with gene name from the latest Ensembl release.
    """
    container: "docker://quay.io/ebigxa/ensembl-update-env:amm1.1.2"
    log: "logs/{accession}-decorate_expression_baseline_{metric}.log"
    input:
        expression="{accession}-{metric}.tsv.undecorated.aggregated"
    params:
        organism=get_organism()
    output:
        decoexpression="{accession}-{metric}.tsv"
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        source {workflow.basedir}/bin/reprocessing_routines.sh

        geneNameFile=$( get_geneNameFile_given_organism {params.organism}  )

        echo $geneNameFile

        test -s "$geneNameFile" || (  >&2 echo "$0 gene name file not found: $geneNameFile" ; exit 1 )
        #test -s "$FPKMexpressionsFile" -o -s "$TPMexpressionsFile" || (  >&2 echo "$0 no data files for $e" ; exit 1 )

        decoratedFile={output.decoexpression} 

        # Ammonite REPL & Script-Runner

        amm -s {workflow.basedir}/bin/decorateFile.sc \
        --geneNameFile "$geneNameFile" \
        --source {input.expression} \
        | awk 'NR == 1; NR > 1 {{print $0 | "sort -n"}}' \
        > $decoratedFile.swp

        decoratedFileLength=$(wc -l "$decoratedFile.swp" | cut -f 1 -d ' ' )
        if [ -s "$decoratedFile.swp" ] && [ "$decoratedFileLength" -gt 1 ]; then
            mv $decoratedFile.swp $decoratedFile
        else
            echo "ERROR: decorate_rnaseq_file baseline for {wildcards.accession} and {wildcards.metric}"
            exit 1
        fi
        """


rule decorate_transcripts_baseline:
    """
    Decorate rna-seq baseline transcripts with transcript name from the latest Ensembl release.
    """
    container: "docker://quay.io/ebigxa/ensembl-update-env:amm1.1.2"
    log: "logs/{accession}-decorate_transcripts_baseline_{metric}.log"
    resources:
        mem_mb=get_mem_mb
    input:
        getagg=rules.summarize_transcripts.output
    params:
        organism=get_organism(),
        transcripts="{accession}-transcripts-{metric}.tsv.undecorated.aggregated",
        decotranscripts="{accession}-transcripts-{metric}.tsv"
    output:
        temp("logs/{accession}-transcripts-{metric}.tsv.done")
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        source {workflow.basedir}/bin/reprocessing_routines.sh

        # transcripts undecorated-aggregated rule done
        echo {input.getagg}

        if [ -s {params.transcripts} ] ; then

            geneNameFile=$( get_geneNameFile_given_organism {params.organism}  )
            transcriptFile=$( get_transcriptFile_given_organism {params.organism}  )

            echo $geneNameFile
            echo $transcriptFile

            test -s "$geneNameFile" || (  >&2 echo "$0 gene name file not found: $geneNameFile" ; exit 1 )
            test -s "$transcriptFile" || (  >&2 echo "$0 transcript to gene mapping file not found: $transcriptFile" ; exit 1 )

            decoratedFile={params.decotranscripts} 

            # Ammonite REPL & Script-Runner
            #export JAVA_OPTS="-Xmx3000M"
            amm -s {workflow.basedir}/bin/decorateFile.sc \
                --geneIdFile "$transcriptFile" \
                --geneNameFile "$geneNameFile" \
                --source {params.transcripts} \
                | awk 'NR == 1; NR > 1 {{print $0 | "sort -n"}}' \
                > $decoratedFile.swp

            decoratedFileLength=$(wc -l "$decoratedFile.swp" | cut -f 1 -d ' ' )
            if [ -s "$decoratedFile.swp" ] && [ "$decoratedFileLength" -gt 1 ]; then
                mv $decoratedFile.swp $decoratedFile
            else
                echo "ERROR: decorate_transcript_baseline_file for {wildcards.accession} and {wildcards.metric}"
	            exit 1
            fi
        else
            echo "File {params.transcripts} not found. Decorate transcript for baseline rna-seq summary not performed "
        fi
        touch {output}
        """


rule create_tracks_symlinks:
    """
    Create bedgraph tracks symlinks during reprocessing (only for tpms).
    """
    log: "logs/{accession}.{assay_id}.create_tracks_symlinks_{metric}.log"
    input:
        bedGraph="{accession}.{assay_id}.genes.expressions_{metric}.bedGraph"
    params:
        output="{accession}.{assay_id}.genes.expressions.bedGraph"
    output:
        temp("logs/{accession}.{assay_id}.create_tracks_symlinks_{metric}.done")
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"

        # remove if any old bedgraph files
        if [ -s {params.output} ] ; then
            rm {params.output} 
        fi

        ln -s {input.bedGraph} {params.output}
        touch {output}
        """



# differential_rnaseq_experiment

rule differential_statistics_rnaseq:
    """
    Calculate RNA-seq differential expression statistics and generate MvA plots.
    """
    conda: "envs/differential-stats.yaml"
    log: "logs/{accession}.differential_statistics_rnaseq.log"
    resources: mem_mb=get_mem_mb
    input:
        config_xml="{accession}-configuration.xml",
        raw_counts_undecorated="{accession}-raw-counts.tsv.undecorated"
    params:
        tmp_dir=get_tmp_dir()
    output:
        differential_expression="{accession}-analytics.tsv.undecorated",
        done=temp("logs/{accession}.differential_statistics_rnaseq.log.done")
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        PATH=$PATH:{workflow.basedir}/bin

        TMPDIR={params.tmp_dir}
        echo $TMPDIR"/tmp"
        if [ ! -d "$TMPDIR"/tmp ]; then
            mkdir $TMPDIR/tmp
        fi

        perl {workflow.basedir}/bin/diffAtlas_DE.pl --experiment {wildcards.accession} --directory ./
        if [ $? -ne 0 ]; then
	        echo "ERROR: Failed to generate differential expression statistics and MvA plots for {wildcards.accession}" >&2
	        rm -rf *.png {wildcards.accession}-analytics.tsv.undecorated
	        exit 1
        fi
        touch {output.done}
        """

rule check_mvaPlot_rnaseq:
    """
    This will check that all mvaPlots were produced by differential rna-seq analysis.
    """
    log: "logs/{accession}.{contrast_id}.check_mvaPlot_rnaseq.log"
    input:
        rules.differential_statistics_rnaseq.output.done
    output:
        temp("logs/{accession}-{contrast_id}-mvaPlot.png.done")
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"

        if [ -s "{wildcards.accession}-{wildcards.contrast_id}-mvaPlot.png" ] ; then
            echo "mvaPlot for {wildcards.accession} and {wildcards.contrast_id} successfully generated during differential rna-seq statistics: {wildcards.accession}-{wildcards.contrast_id}-mvaPlot.png"
        else
            echo "ERROR: mvaPlot for {wildcards.accession} and {wildcards.contrast_id} not produced during differential rna-seq statistics."
            exit 1
        fi
        touch {output}
        """

rule round_log2_fold_changes_rnaseq:
    """
    Round log2fold changes from differential expression analysis to one decimal place.
    It modifies the input and leaves the initial logs in .unrounded
    """
    conda: "envs/perl-math-round.yaml"
    log: "logs/{accession}.round_log2_fold_changes_rnaseq.log"
    resources: mem_mb=get_mem_mb
    input:
        "{accession}-analytics.tsv.undecorated"
    output:
        unrounded="{accession}-analytics.tsv.undecorated.unrounded" #rounded="{accession}-analytics.tsv.rounded"
    params:
        intermediate_rounded="{accession}-analytics.tsv.undecorated.rounded"
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"

        perl {workflow.basedir}/bin/round_log2_fold_changes.pl {input}

        if [ $? -ne 0 ]; then
	        echo "ERROR: Failed to round to one decimal place log2fold changes in {input} " >&2
	        rm -rf {params.intermediate_rounded}
	        exit 1
        fi
        mv {input} {output.unrounded}
        mv {params.intermediate_rounded} {input}
        """

rule generate_methods_differential_rnaseq:
    """
    Fetches metadata about the analysis methods used in ISL/iRap to preprocess the experiment,
    to generate analysis methods.
    """
    conda: "envs/perl-atlas-modules.yaml"
    log: "logs/{accession}-generate_methods_differential_rnaseq.log"
    params:
        organism=get_organism(),
        template=get_methods_template_differential(),
        isl_dir=get_isl_dir()  
    output:
        methods=temp("{accession}-analysis-methods.tsv_differential_rnaseq")
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"

        if [ ! -s {params.template} ] ; then
            echo "Methods template not found "
            exit 1
        fi

        source {workflow.basedir}/bin/reprocessing_routines.sh
        expIslDir={params.isl_dir}/{wildcards.accession}/{params.organism}
        echo "ISL dir: $expIslDir"

        [ ! -z $expIslDir+x ] || (echo "snakemake param exp_isl_dir needs to defined in rule" && exit 1)

        # IRAP methods version file
        if [ ! -s "$expIslDir/irap.versions.tsv" ] ; then
            echo "$expIslDir/irap.versions.tsv not found for {wildcards.accession} "
            exit 1
        fi

        # set env variables for mapper and quantification methods from used in irap.
        get_methods_from_irap "$expIslDir/irap.versions.tsv"
        [ ! -z ${{baseline_mapper+x}} ] || (echo "Env var baseline_mapper not defined." && exit 1)
        [ ! -z ${{baseline_quantMethod+x}} ] || (echo "Env var baseline_quantMethod not defined." && exit 1)
        [ ! -z ${{de_mapper+x}} ] || (echo "Env var de_mapper not defined." && exit 1)
        [ ! -z ${{de_quantMethod+x}} ] || (echo "Env var de_mapper not defined." && exit 1)
        [ ! -z ${{de_deMethod+x}} ] || (echo "Env var de_mapper not defined." && exit 1)

        echo $baseline_mapper
        echo $baseline_quantMethod
        echo $de_mapper
        echo $de_quantMethod
        echo $de_deMethod

        perl {workflow.basedir}/bin/gxa_generate_methods.pl "$expIslDir/irap.versions.tsv" {wildcards.accession} {params.organism} {params.template} "${{baseline_mapper:?}}" "${{baseline_quantMethod:?}}" "${{de_mapper:?}}" "${{de_quantMethod:?}}" "${{de_deMethod:?}}" > {output.methods}

        if [ $? -ne 0 ]; then
            echo "ERROR: Failed to generate analysis methods for {wildcards.accession}" >&2
            exit 1
        fi
        cp {output.methods} {wildcards.accession}-analysis-methods.tsv
        """





rule decorate_expression_differential_rnaseq:



# differential_microarray_experiment

rule get_normalized_expressions:

rule microarray_qc:

rule generate_methods_differential_microarray:

rule merge_probe_ids_per_gene:

rule microarray_calculate_analytics:

rule check_analytics_files:

rule round_log2_fold_changes_microarray:

rule decorate_expression_differential_microarray:





rule delete_experiment:

# mv data to atlas_exps
	

	
	
	




