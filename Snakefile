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

def get_zooma_exclusions():
    if 'zooma_exclusions' in config:
        return config['zooma_exclusions']
    else:
        return None

def get_isl_dir():
    if 'isl_dir' in config:
        return config['isl_dir']
    else:
        return None

def get_isl_genomes():
    if 'isl_genomes' in config:
        return config['isl_genomes']
    else:
        return None

def get_irap_versions():
    if 'irap_versions' in config:
        return config['irap_versions']
    else:
        return None

def get_irap_container():
    if 'irap_container' in config:
        return config['irap_container']
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

def read_run_deconv():
    if 'run_deconv_file' in config:
        global run_deconv
        with open(config['run_deconv_file'], 'r') as stream:
            try:
                run_deconv=yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                print(exc)
        return run_deconv

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
    """
    The logic is based on atlas production files. 
    """
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
    The logic is based on files processed by iRAP/ISL.
    """
    import os
    if 'metric' in config:
        return config['metric'].split(":")
    else:
        metrics_reprocess = []
        isl_dir = get_isl_dir()
        organism = get_organism()
        acc = config['accession']
        files_ils_dir = os.listdir( f"{isl_dir}/{acc}/{organism}" )
        # we only need one match, no need to traverse the full list
        if next((s for s in files_ils_dir if '.tpm.' in s), None):
            metrics_reprocess.append('tpms')
        if next((s for s in files_ils_dir if '.fpkm.' in s), None):
            metrics_reprocess.append('fpkms')

        if not metrics_reprocess:
            sys.exit("No metrics for reprocessing found in isl path.")
        else:
            return metrics_reprocess

def get_meta_config():
    if 'atlas_meta_config' in config:
        return config['atlas_meta_config']
    else:
        print(f" WARNING - atlas_meta_config not provided in config")
        return None

def get_db_params():
    """
    When goal is reprocess, return config parameters to establish connection with isl db.
    """
    isl_vars=[ 'oracle_home', 'python_user', 'python_connect_string', 'python_password' ] 
    db_params = []
    for i in isl_vars:
        if i in config:
            db_params.append( config[i] )
        else:
            print(f" Missing ISL db param: {i}")
            sys.exit(1)
    return db_params


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

def run_deconvolution(acc):
    if acc in run_deconv['deconv']['accession']:
        return True
    else:
        return False

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
    mem_avail = [ 2, 2, 4, 8, 16, 64, 128, 300 ]  
    if attempt > len(mem_avail):
        print(f"Attemps {attempt} exceeds the maximum number of attemps: {len(mem_avail)}")
        print(f"modify value of --restart-times or adjust mem_avail resources accordingly")
        sys.exit(1)
    else:
        return mem_avail[attempt-1] * 1000


def input_percentile_ranks(wildcards):
    """
    Return appropriate input for experiment type and analysis goal
    for rule percentile_ranks.
    """
    if config['goal'] == 'reprocess':
        if experiment_type=='rnaseq_mrna_differential':
            return [ f"logs/{wildcards['accession']}-decorate_differential_rnaseq.done" ]
        elif experiment_type=='proteomics_differential':
            return [ f"logs/{wildcards['accession']}-decorate_differential_proteomics.done" ]
        elif experiment_type == 'microarray_1colour_mrna_differential' or experiment_type =='microarray_2colour_mrna_differential' or experiment_type =='microarray_1colour_microrna_differential':
            inputs = []
            arr_designs=get_array_design_from_xml()
            for s in arr_designs:
                inputs.append( f"logs/{wildcards['accession']}_{s}-decorate_differential_microarray.done" )
            return inputs
        else:
            return None
    if config['goal'] == 'recalculations':
        # No input file needed - trick to force rule execution
        return f"{wildcards['accession']}.metadata_summary.yaml"

def input_differential_tracks_and_gsea(wildcards):
    """
    Return appropriate input for experiment type and analysis goal
    for rule differential_tracks and rule differential_gsea
    """
    if config['goal'] == 'reprocess':
        if experiment_type=='rnaseq_mrna_differential':
            return [ f"logs/{wildcards['accession']}-decorate_differential_rnaseq.done" ]
        elif experiment_type=='proteomics_differential':
            return [ f"logs/{wildcards['accession']}-decorate_differential_proteomics.done" ]
        elif experiment_type == 'microarray_1colour_mrna_differential' or experiment_type =='microarray_2colour_mrna_differential' or experiment_type =='microarray_1colour_microrna_differential':
            inputs = []
            arr_designs=get_array_design_from_xml()
            for s in arr_designs:
                inputs.append( f"logs/{wildcards['accession']}_{s}-decorate_differential_microarray.done" )
            return inputs
        else:
            return None
    if config['goal'] == 'recalculations':
        # No input file needed - trick to force rule execution
        return wildcards['accession']+'.metadata_summary.yaml'


def input_atlas_experiment_summary(wildcards):
    """
    Return appropriate input for experiment type and analysis goal
    for rule atlas_experiment_summary
    """
    if config['goal'] == 'reprocess':
        if experiment_type =='rnaseq_mrna_baseline': 
            return [ wildcards['accession']+'-raw-counts.tsv.undecorated', wildcards['accession']+'-analysis-methods.tsv_baseline_rnaseq' ]
        elif experiment_type=='rnaseq_mrna_differential': 
            return [ wildcards['accession']+'-raw-counts.tsv.undecorated', wildcards['accession']+'-analysis-methods.tsv_differential_rnaseq'  ]
        elif experiment_type == 'microarray_1colour_mrna_differential' or experiment_type =='microarray_1colour_microrna_differential':
            inputs = []
            arr_designs=get_array_design_from_xml()
            for s in arr_designs:
                inputs.append( f"logs/{wildcards['accession']}_{s}-decorate_differential_microarray.done" )
            inputs.append( f"{wildcards['accession']}-analysis-methods.tsv" )
            return inputs
        elif experiment_type =='microarray_2colour_mrna_differential':
            inputs = []
            arr_designs=get_array_design_from_xml()
            for s in arr_designs:
                inputs.append( f"{wildcards['accession']}_{s}-log-fold-changes.tsv.undecorated" )
                inputs.append( f"{wildcards['accession']}_{s}-average-intensities.tsv.undecorated" )
            inputs.append( f"{wildcards['accession']}-analysis-methods.tsv" )
            return inputs
        else:
            return None
    if config['goal'] == 'recalculations':
        return wildcards['accession']+'-configuration.xml'


def get_array_design_from_xml(wildcards):
    """
    Meant to be used within a rule. Parse xml config file here.
    """
    from xml.dom import minidom
    xmldoc = minidom.parse( wildcards['accession']+'-configuration.xml' )
    itemlist = xmldoc.getElementsByTagName('array_design')
    array_designs_grabbed = []
    for s in itemlist:
        array_designs_grabbed.append( " ".join(s.firstChild.nodeValue.split() ) )
    return array_designs_grabbed


def get_checkpoints_cp_atlas_exps(wildcards):
    """
    Only for reprocessing, to enable rule copy_experiment_from_analysis_to_atlas_exps
    """
    if config['goal'] == 'reprocess':
        inputs = get_outputs()
        inputs.remove( f"logs/{wildcards['accession']}-copy_experiment_from_analysis_to_atlas_exps.done" )
        inputs.remove( f"logs/{wildcards['accession']}-get_magetab_for_experiment.done" )
        # remove constraint on gsea.tsv and gsea_list.tsv,  as they could be removed by rule check_differential_gsea
        if experiment_type in ['rnaseq_mrna_differential', 'proteomics_differential', 'microarray_1colour_mrna_differential' , 'microarray_2colour_mrna_differential', 'microarray_1colour_microrna_differential']:
            inputs = [item for item in inputs if 'gsea.tsv' not in item]
            inputs = [item for item in inputs if 'gsea_list.tsv' not in item]
        return inputs
    else:
        return None

def input_round_log2_fold_changes(wildcards):
    """
    Ensure rename files has been run before rounding log2 fold changes, for differential proteomics
    """
    inputs_files = [ f"{wildcards['accession']}-analytics.tsv.undecorated" ]
    if experiment_type == 'proteomics_differential':
        inputs_files.append( f"logs/{wildcards['accession']}-rename_differential_proteomics_files.done" )
    return inputs_files


localrules: check_differential_gsea, link_baseline_coexpression, link_baseline_heatmap, create_tracks_symlinks, check_mvaPlot_rnaseq, check_normalized_expressions_microarray, delete_intermediate_files_microarray, touch_inputs_baseline

ruleorder: decorate_differential_rnaseq > decorate_differential_proteomics

wildcard_constraints:
    accession="E-\D+-\d+",
    metric="tpms|fpkms"


rule percentile_ranks:
    conda: "envs/atlas-internal.yaml"
    log: "logs/{accession}-percentile_ranks.log"
    resources: mem_mb=get_mem_mb
    input: input_percentile_ranks
    output:
        percentile_ranks_merged="{accession}-percentile-ranks.tsv"
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        rm -f {wildcards.accession}*-percentile-ranks.tsv
        for analytics in $(ls {wildcards.accession}*-analytics.tsv.unrounded); do
            {workflow.basedir}/atlas-analysis/calculate_percentile_ranks.R $analytics
        done
        # multiple ranks file will be generated for the microarray case,
        # in that case these need to be merged by gene id (first column),
        # each contrast on a different columns, all tab separated.
        # NA should be added for Gene x contrast for which there is no value
        percentile_ranks=( $(ls {wildcards.accession}*-percentile-ranks.tsv) )
        if [ ${{#percentile_ranks[@]}} -gt 1 ]; then
            # more than one file, requires merging by Gene.ID
            {workflow.basedir}/atlas-analysis/merge_by_gene_id.R {output.percentile_ranks_merged} {wildcards.accession}_*-percentile-ranks.tsv
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
        gff=get_gff(),
        check_point=input_differential_tracks_and_gsea
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
    input: input_differential_tracks_and_gsea
    params:
        organism=get_organism(),
        BIOENTITIES_PROPERTIES_PATH=config['bioentities_properties'],
        contrast_label=get_contrast_label,
        ext_db_label=get_ext_db_label,
        exp_type=get_from_config_or_metadata_summary('experiment_type')
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

        expType={params.exp_type}
        if [ "$expType" == "proteomics_differential" ]; then
            analyticsFile={wildcards.accession}-analytics.tsv
        else
            analyticsFile=$(grep -l -P "\\t{wildcards.contrast_id}\." {wildcards.accession}_A-*-analytics.tsv)
            if [ $? -ne 0 ]; then
                # rnaseq case
                analyticsFile={wildcards.accession}-analytics.tsv
            fi
        fi
        set -e
        annotationFile=$(find_properties_file_gsea {params.organism} {wildcards.ext_db})
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
        assay_label=get_assay_label,
        analytics="{accession}-{metric}.tsv"
    input:
        gff=get_gff(),
        analytics=lambda wildcards: f"{wildcards.accession}-{wildcards.metric}.tsv" if 'reprocess' in config['goal'] else f"{wildcards.accession}-{wildcards.metric}-touch_inputs_baseline.done"
    output:
        bedGraph="{accession}.{assay_id}.genes.expressions_{metric}.bedGraph"
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> {log:q}
        source {workflow.basedir}/bin/tracks_functions.sh
        echo "Past sourcing"
        generate_baseline_tracks {wildcards.accession} {wildcards.assay_id} {params.analytics} {input.gff} ./ {params.assay_label:q}
        """

rule touch_inputs_baseline:
    """
    Rule to avoid execution of upstream rules when forceall=true in baseline rna-seq recalculations.
    """
    log: 
        "logs/{accession}-{metric}-touch_inputs_baseline.log"
    output:
        temp("{accession}-{metric}-touch_inputs_baseline.done")
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"

        # Do not create the file if it does not exist (-c)
        touch -c -m -a {wildcards.accession}-{wildcards.metric}.tsv
        touch -c -m -a {wildcards.accession}-{wildcards.metric}.tsv.undecorated.aggregated

        touch {output}
        """

rule baseline_coexpression:
    conda: "envs/clusterseq.yaml"
    log: "logs/{accession}-{metric}-baseline_coexpression.log"
    resources: mem_mb=get_mem_mb
    params:
        num_retries=5,
        expression="{accession}-{metric}.tsv.undecorated.aggregated"
    input:  
        expression=lambda wildcards: f"{wildcards.accession}-{wildcards.metric}.tsv.undecorated.aggregated" if 'reprocess' in config['goal'] else f"{wildcards.accession}-{wildcards.metric}-touch_inputs_baseline.done"
    threads: 16
    output:
        coexpression_comp="{accession}-{metric}-coexpressions.tsv.gz"
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        {workflow.basedir}/atlas-analysis/run_coexpression_for_experiment.R {params.expression} {output.coexpression_comp} {workflow.basedir}/atlas-analysis/ {threads} {params.num_retries}
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
        expression=lambda wildcards: f"{wildcards.accession}-{wildcards.metric}.tsv" if 'reprocess' in config['goal'] else f"{wildcards.accession}-{wildcards.metric}-touch_inputs_baseline.done"
    output:
        heatmap="{accession}-heatmap-{metric}.pdf"
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        {workflow.basedir}/atlas-analysis/generateBaselineHeatmap.R --configuration {wildcards.accession}-configuration.xml \
		--input  {wildcards.accession}-{wildcards.metric}.tsv \
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
    """
    During (re)processing, some input files are necessary in addition to the XML and SDRF
    for internal-ExpressionAtlas, and it depends on the experiment type.
    """
    conda: "envs/atlas-internal.yaml"
    log: "logs/{accession}-atlas_experiment_summary.log"
    resources: mem_mb=get_mem_mb
    input:
        sdrf=get_sdrf(),
        input_files=input_atlas_experiment_summary
    output:
        rsummary="{accession}-atlasExperimentSummary.Rdata"
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        export SDRF_PATH={input.sdrf}
        {workflow.basedir}/atlas-analysis/atlasExperimentInR/createAtlasExperimentSummary.R \
	          --source ./ \
	          --accession {wildcards.accession} \
	          --output {output.rsummary}
        """


# rules below are specific for reprocessing

# baseline_rnaseq_experiment

rule add_runs_to_db:
    """
    Get run ids from config file and add them to isl db.
    """
    conda: "envs/isl-db.yaml"
    log: "logs/{accession}-add_runs_to_db.log"
    input:
        config_xml="{accession}-configuration.xml"
    params:
        db_params=get_db_params()
    output:
        temp("logs/{accession}-add_runs_to_db.done")
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        export TNS_ADMIN={params.db_params[0]}/network/admin
        export LD_LIBRARY_PATH={params.db_params[0]}/lib:$LD_LIBRARY_PATH
        export PATH={params.db_params[0]}/bin:$PATH
        export PYTHON_USER={params.db_params[1]}
        export PYTHON_CONNECT_STRING={params.db_params[2]}
        export PYTHON_PASSWORD={params.db_params[3]}

        python {workflow.basedir}/isl/db/scripts/get_run_ids_atlas_prod.py {input.config_xml}
        if [ $? -ne 0 ]; then
	        echo "ERROR: Failed to parse atlas config file and get run ids for {wildcards.accession} " >&2
	        exit 1
        fi
        touch {output}
        """

rule copy_raw_gene_counts_from_isl:
    """
    Copy raw gene counts file.
    """
    log: "logs/{accession}-copy_raw_gene_counts_from_isl.log"
    input:
        rules.add_runs_to_db.output
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
    QC step for rnaseq experiments.
    Non-standard experiments for which there is no QC information stored in the database
    should be added to 'skip_steps_file' to skip this rule.
    """
    conda: "envs/perl-atlas-modules.yaml"
    input: "{accession}-raw-counts.tsv.undecorated"
    log: "logs/{accession}-rnaseq_qc.log"
    output: "qc/{accession}-irap-single-lib-report.tsv"
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        echo "Running QC, raw counts file found {input}"

        {workflow.basedir}/atlas-analysis/qc/rnaseqQC.sh {wildcards.accession}
        qcExitCode=$?

        if [ "$qcExitCode" -eq 2 ]; then
            echo "Experiment {wildcards.accession} has been disqualified due to insufficient quality, exiting"
            exit 1
        elif [ "$qcExitCode" -ne 0 ]; then
            echo "ERROR: QC for {wildcards.accession} failed" >&2
            exit "$qcExitCode"
        fi

        if [ -d .qc ]; then
            # TODO: Temporary workaround until irap_single_lib is able to aggregate quality reports per study
            mv .qc qc
        fi
        """

rule quantile_normalise_expression:
    """
    Quantile normalize and summarize expression in tpms and fpkms.
    To maintain previous logic, the output file is marked as temporal
    """
    conda: "envs/quantile.yaml"
    log: "logs/{accession}-quantile_normalise_expression_{metric}.log"
    resources: mem_mb=get_mem_mb
    input:
        config_xml="{accession}-configuration.xml",
        expression="{accession}-{metric}.tsv.undecorated"
    output:
        qn_expression=temp("{accession}-{metric}.tsv.undecorated.quantile_normalized")
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        {workflow.basedir}/atlas-analysis/norm/quantile_normalize.sh  -c {input.config_xml} -s {input.expression} -d {output.qn_expression} -b {workflow.basedir}
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
    resources: mem_mb=get_mem_mb
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
        """

rule transcripts_na_check:
    """
    Replace NAs with 0 in Kallisto TPM transcripts, if the file exists
    (input and output file here is the same)
    """
    conda: "envs/atlas-internal.yaml"
    log: "logs/{accession}-rule-transcripts_na_check_{metric}.log"
    resources: mem_mb=get_mem_mb
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
            {workflow.basedir}/atlas-analysis/transcripts_expr_values_check.R {input.transcripts} $expIslDir/transcripts.raw.kallisto.tsv
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
    resources: mem_mb=get_mem_mb
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
            {workflow.basedir}/atlas-analysis/norm/quantile_normalize.sh -c {input.xml} -s {params.transcripts} -d {params.qntranscripts} -b {workflow.basedir}
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
    resources: mem_mb=get_mem_mb
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

            # maintain previous behaviour - quantile_normalized files are temporary
            rm {params.qn_transcripts}
        else
            echo "File {params.qn_transcripts} not found. Transcript summary not performed "
        fi
        touch {output}
        """

rule get_irap_versions_file:
    """
    If iRAP versions file not present, get it from the container.
    """
    log: 
        "logs/{accession}-get_irap_versions_file.log"
    params:
        irap_versions=get_irap_versions(),
        irap_container=get_irap_container()
    output:
        temp("logs/{accession}-get_irap_versions_file.done")
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"

        if [ ! -s {params.irap_versions} ] ; then
            echo "Attempting to transfer the file {params.irap_versions} from the singularity container {params.irap_container}"
            singularity exec {params.irap_container} cp /opt/irap/aux/mk/irap_versions.mk {params.irap_versions}
            if [ ! -s {params.irap_versions} ] ; then
                echo "ERROR: Failed to retrieve {params.irap_versions} from the container" >&2
                exit 1
            fi
        else
            echo "The file {params.irap_versions} already exists."
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
    input:
        rules.get_irap_versions_file.output
    params:
        organism=get_organism(),
        template=get_methods_template_baseline(),
        isl_dir=get_isl_dir(),
        isl_genomes=get_isl_genomes(),
        irap_versions=get_irap_versions()
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
        echo "ISL genome references: {params.isl_genomes}"

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

        # not used, only for differential rnaseq
        deseq2version='none'
        echo $deseq2version

        perl {workflow.basedir}/bin/gxa_generate_methods.pl "$expIslDir/irap.versions.tsv" {wildcards.accession} {params.organism} {params.template} "${{baseline_mapper:?}}" "${{baseline_quantMethod:?}}" "${{de_mapper:?}}" "${{de_quantMethod:?}}" "${{de_deMethod:?}}" "${{deseq2version:?}}" {params.isl_genomes} {params.irap_versions} > {output.methods}

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
    resources: mem_mb=get_mem_mb   
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
        source {workflow.basedir}/atlas-bash-util/generic_routines.sh

        geneNameFile=$( get_geneNameFile_given_organism {params.organism}  )

        echo $geneNameFile

        test -s "$geneNameFile" || (  >&2 echo "$0 gene name file not found: $geneNameFile" ; exit 1 )
   
        decoratedFile={output.decoexpression} 


        # pass avail custom memory to JVM for Ammonite REPL
        export JAVA_OPTS="-Xmx{resources.mem_mb}M"
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
        source {workflow.basedir}/atlas-bash-util/generic_routines.sh

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

            # pass avail custom memory to JVM for Ammonite REPL
            export JAVA_OPTS="-Xmx{resources.mem_mb}M"
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
    log: "logs/{accession}.{assay_id}.create_tracks_symlinks_tpms.log"
    input:
        bedGraph="{accession}.{assay_id}.genes.expressions_tpms.bedGraph"
    output: 
        "{accession}.{assay_id}.genes.expressions.bedGraph"
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"

        # remove if any old bedgraph files
        rm -f {output} 

        ln -s {input.bedGraph} {output}
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
        raw_counts_undecorated=lambda wildcards: f"{wildcards.accession}-raw-counts.tsv.undecorated" if experiment_type != 'proteomics_differential' else 'none_necessary'
    params:
        tmp_dir=get_tmp_dir(),
        exp_type=get_from_config_or_metadata_summary('experiment_type')
    output:
        differential_expression="{accession}-analytics.tsv.undecorated",
        done=temp("logs/{accession}.differential_statistics_rnaseq.done"),
        deseq2version=temp("{accession}.differential_statistics_rnaseq.deseq2version")
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"

        PATH=$PATH:{workflow.basedir}/atlas-analysis/differential
        
        expType={params.exp_type}
        if [ "$expType" == "proteomics_differential" ]; then
            # change only modification time of -analytics.tsv.undecorated; touch other temp outputs
            touch -m {output.differential_expression}
            touch {output.done}
            touch {output.deseq2version}
        else
            PATH=$PATH:{workflow.basedir}/atlas-analysis/differential
        
            source {workflow.basedir}/bin/reprocessing_routines.sh
            mktemp_dir {params.tmp_dir}

            rm -rf *.png {wildcards.accession}-analytics.tsv.undecorated

            perl {workflow.basedir}/atlas-analysis/differential/diffAtlas_DE.pl --experiment {wildcards.accession} --directory ./
        
            Rscript -e "library('DESeq2'); write.table(packageVersion('DESeq2'), file='{output.deseq2version}', quote=FALSE, col.names=FALSE, row.names = FALSE)" 

            touch {output.done}
        fi
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
    conda: "envs/atlas-internal.yaml"
    log: "logs/{accession}.round_log2_fold_changes_rnaseq.log"
    resources: mem_mb=get_mem_mb
    input:
        input_round_log2_fold_changes
    output:
        unrounded="{accession}-analytics.tsv.undecorated.unrounded",
        done=temp("logs/{accession}-round_log2_fold_changes_rnaseq.done")
    params:
        exp_type=get_from_config_or_metadata_summary('experiment_type'),
        intermediate_rounded="{accession}-analytics.tsv.undecorated.rounded"
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"

        rm -rf {params.intermediate_rounded}

        {workflow.basedir}/atlas-analysis/differential/round_log2_fold_changes.R \
            --experiment_type {params.exp_type} \
            --input_to_round {input[0]} \
            --intermediate_output {params.intermediate_rounded}  

        mv {input[0]} {output.unrounded}
        mv {params.intermediate_rounded} {input[0]}

        touch {output.done}
        """

rule generate_methods_differential_rnaseq:
    """
    Fetches metadata about the analysis methods used in ISL/iRap to preprocess the experiment,
    to generate analysis methods.
    """
    conda: "envs/perl-atlas-modules.yaml"
    log: "logs/{accession}-generate_methods_differential_rnaseq.log"
    input:
        deseq2version=rules.differential_statistics_rnaseq.output.deseq2version,
        check_irap_done=rules.get_irap_versions_file.output
    params:
        organism=get_organism(),
        template=get_methods_template_differential(),
        isl_dir=get_isl_dir(),
        isl_genomes=get_isl_genomes(),
        irap_versions=get_irap_versions()
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
        echo "ISL genome references: {params.isl_genomes}"

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

        deseq2version=`cat {input.deseq2version}`
        echo $deseq2version

        perl {workflow.basedir}/bin/gxa_generate_methods.pl "$expIslDir/irap.versions.tsv" {wildcards.accession} {params.organism} {params.template} "${{baseline_mapper:?}}" "${{baseline_quantMethod:?}}" "${{de_mapper:?}}" "${{de_quantMethod:?}}" "${{de_deMethod:?}}" "${{deseq2version:?}}" {params.isl_genomes} {params.irap_versions} > {output.methods}

        if [ $? -ne 0 ]; then
            echo "ERROR: Failed to generate analysis methods for {wildcards.accession}" >&2
            exit 1
        fi
        cp {output.methods} {wildcards.accession}-analysis-methods.tsv
        """

rule deconvolution:
    """
    Runs three deconvolution tools (DWLS, FARDEEP, EpiDIS) for experiments selected in deconvolution.yaml.
    """
    conda: "envs/deconvolution.yaml"
    log: "logs/{accession}.deconvolution.log"
    #resources: mem_mb=get_mem_mb
    resources: mem_mb=64000
    threads: 16
    input: 
        fpkms="{accession}-fpkms.tsv.undecorated",
        methods="{accession}-analysis-methods.tsv",
        sdrf=get_sdrf()
        #methods="{accession}-analysis-methods.tsv"
    params:
        organism=get_organism(),
        #exp_type=get_from_config_or_metadata_summary('experiment_type'),
        exp_type="rnaseq_mrna_differential",
        signature_dir=config["deconv_ref"]
    output:
        proportions="{accession}-summarized_proportions.tsv", 
        methods="{accession}-analysis-methods.updated.tsv",
        #results=temp(directory('Output/{accession}')),
        #splits=temp(directory('Tissue_splits/{accession}')),
        scratch=temp(directory('scratch/{accession}'))
    shell:
        """
        exec &> "logs/{wildcards.accession}.deconvolution.log"
        set -e # exit when any command fails
        echo "starting..."
        if [ ! -d "Tissue_splits/{wildcards.accession}" ]; then
            mkdir -p Tissue_splits/{wildcards.accession}
            Rscript {workflow.basedir}/atlas-analysis/deconvolution/splitAndScale.R {input.fpkms} {input.sdrf} {wildcards.accession}
        fi
        # list all files that FPKMs were split into
        files=$(ls Tissue_splits/{wildcards.accession}/{wildcards.accession}*-fpkms_scaled.rds)
        # iterate through tissues 
        for file in ${{files[@]}}; do
            # get tissue name from filename
            file=$(basename "$file")
            tissue="${{file#{wildcards.accession}-}}"
            tissue="${{tissue%-fpkms_scaled.rds}}"
            echo $tissue
            # search for suitable reference in deconvolution library
            REFERENCE_FOUND=$(Rscript {workflow.basedir}/atlas-analysis/deconvolution/findReference.R $tissue {params.signature_dir})
            # check if deconvolution for tissue was already completed
            if [ "$REFERENCE_FOUND" == "noref" ]; then
                echo "no reference for $tissue found"
                sc_reference_C1="noref"
            else
                # check if reference library is correct for this tissue
                number=$(ls {params.signature_dir}/${{REFERENCE_FOUND}}* | wc -l)
                if [ "$number" != 4 ]; then
                    echo "Error in reference library, check that there are no duplicated or missing references for $REFERENCE_FOUND!" 
                    exit 125
                fi 
                # find the different reference files
                sc_reference_C1=$(ls {params.signature_dir}/${{REFERENCE_FOUND}}_*_C1.rds | head -1)
                sc_reference_C0=$(ls {params.signature_dir}/${{REFERENCE_FOUND}}_*_C0_scaled.rds | head -1)
                sc_reference_phen=$(ls {params.signature_dir}/${{REFERENCE_FOUND}}_*_phenData.rds | head -1)
                # check if DWLS output already exists and results are more recent than the reference
                if [ ! -f "Output/{wildcards.accession}/{wildcards.accession}-${{tissue}}_res_DWLS.rds" ] || [ "Output/{wildcards.accession}/{wildcards.accession}-${{tissue}}_res_DWLS.rds" -ot "$sc_reference_C1" ]; then
                    echo "$REFERENCE_FOUND for $tissue found, running deconvolution"
                    # run deconvlution for this tisssue with FARDEEP, DWLS and EpiDISH
                    mkdir -p Output/{wildcards.accession}
                    {workflow.basedir}/scripts/run_deconvolution.sh $tissue {wildcards.accession} $sc_reference_C1 $sc_reference_C0 $sc_reference_phen {workflow.basedir}
                else
                    echo "$REFERENCE_FOUND for $tissue found, Skipping deconolution as for $tissue results already exist"
                fi
                mkdir -p ConsensusPlot/{wildcards.accession}
                Rscript {workflow.basedir}/scripts/getConsensus.R {wildcards.accession} $tissue
            fi
            # produce output files
            Rscript {workflow.basedir}/scripts/summarizeDeconvolutionResults.R {input.sdrf} {wildcards.accession} $tissue $sc_reference_C1 {output.proportions}
            Rscript {workflow.basedir}/scripts/getDeconvolutionInfo.R $tissue {wildcards.accession} $sc_reference_C1
        done
        Rscript {workflow.basedir}/scripts/appendAnalysisMethods.R {input.methods} {wildcards.accession}
        """

rule decorate_differential_rnaseq:
    """
    Decorate differential rna-seq experiment with gene name from the latest Ensembl release.
    """
    container: "docker://quay.io/ebigxa/ensembl-update-env:amm1.1.2"
    log: "logs/{accession}-decorate_differential_rnaseq.log"
    resources: mem_mb=get_mem_mb   
    input:
        "{accession}-analytics.tsv.undecorated",
        "{accession}-analytics.tsv.undecorated.unrounded",
        "{accession}-raw-counts.tsv.undecorated"
    params:
        organism=get_organism()
    output:
        "{accession}-analytics.tsv",
        "{accession}-analytics.tsv.unrounded",
        "{accession}-raw-counts.tsv",
        temp("logs/{accession}-decorate_differential_rnaseq.done")
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        source {workflow.basedir}/bin/reprocessing_routines.sh
        source {workflow.basedir}/atlas-bash-util/generic_routines.sh

        geneNameFile=$( get_geneNameFile_given_organism {params.organism}  )

        echo $geneNameFile

        test -s "$geneNameFile" || (  >&2 echo "$0 gene name file not found: $geneNameFile" ; exit 1 )

        for i in {input}
        do  
            echo "${{i}}"
            decoratedFile=`echo "${{i}}" | sed 's/\.undecorated//'`
            echo "$decoratedFile"

            # pass avail custom memory to JVM for Ammonite REPL
            export JAVA_OPTS="-Xmx{resources.mem_mb}M"
            amm -s {workflow.basedir}/bin/decorateFile.sc \
                --geneNameFile "$geneNameFile" \
                --source "${{i}}" \
                | awk 'NR == 1; NR > 1 {{print $0 | "sort -n"}}' \
                > $decoratedFile.swp

            decoratedFileLength=$(wc -l "$decoratedFile.swp" | cut -f 1 -d ' ' )
            if [ -s "$decoratedFile.swp" ] && [ "$decoratedFileLength" -gt 1 ]; then
                mv $decoratedFile.swp $decoratedFile
            else                                                                                                                                                                                                                    
                echo "ERROR: decorate_differential_rnaseq for {wildcards.accession} and undecorated input ${{i}} "
                exit 1                                                                                                                                                                                                              
            fi
        done
        touch {output[3]}                                                                                                                                                                                                       
        """                                                                                                                                                                                                                         
                                                                                                                                                                                                                  

# differential_microarray_experiment

rule get_normalized_expressions_microarray:
    """
    Get normalized expressions for differential microarray analysis.
    """
    conda: "envs/quantile.yaml"
    log: "logs/{accession}-get_normalized_expressions_microarray.log"
    resources: mem_mb=get_mem_mb
    params:
        tmp_dir=get_tmp_dir()
    output:
        temp("logs/{accession}-get_normalized_expressions_microarray.done")
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"

        source {workflow.basedir}/bin/reprocessing_routines.sh
        mktemp_dir {params.tmp_dir}

        ## get path to IDF file name and Array Express load directory
        # -i flag will retrieve path to idf filename 
        # -d flag will retrieve path to Array Express load directory
        # -m flag will retrieve path to mirbase directory
        idf_filename=$(perl {workflow.basedir}/bin/get_magetab_paths.pl -e {wildcards.accession} -i) 
        ae_dir=$(perl {workflow.basedir}/bin/get_magetab_paths.pl -e {wildcards.accession} -d)
        mirbase_dir=$(perl {workflow.basedir}/bin/get_magetab_paths.pl -e {wildcards.accession} -m)  

        echo $idf_filename
        echo $ae_dir
        echo $mirbase_dir

        # Get normalized expressions
        perl {workflow.basedir}/atlas-analysis/norm/arrayNormalization.pl {wildcards.accession} $idf_filename $ae_dir $mirbase_dir {workflow.basedir} $(pwd)

        touch {output}
        """ 

rule check_normalized_expressions_microarray:
    """
    Check that normalized expressions have been created for each array_design.
    """
    log: "logs/{accession}-check_normalized_expressions_microarray.log"
    input:
        rule_get_done=rules.get_normalized_expressions_microarray.output
    params:
        array_designs=get_array_design_from_xml
    output:
        temp("logs/{accession}-check_normalized_expressions_microarray.done")
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"

        for p in {params.array_designs}
        do  
            echo "Array designs found in the xml config file: ${{p}}"
            if [ -s {wildcards.accession}_"${{p}}-normalized-expressions.tsv.undecorated" ] ; then
                echo "File {wildcards.accession}_${{p}}-normalized-expressions.tsv.undecorated exists for {wildcards.accession} and array_design ${{p}}"
                cp {wildcards.accession}_"${{p}}-normalized-expressions.tsv.undecorated" {wildcards.accession}_"${{p}}-normalized-expressions.tsv.undecorated.unmerged"
            else
                echo "ERROR: File does not exist for {wildcards.accession} and array_design ${{p}}" >&2
                exit 1
            fi
        done
        touch {output}
        """



rule microarray_qc:
    """
    Run array quality control (and modify the experiment configuration file as necessary).
    """
    conda: "envs/quantile.yaml"
    log: "logs/{accession}-microarray_qc.log"
    input:
        rule_check_done=rules.check_normalized_expressions_microarray.output
    resources: mem_mb=get_mem_mb
    params:
        array_designs=get_array_design_from_xml
    output:
        temp("logs/{accession}-microarray_qc.done")
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"

        ## get path to IDF file name and Array Express load directory
        # -i flag will retrieve path to idf filename 
        # -d flag will retrieve path to Array Express load directory
        # -m flag will retrieve path to mirbase directory
        idf_filename=$(perl {workflow.basedir}/bin/get_magetab_paths.pl -e {wildcards.accession} -i) 
        ae_dir=$(perl {workflow.basedir}/bin/get_magetab_paths.pl -e {wildcards.accession} -d)
        mirbase_dir=$(perl {workflow.basedir}/bin/get_magetab_paths.pl -e {wildcards.accession} -m)  

        echo $idf_filename
        echo $ae_dir
        echo $mirbase_dir

        # arrayQualityMetrics should be > 3.32.0

        {workflow.basedir}/atlas-analysis/arrays/arrayQC.sh $(pwd) $idf_filename $ae_dir $mirbase_dir {workflow.basedir}/atlas-analysis/arrays
        qcExitCode=$?

        if [ "$qcExitCode" -eq 2 ]; then
	        echo "Experiment $expAcc has been disqualified due to insufficient quality, exiting"
            exit 0
        elif [ "$qcExitCode" -ne 0 ]; then
	        echo "ERROR: QC for {wildcards.accession} failed" >&2
            exit "$qcExitCode"
        fi

        # check that qc dirs exist for each array design
        for p in {params.array_designs}
        do
            if [ -d "./qc/{wildcards.accession}_${{p}}_QM" ]; then
                echo " ./qc/{wildcards.accession}_${{p}}_QM found"
            else
                echo "ERROR: ./qc/{wildcards.accession}_${{p}}_QM not found"
                exit 1
            fi
        done
        touch {output}
        """


rule generate_methods_differential_microarray:
    """
    Populate analysis methods for microarrays.
    """
    conda: "envs/perl-atlas-modules.yaml"
    log: "logs/{accession}-generate_methods_differential_microarray.log"
    input:
        config_xml="{accession}-configuration.xml"
    params:
        atlas_meta_config=get_meta_config(),
        exp_type=get_from_config_or_metadata_summary('experiment_type')
    output:
        "{accession}-analysis-methods.tsv"
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
	
        if [ -d "{params.atlas_meta_config}" ]; then
            export ATLAS_META_CONFIG={params.atlas_meta_config}
        else
            echo "ERROR: dir {params.atlas_meta_config} does not exist" >&2
            exit 1
        fi

        # Populate analysis methods
        target=""
        expType={params.exp_type}
        if [ "$expType" == "microarray_1colour_mrna_differential" ]; then
            arrayDataType=$(perl {workflow.basedir}/bin/get_experiment_info.pl --experiment {wildcards.accession} --xmlfile {input.config_xml} --rawdatafiles | head -1 | awk -F"." '{{print $2}}' | tr '[:upper:]' '[:lower:]')
            if [ "$arrayDataType" == "cel" ]; then
                target=../../affymetrix-differential-analytics-methods.tsv
            else
                target=../../onecolour-microarray-differential-analytics-methods.tsv
            fi
        elif [ "$expType" == "microarray_2colour_mrna_differential" ]; then
            target=../../twocolour-microarray-differential-analytics-methods.tsv
        elif [ "$expType" == "microarray_1colour_microrna_differential" ]; then
            target=../../onecolour-mirna-microarray-differential-analytics-methods.tsv
        elif [ "$expType" == "microarray_2colour_microrna_differential" ]; then
            target=../../twocolour-mirna-microarray-differential-analytics-methods.tsv
        else
            echo "ERROR: Unrecognised type: $expType for experiment: {wildcards.accession}" >&2
            exit 1
        fi

        if [ ! -f "$target" ]; then
            echo "ERROR: Failed to obtain analysis methods for {wildcards.accession}, file not found: $target" >&2
            exit 1
        fi
        ln -s $target {output}
        """



rule decorate_temp_norm_expr_microarray:
    """
    Generate temp decorated file.
    """
    container: "docker://quay.io/ebigxa/ensembl-update-env:amm1.1.2"
    log: "logs/{accession}-decorate_temp_norm_expr_microarray.log"
    resources: mem_mb=get_mem_mb
    input:
        rules.check_normalized_expressions_microarray.output
    output:
        temp("logs/{accession}-decorate_temp_norm_expr_microarray.done")
    params:
        organism=get_organism()
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"

        source {workflow.basedir}/bin/reprocessing_routines.sh
        source {workflow.basedir}/bin/decorate_microarray_routines.sh
        source {workflow.basedir}/atlas-bash-util/generic_routines.sh

        # pass avail custom memory to JVM for Ammonite REPL
        export JAVA_OPTS="-Xmx{resources.mem_mb}M"

        expPath=$(pwd)
        organism={params.organism}
        echo $expPath
        echo $organism

        export WF_BASEDIR={workflow.basedir}/bin
        export LC_ALL=C # avoid Perl warning

        # get normalized-expressions.tsv file
        find ${{expPath}} -maxdepth 1 -name "{wildcards.accession}_A-*-normalized-expressions.tsv.undecorated" \
            | xargs -n1 basename \
            | sed "s/{wildcards.accession}_//" \
            | sed "s/-normalized-expressions.tsv.undecorated//" \
            | while read -r arrayDesign ; do
            echo "Finding array design file..."
            arrayDesignFile=$(get_arraydesign_file ${{arrayDesign}} {params.organism}  )
            echo "...array design file found"
            # previously, however this has issues with organisms that share the array design
            # organism=$(get_organism_given_arraydesign_file ${{arrayDesignFile}} )
            if [ ! -z `echo $arrayDesignFile | grep mirbase` ]; then
                # This is a miRNA microarray experiment
                geneNameFile="${{expPath}}/mature.accession.tsv.aux"
                tail -n +2 ${{ATLAS_PROD}}/bioentity_properties/mirbase/${{organism}}.mature.tsv | awk -F"\t" '{{print $2"\t"$1}}' | sort -k 1,1 > $geneNameFile
            else
                geneNameFile=$(get_geneNameFile_given_organism {params.organism} )
            fi

            # generate temp decorated file normalized-expressions.tsv
            decorate_if_exists_norm "${{expPath}}/{wildcards.accession}_${{arrayDesign}}-normalized-expressions.tsv.undecorated" "$arrayDesignFile" "$geneNameFile"
        done
        touch {output}
        """

rule merge_probe_ids_microarray:
    """
    Merge probe ids with highest mean.
    """
    conda: "envs/atlas-internal.yaml"
    log: "logs/{accession}-merge_probe_ids_microarray.log"
    resources: mem_mb=get_mem_mb
    input:
        rules.decorate_temp_norm_expr_microarray.output
    params:
        array_designs=get_array_design_from_xml
    output:
        temp("logs/{accession}-merge_probe_ids_microarray.done")
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"

        # use the temp decorated file to find highest mean of probe ids per gene
        for p in {params.array_designs}
        do
            if [ -s "{wildcards.accession}_${{p}}-normalized-expressions.tsv.decorated.tmp" ]; then
                echo "Merging probe ids with highest mean per gene for {wildcards.accession} and array design ${{p}}"
                {workflow.basedir}/atlas-analysis/arrays/highestMeanProbeIdsPerGene.R "{wildcards.accession}_${{p}}-normalized-expressions.tsv.decorated.tmp"
            else
                echo "ERROR: {wildcards.accession}_${{p}}-normalized-expressions.tsv.decorated.tmp doesn't exist"
                exit 1
            fi
        done
        touch {output}
        """


rule differential_statistics_microarray:
    """
    Calculate differential expression statistics for microarrays - after probe ids have been merged.
    """
    conda: "envs/differential-stats.yaml"
    log: "logs/{accession}-differential_statistics_microarray.log"
    resources: mem_mb=get_mem_mb
    input:
        config_xml="{accession}-configuration.xml",
        merged_done=rules.merge_probe_ids_microarray.output
    params:
        tmp_dir=get_tmp_dir()
    output:
        done=temp("logs/{accession}-differential_statistics_microarray.done")
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        PATH=$PATH:{workflow.basedir}/atlas-analysis/differential
        
        source {workflow.basedir}/bin/reprocessing_routines.sh
        mktemp_dir {params.tmp_dir}

        # Calculate analytics
        rm -rf *.png *-analytics.tsv.undecorated
        perl {workflow.basedir}/atlas-analysis/differential/diffAtlas_DE.pl --experiment {wildcards.accession} --directory ./

        touch {output.done}
        """

rule check_nas_microarray:
    """
    Check that the analytics files have some p-values that are not NA.
    """
    conda: "envs/atlas-internal.yaml"
    log: "logs/{accession}-check_nas_microarray.log"
    input:
        rules.differential_statistics_microarray.output
    params:
        array_designs=get_array_design_from_xml
    output:
        temp("logs/{accession}-check_nas_microarray.done")
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"

        for p in {params.array_designs}
        do
            if [ -s "{wildcards.accession}_${{p}}-analytics.tsv.undecorated" ]; then
                {workflow.basedir}/atlas-analysis/differential/check_na_pvals.R "{wildcards.accession}_${{p}}-analytics.tsv.undecorated"
            else
                echo "ERROR: {wildcards.accession}_${{p}}-analytics.tsv.undecorated does not exist"
                exit 1
            fi
        done
        touch {output}
        """

rule round_log2_fold_changes_microarray:                                                                                                                                                                       
    """                                                                                                                                                                                                        
    Round log2fold changes to one decimal place.                                                                                                                                                               
    It modifies the input and leaves the initial logs in .unrounded                                                                                                                                            
    """                                                                                                                                                                                                        
    conda: "envs/atlas-internal.yaml"                                                                                                                                                                         
    log: "logs/{accession}_{array_design}-round_log2_fold_changes_microarray.log"                                                                                                                              
    resources: mem_mb=get_mem_mb                                                                                                                                                                               
    input:                                                                                                                                                                                                     
        rules.check_nas_microarray.output                                                                                                                                                                   
    output:                                                                                                                                                                                                    
        unrounded="{accession}_{array_design}-analytics.tsv.undecorated.unrounded",
        tmp=temp("logs/{accession}_{array_design}-round_log2_fold_changes_microarray.done")                                                                                                                            
    params:
        exp_type=get_from_config_or_metadata_summary('experiment_type'),                                                                                                                                                  
        analytics="{accession}_{array_design}-analytics.tsv.undecorated",                                                                                                                                      
        intermediate_rounded="{accession}_{array_design}-analytics.tsv.undecorated.rounded"                                                                                                                    
    shell:                                                                                                                                                                                                     
        """                                                                                                                                                                                                    
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set                                                                                                                       
        exec &> "{log}"                                                                                                                                                                                        

        rm -rf {params.intermediate_rounded}  
                                                                                                                                                                                                         
        {workflow.basedir}/atlas-analysis/differential/round_log2_fold_changes.R \
            --experiment_type {params.exp_type} \
            --input_to_round {params.analytics} \
            --intermediate_output {params.intermediate_rounded}                                                                                                                               
                                                                                                                                                                                                  
        mv {params.analytics} {output.unrounded}                                                                                                                                                                         
        mv {params.intermediate_rounded} {params.analytics} 
        touch {output.tmp}                                                                                                                                                 
        """                                                                                                                                                                                                    


rule decorate_differential_microarray:
    """
    Decorate a microarray experiment with gene name and identifier from the latest 
    Ensembl (or miRBase - as applicable) release.
    """
    container: "docker://quay.io/ebigxa/ensembl-update-env:amm1.1.2"
    log: "logs/{accession}_{array_design}-decorate_differential_microarray.log"
    resources: mem_mb=get_mem_mb    
    input:
        rules.round_log2_fold_changes_microarray.output.tmp       
    params:
        organism=get_organism()
    output:
        "{accession}_{array_design}-analytics.tsv",
        temp("logs/{accession}_{array_design}-decorate_differential_microarray.done")
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"

        source {workflow.basedir}/bin/reprocessing_routines.sh
        source {workflow.basedir}/bin/decorate_microarray_routines.sh
        source {workflow.basedir}/atlas-bash-util/generic_routines.sh

        # pass avail custom memory to JVM for Ammonite REPL
        export JAVA_OPTS="-Xmx{resources.mem_mb}M"

        expPath=$(pwd)
        organism={params.organism}
        echo $expPath
        echo $organism

        export WF_BASEDIR={workflow.basedir}/bin
        #export LC_ALL=C # avoid Perl warning

        # get normalized-expressions.tsv file
        find ${{expPath}} -maxdepth 1 -name "{wildcards.accession}_{wildcards.array_design}-analytics.tsv.undecorated" \
            | xargs -n1 basename \
            | sed "s/{wildcards.accession}_//" \
            | sed "s/-analytics.tsv.undecorated//" \
            | while read -r arrayDesign ; do
            echo "Finding array design file..."
            arrayDesignFile=$(get_arraydesign_file ${{arrayDesign}} {params.organism}  )
            echo "...array design file found"

            # previously, however this has issues with organisms that share the array design
            # organism=$(get_organism_given_arraydesign_file ${{arrayDesignFile}} )
            if [ ! -z `echo $arrayDesignFile | grep mirbase` ]; then
                # This is a miRNA microarray experiment
                geneNameFile="${{expPath}}/mature.accession.tsv.aux"
                tail -n +2 ${{ATLAS_PROD}}/bioentity_properties/mirbase/${{organism}}.mature.tsv | awk -F"\t" '{{print $2"\t"$1}}' | sort -k 1,1 > $geneNameFile
            else
                geneNameFile=$(get_geneNameFile_given_organism {params.organism} )
            fi
            # ${{arrayDesign}} == {wildcards.array_design}
            # proceed with decorations, *-analytics.tsv is mandatory output for this rule
            decorate_if_exists "${{expPath}}/{wildcards.accession}_${{arrayDesign}}-analytics.tsv.undecorated" "$arrayDesignFile" "$geneNameFile"
            decorate_if_exists "${{expPath}}/{wildcards.accession}_${{arrayDesign}}-normalized-expressions.tsv.undecorated" "$arrayDesignFile" "$geneNameFile"
            decorate_if_exists "${{expPath}}/{wildcards.accession}_${{arrayDesign}}-average-intensities.tsv.undecorated" "$arrayDesignFile" "$geneNameFile"
            decorate_if_exists "${{expPath}}/{wildcards.accession}_${{arrayDesign}}-log-fold-changes.tsv.undecorated" "$arrayDesignFile" "$geneNameFile"
            decorate_if_exists "${{expPath}}/{wildcards.accession}_${{arrayDesign}}-analytics.tsv.undecorated.unrounded" "$arrayDesignFile" "$geneNameFile"
        done
        touch {output}                                                                                                                                                                                                   
        """                                                                                                                                                                                                                         
         

rule delete_intermediate_files_microarray:
    """
    Delete intermediate files after microarray reprocessing. It runs after decoration
    """
    log: "logs/{accession}-delete_intermediate_files_microarray.log"
    input: input_differential_tracks_and_gsea
    params:
        array_designs=get_array_design_from_xml
    output:
        temp("logs/{accession}-delete_intermediate_files_microarray.done")
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"

        rm -rf intensities-*.ff
        rm -rf rma-*.ff

        [ -e mature.accession.tsv.aux ] && rm mature.accession.tsv.aux

        for p in {params.array_designs}
        do
            [ -e "{wildcards.accession}_${{p}}-normalized-expressions.tsv.undecorated.unmerged" ] && rm "{wildcards.accession}_${{p}}-normalized-expressions.tsv.undecorated.unmerged" 
            [ -e "{wildcards.accession}_${{p}}-normalized-expressions.tsv.decorated.tmp" ] && rm "{wildcards.accession}_${{p}}-normalized-expressions.tsv.decorated.tmp"

        done
        touch {output}   
        """   



######################################################
# Final reprocessing rules, for all experiments

rule copy_experiment_from_analysis_to_atlas_exps:
    """
    Copy data to atlas_exps. Its execution timing depends on experiment_type.
    NOTE: Use target_dir = config['atlas_exps'] for production
    """
    conda: "envs/perl-atlas-modules.yaml"
    log: "logs/{accession}-copy_experiment_from_analysis_to_atlas_exps.log"
    input: get_checkpoints_cp_atlas_exps
    params:
        target_dir=get_tmp_dir(),
        privacy_status_file=config['priv_stat_file']
    output:
        temp("logs/{accession}-copy_experiment_from_analysis_to_atlas_exps.done")
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        export ATLAS_EXPS={params.target_dir} #"/tmp" # {params.target_dir} for production
        source {workflow.basedir}/bin/reprocessing_routines.sh
        source {workflow.basedir}/atlas-bash-util/generic_routines.sh

        echo "Copying data to stage for: {wildcards.accession} to $ATLAS_EXPS"

        copy_experiment_from_analysis_to_atlas_exps {wildcards.accession} {params.privacy_status_file}

        echo "Copied data to stage"

        touch {output} 
        """

rule get_magetab_for_experiment:
    """
    Generate condensed SDRF with Zooma mappings for the experiment - in atlas_exps.
    NOTE: Use target_dir = config['atlas_exps'] for production
    """
    conda: "envs/perl-atlas-modules.yaml"
    log: "logs/{accession}-get_magetab_for_experiment.log"
    input: rules.copy_experiment_from_analysis_to_atlas_exps.output
    params:
        target_dir=get_tmp_dir(),
        exp_type=get_from_config_or_metadata_summary('experiment_type'),
        zooma_exclusions=get_zooma_exclusions()
    output:
        temp("logs/{accession}-get_magetab_for_experiment.done")
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        export ATLAS_EXPS={params.target_dir} #"/tmp"  # edit {params.target_dir} for production
        source {workflow.basedir}/bin/reprocessing_routines.sh
        source {workflow.basedir}/atlas-bash-util/generic_routines.sh

        echo "Retrieving magetab files for {wildcards.accession}"

        idf_filename=$(perl {workflow.basedir}/bin/get_magetab_paths.pl -e {wildcards.accession} -i) 

        get_magetab_for_experiment {wildcards.accession} {params.exp_type} {workflow.basedir} {params.zooma_exclusions} $idf_filename
        
        echo "Retrieved magetab files for {wildcards.accession}"

        touch {output} 
        """

