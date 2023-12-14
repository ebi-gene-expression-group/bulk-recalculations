#!/usr/bin/env bash

GOAL='recalculations'
# recalculations not currently implemented for proteomics experiments

CLUSTER=${CLUSTER:-"SLURM"}

NEXPS=${NEXPS:-30}
NJOBS=${NJOBS:-10}
GTF=$( pwd )/test-data/gff
FORCEALL=${FORCEALL:-true}
RESTART_TIMES=3
SKIP_STEPS=$( pwd )/step_skip.yaml
SKIP_ACC_REPROC=$( pwd )/accession_skip_reprocess.yaml
TEMPLATE_METHODS_BASELINE=$( pwd )/baseline_atlas_methods_template.conf
TEMPLATE_METHODS_DIFFERENTIAL=$( pwd )/differential_atlas_methods_template.conf
ZOOMA_EXCLUSIONS=$( pwd )/zooma_exclusions.yml
RUN_DECONVOLUTION=$( pwd )/accession_deconvolution.yaml
DECONV_REF=$atlas_prod/deconvolution_references/YYYY-MM-DD/

ISL_DIR=path/to/isl_dir
ISL_GENOMES_REFERENCES=$ISL_GENOMES
IRAP_VERSIONS=path/to/atlasprod/irap_versions.mk
IRAP_CONTAINER=path/to/singularity/irap_container:vx.y.z.sif
TMP_DEFINED=path/to/tmp_dir
TMP_DIR=${TMP_DEFINED:-$(mktemp -d)}
# CHECK_SPECIES file should come from a clone of the atlas-config repo in Jenkins:
CHECK_SPECIES=$( pwd )/atlas-species-name-mapping.yaml
BIOENTITIES_PROPERTIES=$( pwd )/test-data/bioentity_properties
SORTING_HAT=${SORTING_HAT:-$( pwd )/Snakefile-sorting-hat}
LOG_HANDLER=${LOG_HANDLER:-$( pwd )/log_handler.py}
SN_CONDA_PREFIX=${SN_CONDA_PREFIX:-$( pwd )/conda_installs}


if [ "$CLUSTER" = "LSF" ]; then
        PROFILE_LINE="--profile lsf-profilename"
        SLURM_OR_NONE=""
elif [ "$CLUSTER" = "SLURM" ]; then
        PROFILE_LINE="--profile slurm-profilename"
        SLURM_OR_NONE="--slurm"
else
        echo "'$CLUSTER' must be either LSF or SLURM"
        exit 1
fi
CLUSTER_CONFIG=${CLUSTER_CONFIG:-$( pwd )/cluster.yaml}


[ ! -z ${BIOSTUDIES_AE_PRIVACY_STATUS_FILE+x} ] || (echo "Env var BIOSTUDIES_AE_PRIVACY_STATUS_FILE not defined." && exit 1)

if [ ! -f "$BIOSTUDIES_AE_PRIVACY_STATUS_FILE" ]; then
        echo "$BIOSTUDIES_AE_PRIVACY_STATUS_FILE does not exist"
        exit 1
fi


CONDA_PREFIX_LINE="--conda-prefix $SN_CONDA_PREFIX"
export LOG_PATH=${LOG_PATH:-$( pwd )/sorting.log}
USUAL_SM_ERR_OUT=${USUAL_SM_ERR_OUT:-$( pwd )/snakemake.log}
# isl db
ORACLE_BASE=path/to/oracle_base
ORACLE_HOME=path/to/oracle_home
PYTHON_USER=
PYTHON_CONNECT_STRING=
PYTHON_PASSWORD=
PROT_MAGETABFILES=
# Optionally one, either ACCESSIONS or SPECIES can be defined
# otherwise ACCESSIONS supersedes SPECIES
# SPECIES=homo_sapiens:rattus_norvegicus
# ACCESSIONS=E-MTAB-2812
TOMCAT_HOST=
TOMCAT_HOST_USERNAME=
TOMCAT_HOST_PASSWORD=


###

if [ "$FORCEALL" = true ]; then FORCE_ALL="--forceall"; else FORCE_ALL=""; fi

if [ "$GOAL" != 'reprocess' ] && [ "$GOAL" != 'recalculations' ]; then
        echo "'$GOAL' is not a valid analysis for GXA"
        exit 1
fi

if [ -z ${ACCESSIONS+x} ]; then 
        ACC=""
else 
        ACC="accessions="$ACCESSIONS
fi

if [ -z ${SPECIES+x} ]; then 
        SPE=""
else 
        if [ -z ${ACCESSIONS+x} ]; then 
                SPE="species="$SPECIES
        else 
                SPE=""
        fi
fi



#create conda envs only
for yaml_file in $(ls $( pwd )/envs); do
  echo $yaml_file
  snakemake --use-conda --conda-create-envs-only $CONDA_PREFIX_LINE $PROFILE_LINE  \
            --conda-frontend mamba \
            --config \
            yamlFile=$( pwd )/envs/$yaml_file \
            -s Snakefile-create-conda-envs
done

pushd $1
rm -f $LOG_PATH
touch $LOG_PATH
tail -f $LOG_PATH &

start=`date +%s`
echo 'starting bulk '$GOAL'...'

snakemake $SLURM_OR_NONE --use-conda --conda-frontend mamba \
        --log-handler-script $LOG_HANDLER \
        $PROFILE_LINE \
        $FORCE_ALL \
        $CONDA_PREFIX_LINE \
        --latency-wait 10 \
        --keep-going \
        --config $ACC $SPE gtf_dir=$GTF \
        atlas_prod=path/to/atlasprod \
        atlas_exps=path/to/atlasexps \
        cluster_config=$CLUSTER_CONFIG \
        deconv_ref=$DECONV_REF \
        goal=$GOAL \
        atlas_meta_config=path/to/supporting_files \
        skip_steps_file=$SKIP_STEPS \
        skip_accessions_reproc_file=$SKIP_ACC_REPROC \
        run_deconv_file=$RUN_DECONVOLUTION \
        methods_base=$TEMPLATE_METHODS_BASELINE \
        methods_dif=$TEMPLATE_METHODS_DIFFERENTIAL \
        zooma_exclusions=$ZOOMA_EXCLUSIONS \
        isl_genomes=$ISL_GENOMES_REFERENCES \
        isl_dir=$ISL_DIR \
        irap_versions=$IRAP_VERSIONS \
        irap_container=$IRAP_CONTAINER \
        tmp_dir=$TMP_DIR \
        priv_stat_file=$BIOSTUDIES_AE_PRIVACY_STATUS_FILE \
        check_sp_file=$CHECK_SPECIES \
        oracle_home=$ORACLE_HOME \                                                                                                                        
        python_user=$PYTHON_USER \                                                                                                                        
        python_connect_string=$PYTHON_CONNECT_STRING \                                                                                                    
        python_password=$PYTHON_PASSWORD \
        tomcat_host_username=$TOMCAT_HOST_USERNAME \     
        tomcat_host_password=$TOMCAT_HOST_PASSWORD \
        tomcat_host=$TOMCAT_HOST \
        prot_magetabfiles=$PROT_MAGETABFILES \
        sm_options="$SLURM_OR_NONE --use-conda --conda-frontend mamba --keep-going $PROFILE_LINE -j $NJOBS $CONDA_PREFIX_LINE $FORCE_ALL --restart-times $RESTART_TIMES " \
        bioentities_properties=$BIOENTITIES_PROPERTIES -j $NEXPS -s $SORTING_HAT &> $USUAL_SM_ERR_OUT

end=`date +%s`
echo "bulk '$GOAL' took: "`expr $end - $start`" s"

sleep 5
# remove tail process
kill $(jobs -p)
