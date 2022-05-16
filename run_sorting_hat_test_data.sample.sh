#!/usr/bin/env bash

GOAL='recalculations'

if [ "$GOAL" != 'reprocess' ] && [ "$GOAL" != 'recalculations' ]; then
        echo "'$GOAL' is not a valid analysis for GXA"
        exit 1
fi

# only one, either ACCESSIONS or SPECIES can be defined
# SPECIES=homo_sapiens:rattus_norvegicus

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


NEXPS=${NEXPS:-30}
NJOBS=${NJOBS:-10}
GTF=$( pwd )/test-data/gff
FORCEALL=${FORCEALL:-true}
RESTART_TIMES=3
SKIP_STEPS=$( pwd )/step_skip.yaml
TEMPLATE_METHODS_BASELINE=$( pwd )/baseline_atlas_methods_template.conf
TEMPLATE_METHODS_DIFFERENTIAL=$( pwd )/differential_atlas_methods_template.conf
ZOOMA_EXCLUSIONS=$( pwd )/zooma_exclusions.yml
ISL_DIR=path/to/isl_dir
ISL_GENOMES_REFERENCES=$ISL_GENOMES
IRAP_VERSIONS=path/to/atlasprod/irap_versions.mk
IRAP_CONTAINER=path/to/singularity/*.sif
TMP_DIR=path/to/tmp_dir
# CHECK_SPECIES file should come from a clone of the atlas-config repo in Jenkins:
CHECK_SPECIES=$( pwd )/atlas-species-name-mapping.yaml
BIOENTITIES_PROPERTIES=$( pwd )/test-data/bioentity_properties
SORTING_HAT=${SORTING_HAT:-$( pwd )/Snakefile-sorting-hat}
LOG_HANDLER=${LOG_HANDLER:-$( pwd )/log_handler.py}
SN_CONDA_PREFIX=${SN_CONDA_PREFIX:-$( pwd )/conda_installs}
PROFILE_LINE="--profile profilename"
LSF_CONFIG=${LSF_CONFIG:-$( pwd )/lsf.yaml}

if [ "$FORCEALL" = true ]; then FORCE_ALL="--forceall"; else FORCE_ALL=""; fi

CONDA_PREFIX_LINE="--conda-prefix $SN_CONDA_PREFIX"
export LOG_PATH=${LOG_PATH:-$( pwd )/sorting.log}
USUAL_SM_ERR_OUT=${USUAL_SM_ERR_OUT:-$( pwd )/snakemake.log}

# isl db
ORACLE_BASE=
ORACLE_HOME=
PYTHON_USER=
PYTHON_CONNECT_STRING=
PYTHON_PASSWORD=

#create conda envs only
for yaml_file in $(ls $( pwd )/envs); do
  echo $yaml_file
  snakemake --use-conda --conda-create-envs-only $CONDA_PREFIX_LINE  \
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
echo 'starting bulk-recalculations...'

snakemake --use-conda --conda-frontend mamba \
        --log-handler-script $LOG_HANDLER \
        $PROFILE_LINE \
        $CONDA_PREFIX_LINE \
        --latency-wait 10 \
        --keep-going \
        --config $ACC $SPE gtf_dir=$GTF \
        atlas_prod=path/to/atlasprod \
        atlas_exps=path/to/atlasexps \
        lsf_config=$LSF_CONFIG \
        goal=$GOAL \
        atlas_meta_config=path/to/supporting_files \
        skip_steps_file=$SKIP_STEPS \
        methods_base=$TEMPLATE_METHODS_BASELINE \
        methods_dif=$TEMPLATE_METHODS_DIFFERENTIAL \
        zooma_exclusions=$ZOOMA_EXCLUSIONS \
        isl_genomes=$ISL_GENOMES_REFERENCES \
        irap_versions=$IRAP_VERSIONS \
        irap_container=$IRAP_CONTAINER \
        check_sp_file=$CHECK_SPECIES \
        oracle_home=$ORACLE_HOME \                                                                                                                        
        python_user=$PYTHON_USER \                                                                                                                        
        python_connect_string=$PYTHON_CONNECT_STRING \                                                                                                    
        python_password=$PYTHON_PASSWORD \     
        sm_options="--use-conda --conda-frontend mamba --keep-going $PROFILE_LINE -j $NJOBS $CONDA_PREFIX_LINE $FORCE_ALL --restart-times $RESTART_TIMES " \
        bioentities_properties=$BIOENTITIES_PROPERTIES -j $NEXPS -s $SORTING_HAT &> $USUAL_SM_ERR_OUT

end=`date +%s`
echo "bulk-recalculations took: "`expr $end - $start`" s"

sleep 5
# remove tail process
kill $(jobs -p)
