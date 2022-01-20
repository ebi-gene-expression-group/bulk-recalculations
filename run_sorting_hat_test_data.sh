#!/usr/bin/env bash

if [ -z ${ACCESSIONS+x} ]; then 
        ACC=""
else 
        ACC="accessions="$ACCESSIONS
fi

NEXPS=${NEXPS:-30}
NJOBS=${NJOBS:-10}
GTF=$( pwd )/test-data/gff
FORCEALL=${FORCEALL:-true}
RESTART_TIMES=3
SKIP_STEPS=$( pwd )/step_skip.yaml
TEMPLATE_METHODS_BASELINE=$( pwd )/baseline_atlas_methods_template.conf
TEMPLATE_METHODS_DIFFERENTIAL=$( pwd )/differential_atlas_methods_template.conf
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
        --config $ACC gtf_dir=$GTF \
        atlas_prod=path/to/atlasprod \
        atlas_exps=path/to/atlasexps \
        lsf_config=$LSF_CONFIG \
        goal='recalculations' \
        atlas_meta_config=path/to/supporting_files \
        skip_steps_file=$SKIP_STEPS \
        methods_base=$TEMPLATE_METHODS_BASELINE \
        methods_dif=$TEMPLATE_METHODS_DIFFERENTIAL \
        check_sp_file=$CHECK_SPECIES \
        sm_options="--use-conda --conda-frontend mamba --keep-going $PROFILE_LINE -j $NJOBS $CONDA_PREFIX_LINE $FORCE_ALL --restart-times $RESTART_TIMES " \
        bioentities_properties=$BIOENTITIES_PROPERTIES -j $NEXPS -s $SORTING_HAT &> $USUAL_SM_ERR_OUT

end=`date +%s`
echo "bulk-recalculations took: "`expr $end - $start`" s"

sleep 5
# remove tail process
kill $(jobs -p)
