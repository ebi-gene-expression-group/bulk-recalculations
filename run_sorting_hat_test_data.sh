#!/usr/bin/env bash

GTF=$( pwd )/test-data/gff
DELETE_PREV_OUTPUT=True
BIOENTITIES_PROPERTIES=$( pwd )/test-data/bioentity_properties
SORTING_HAT=$( pwd )/Snakefile-sorting-hat
LOG_HANDLER=$( pwd )/log_handler.py
SN_CONDA_PREFIX=${SN_CONDA_PREFIX:-$( pwd )/conda_installs}

CONDA_PREFIX_LINE="--conda-prefix $SN_CONDA_PREFIX"
export LOG_PATH=${LOG_PATH:-$( pwd )/sorting.log}
USUAL_SM_ERR_OUT=${USUAL_SM_ERR_OUT:-$( pwd )/snakemake.log}

#create conda envs only
for yaml_file in $(ls $( pwd )/envs); do
  echo $yaml_file
  snakemake --dry-run --use-conda --conda-create-envs-only $CONDA_PREFIX_LINE  \
            --conda-frontend mamba \
            --config \
            yamlFile=$( pwd )/envs/$yaml_file \
            -s Snakefile-create-conda-envs
done

pushd $1
rm -f $LOG_PATH
touch $LOG_PATH
tail -f $LOG_PATH &

snakemake --use-conda --conda-frontend mamba \
        --log-handler-script $LOG_HANDLER \
        $CONDA_PREFIX_LINE \
        --keep-going \
        --config gtf_dir=$GTF \
        atlas_prod=path/to/atlasprod \
        atlas_exps=path/to/atlasexps \
        atlas_meta_config=path/to/supporting_files \
        sm_options="--use-conda --conda-frontend mamba --keep-going -j 2 $CONDA_PREFIX_LINE " \
        delete_previous_output=$DELETE_PREV_OUTPUT \  
        bioentities_properties=$BIOENTITIES_PROPERTIES -j 1 -s $SORTING_HAT &> $USUAL_SM_ERR_OUT

sleep 5
# remove tail process
kill $(jobs -p)
