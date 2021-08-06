#!/usr/bin/env bash

GTF=$( pwd )/test-data/gff
BIOENTITIES_PROPERTIES=$( pwd )/test-data/bioentity_properties
SORTING_HAT=$( pwd )/Snakefile-sorting-hat
LOG_HANDLER=$( pwd )/log_handler.py
CONDA_PREFIX=${CONDA_PREFIX:-$( pwd )/conda_installs}

CONDA_PREFIX_LINE="--conda-prefix $CONDA_PREFIX"
export LOG_PATH=${LOG_PATH:-$( pwd )/sorting.log}
USUAL_SM_ERR_OUT=${USUAL_SM_ERR_OUT:-$( pwd )/snakemake.log}
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
        bioentities_properties=$BIOENTITIES_PROPERTIES -j 1 -s $SORTING_HAT &> $USUAL_SM_ERR_OUT

sleep 5
# remove tail process
kill $(jobs -p)
