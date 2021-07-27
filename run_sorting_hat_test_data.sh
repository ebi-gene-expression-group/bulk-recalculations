#!/usr/bin/env bash

GTF=$( pwd )/test-data/gff
BIOENTITIES_PROPERTIES=$( pwd )/test-data/bioentity_properties
SORTING_HAT=$( pwd )/Snakefile-sorting-hat
LOG_HANDLER=$( pwd )/log_handler.py

export LOG_PATH=${LOG_PATH:-$( pwd )/sorting.log}
pushd $1

snakemake --use-conda --conda-frontend mamba \
        --log-handler-script $LOG_HANDLER \
        --keep-going \
        --config gtf_dir=$GTF sm_options="--use-conda --conda-frontend mamba -j 2" \
        bioentities_properties=$BIOENTITIES_PROPERTIES -j 1 -s $SORTING_HAT &>/dev/null

cat $LOG_PATH
