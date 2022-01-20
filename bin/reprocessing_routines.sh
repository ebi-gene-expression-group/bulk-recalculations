# Old comment:
# ******** TEMPORARY **********
# Until we have "alternative views" in place and capability to access
# appropriate iRAP single-lib directory based on accession and organism
# name from experiment config XML, just select _all_ subdirectories of the
# iRAP single-lib directory for this accession.
# FIXME: Sanity check! If we do this there should only be ONE species
# directory in the iRAP single-lib directory! Fail if not.
# Update to comment: there are no multi-species experiments being processed for now
find_exp_isl_dir(){
    expAcc=$1
    result=$(find ${IRAP_SINGLE_LIB}/studies/$expAcc -mindepth 1 -maxdepth 1 -type d)

    count=$(wc -l <<< "$result")
    if [ "$count" -gt 1 ]; then
        echo "ERROR: Too many subdirectories in ${IRAP_SINGLE_LIB}/studies/$expAcc" >&2
        exit 1
    elif [ "$count" -eq 0 ]; then
        echo "Not found in ${IRAP_SINGLE_LIB}/studies : $expAcc" >&2
        exit 1
    else
        echo "$result"
        return 0
    fi
}

get_methods_from_irap(){
    irap_versions_file=$1
    if [ ! -s "$irap_versions_file" ]; then
        echo "irap.version file doesnt exist for $expAcc"
    fi

    ## Baseline
    # Reads alignment
    export baseline_mapper=$(cat $irap_versions_file | grep -P 'Reads alignment' | tail -1 | awk -F'\t' '{print $2}')

    # Gene, transcript and exon quantification
    export baseline_quantMethod=$(cat $irap_versions_file | grep -P 'Gene, transcript and exon quantification|Quantification' | awk -F'\t' '{print $2}')

    ## Differential
    # Reads alignment
    export de_mapper=$(cat $irap_versions_file | grep -P 'Reads alignment' | tail -1 | awk -F'\t' '{print $2}')

    # Gene  quantification
    export de_quantMethod=$(cat $irap_versions_file | grep -P 'Gene, transcript and exon quantification|Quantification' | awk -F'\t' '{print $2}')

    # Differential gene expression Method
    export de_deMethod=$(cat $irap_versions_file | grep -P 'Differential gene expression' | awk -F'\t' '{print $2}')
}

