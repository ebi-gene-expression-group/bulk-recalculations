#!/bin/bash

get_methods_from_irap(){
    irap_versions_file=$1
    if [ ! -s "$irap_versions_file" ]; then
        echo "irap.version file does not exist"
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

# to decorate rnaseq baseline expression:

get_geneNameFile_given_organism() {
  find_properties_file $1 "symbol"
}

find_properties_file() {
    organism=$1
    property=$2

    ensFile="${ATLAS_PROD}/bioentity_properties/ensembl/${organism}.ensgene.${property}.tsv"
    if [ -s "$ensFile" ]; then
        echo $ensFile
    else
        wbpsFile="${ATLAS_PROD}/bioentity_properties/wbps/${organism}.wbpsgene.${property}.tsv"
        if [ -s "$wbpsFile" ]; then
            echo $wbpsFile
        else
            >&2 echo "No annotation file found for organism $organism and property $property"
            exit 1
        fi
    fi
}

get_transcriptFile_given_organism() {
    organism=$1

    ensFile="${ATLAS_PROD}/bioentity_properties/ensembl/${organism}.ensgene.enstranscript.tsv"
    if [ -s "$ensFile" ]; then
        echo $ensFile
    else
        wbpsFile="${ATLAS_PROD}/bioentity_properties/wbps/${organism}.wbpsgene.wbpstranscript.tsv"
        if [ -s "$wbpsFile" ]; then
            echo $wbpsFile
        else
            >&2 echo "No transcript annotation file found for organism $organism"
            exit 1
        fi
    fi
}


