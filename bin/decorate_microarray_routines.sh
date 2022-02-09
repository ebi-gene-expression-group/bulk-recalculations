#!/bin/bash

# /bioentity_annotations/decorate_routines.sh
decorate_microarray_file() {
    f=$1
    arrayDesignFile=$2
    geneNameFile=$3

    decoratedFile=`echo $f | sed 's/\.undecorated//'`
    amm -s "$WF_BASEDIR"/decorateFile.sc \
        --geneIdFile "$arrayDesignFile" \
        --geneNameFile "$geneNameFile" \
        --source "$f" \
        | awk 'NR == 1; NR > 1 {print $0 | "sort -n"}' \
        > $decoratedFile.swp
    decoratedFileLength=$(wc -l "$decoratedFile.swp" | cut -f 1 -d ' ' )
    if [ -s "$decoratedFile.swp" ] && [ "$decoratedFileLength" -gt 1 ]; then
        mv $decoratedFile.swp $decoratedFile
        return 0
    else
      echo "ERROR: decorate_microarray_file $@"
	    return 1
    fi
}

# microarray-import/bin/mergeProbeIdsPerGene.sh
decorate_if_exists() {
    fileToDecorate=$1
    if [ -e "$fileToDecorate" ]; then
        decorate_normalized_file $@
        if [ $? -ne 0 ]; then
                echo "ERROR: FAILED decorate_normalized_microarray_file $@" >&2
                exit 1
        fi
    fi
}

# microarray-import/bin/mergeProbeIdsPerGene.sh
decorate_normalized_file() {
    f=$1
    arrayDesignFile=$2
    geneNameFile=$3

    echo decorate_microarray_file "$@"
    decoratedFile=`echo $f | sed 's/\.undecorated//'`
    amm -s "$WF_BASEDIR"/decorateFile.sc \
        --geneIdFile "$arrayDesignFile" \
        --geneNameFile "$geneNameFile" \
        --source "$f" \
        | perl -e 'print scalar <>, sort <>;' \
        > $decoratedFile.decorated.tmp
    decoratedFileLength=$(wc -l "$decoratedFile.decorated.tmp" | cut -f 1 -d ' ' )
    if [ -s "$decoratedFile.decorated.tmp" ] && [ "$decoratedFileLength" -gt 1 ]; then
        return 0
    else
     return 1
    fi
}

# /bash_util/generic_routines.sh
get_arraydesign_file() {
  arraydesign=$1
  organism=$2
  if [ -z ${2+x} ]; then
    find -L ${ATLAS_PROD}/bioentity_properties/array_designs -type f -name "*.${arraydesign}.tsv" | head -n1
  else
    find -L ${ATLAS_PROD}/bioentity_properties/array_designs -type f -name "${organism}.${arraydesign}.tsv" | head -n1
  fi
}


