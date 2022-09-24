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


# Generate condensed SDRF with Zooma mappings for the experiment
get_magetab_for_experiment() {
    expAcc=$1
    if [ $? -ne 0 ]; then
        echo "ERROR: failed to get $expAcc experiment type from XML config. Cannot generate condensed SDRF."
        exit 1
    fi
    pushd ${ATLAS_EXPS}

    # Get the experiment type from the XML config.
    expType=$2 #`${ATLAS_PROD}/sw/atlasinstall_prod/atlasprod/db/scripts/get_experiment_type_from_xml.pl $expAcc/$expAcc-configuration.xml`
    scriptsDir=$3/bin
    zooma_exclusions_file=$4
    idf_file=$5
    fixesFileDir=$ATLAS_PROD/sw/atlasinstall_prod/atlasprod/experiment_metadata
    echo "Using Zooma exclusions file $zooma_exclusions_file to generate condensed SDRF."
    echo "Experiment type: $expType"
    # Now generate condensed sdrf containing ontology mappings from Zooma. This
    # will also copy IDF from ArrayExpress load directory (using "-i" option).
    # If this is a baseline experiment, pass the factors XML filename as well to ensure factors match in condensed SDRF.
    if [[ $expType == *baseline ]] || [[ $expType == *baseline_dia ]]; then

        $scriptsDir/condense_sdrf.pl -e $expAcc -f $expAcc/$expAcc-factors.xml -z -i -o $expAcc -x $zooma_exclusions_file -fi $idf_file
        if [ $? -ne 0 ]; then
            echo "ERROR: Failed to generate $expAcc/${expAcc}.condensed-sdrf.tsv with Zooma mappings, trying without..."
            $scriptsDir/condense_sdrf.pl -e $expAcc -f $expAcc/$expAcc-factors.xml -i -o $expAcc -fi $idf_file
        fi
        if [ $? -ne 0 ]; then
            echo "ERROR: Failed to generate $expAcc/${expAcc}.condensed-sdrf.tsv"
            return 1
        fi
    else

        $scriptsDir/condense_sdrf.pl -e $expAcc -z -i -o $expAcc -x $zooma_exclusions_file -fi $idf_file
        if [ $? -ne 0 ]; then
            echo "ERROR: Failed to generate $expAcc/${expAcc}.condensed-sdrf.tsv with Zooma mappings, trying without..."
            $scriptsDir/condense_sdrf.pl -e $expAcc -i -o $expAcc -fi $idf_file
        fi
        if [ $? -ne 0 ]; then
            echo "ERROR: Failed to generate $expAcc/${expAcc}.condensed-sdrf.tsv"
            return 1
        fi
    fi

    if [ ! -s "$expAcc/${expAcc}.condensed-sdrf.tsv" ]; then
	    echo "ERROR: Failed to generate $expAcc/${expAcc}.condensed-sdrf.tsv"
	    return 1
    fi

    applyAllFixesForExperiment $expAcc $fixesFileDir
    if [ $? -ne 0 ]; then
	    echo "ERROR: Applying fixes for experiment $e failed" >&2
        return 1
    fi

    rm -rf $expAcc/$expAcc-zoomifications-log.tsv
    popd
}

# copy_experiment_from_analysis_to_atlas_exps
# and related functions

# 7: see directory, see its contents, can write (fg_atlas only)
# 5: see directory, see its contents: public read-only directory
experiment_directory_permissions_from_peach_api_privacy() {
    expAcc=$1
    response=`peach_api_privacy_status $expAcc`

    case $response in
        *public*)
            echo 755
            ;;
        *private*)
            echo 750
            ;;
        ?)
            >&2 echo "experiment_directory_permissions_from_peach_api_privacy could not determine privacy status for $1, received: $response"
            return 1
            ;;
    esac
}

copy_experiment_from_analysis_to_atlas_exps(){
    expAcc=$1
    sourceDir=$(get_analysis_path_for_experiment_accession "$expAcc" )
    if [ ! -d "$sourceDir" ] ; then
        echo "copy_experiment_from_analysis_to_atlas_exps ERROR: Could not find in analysis directory: $expAcc" >&2
        exit 1
    fi
    mode=$(experiment_directory_permissions_from_peach_api_privacy "$expAcc" )

    if [ ! "$mode" ]; then
      echo "copy_experiment_from_analysis_to_atlas_exps ERROR: Failed to retrieve public/private status for $expAcc" >&2
      exit 1
    fi
    echo "source dir is: $sourceDir"
    copy_experiment -c "$mode" -s "$sourceDir" -t "${ATLAS_EXPS}/$expAcc"
    if [ $? -ne 0 ]; then
        echo "copy_experiment_from_analysis_to_atlas_exps ERROR: Command failed: copy_experiment -c $mode -s $sourceDir -t ${ATLAS_EXPS}/$expAcc" >&2
        exit 1
    fi

}

#Copy experiment into the target folder
# - archive previous content of the target in target_dir/archive
# - preserve timestamps (through rsync -a )
# - copy the subset of the data
copy_experiment() {
    rsyncExperimentFolders(){
    	rsync -a --copy-links --out-format="%n%L" \
    	    --exclude '*archive/**' \
    	    --exclude '*condensed-sdrf*' \
    	    --exclude '*lsf*' \
            --exclude 'logs' \
    	    --include '*/' \
    	    --include '*.tsv' \
    	    --include 'qc/**' \
    	    --include '*.xml' \
    	    --include '*.txt' \
    	    --include '*.png' \
    	    --include '*.bedGraph'\
    	    --include '*.Rdata' \
    	    --include '*.pdf' \
    	    --include '*.tsv.gz' \
    	    --exclude '*' \
    		$@
    }
    copy_experiment_usage(){
        echo "copy_experiment: " $@ " Usage: [-c=<mode> create if needed and set permissions] [-s=<source directory>] -t=<target directory>"
        return 2
    }

    mode=""
    source_dir=""
    target_dir=""
    local OPTARG OPTIND opt
    while getopts ":c:s:t:" opt; do
    	case $opt in
    		c)
          		mode=$OPTARG;
          		;;
    		s)
    			source_dir=$OPTARG;
    			if [ "$source_dir" -a ! -d "$source_dir" ]; then
    				copy_experiment_usage "-s should be a directory: $source_dir"
    			fi
    			;;
    		t)
    			target_dir=$OPTARG;
    			;;
    		?)
    			copy_experiment_usage "Unknown option: $OPTARG"
    			;;
    	esac
    done

    if [ ! "$target_dir" ] ; then
        copy_experiment_usage "Target directory not provided!"
    fi

    if [ "$mode" ]; then
    	mkdir -p "$target_dir"
    	chmod "$mode" "$target_dir"
    fi

    if [ ! -d "$target_dir" ] ; then
        copy_experiment_usage "Not a directory: $target_dir"
    fi

    if [ -d "$source_dir" ]; then
    	rsyncExperimentFolders --prune-empty-dirs -b --backup-dir "archive" --suffix ".1" "$source_dir/*" "$target_dir"
    fi
}

mktemp_dir() {
  TMPDIR=$1
  echo $TMPDIR"/tmp"
  if [ ! -d "$TMPDIR"/tmp ]; then
    mkdir $TMPDIR/tmp
  fi
}

