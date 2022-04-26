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





# Generate condensed SDRF with Zooma mappings for the experiment
get_magetab_for_experiment() {
    expAcc=$1
    pushd  ${ATLAS_EXPS}

    # Get the experiment type from the XML config.
    expType=`${ATLAS_PROD}/sw/atlasinstall_prod/atlasprod/db/scripts/get_experiment_type_from_xml.pl $expAcc/$expAcc-configuration.xml`
    if [ $? -ne 0 ]; then
        echo "ERROR: failed to get $expAcc experiment type from XML config. Cannot generate condensed SDRF."
        exit 1
    fi

    # Now generate condensed sdrf containing ontology mappings from Zooma. This
    # will also copy IDF from ArrayExpress load directory (using "-i" option).
    # If this is a baseline experiment, pass the factors XML filename as well to ensure factors match in condensed SDRF.
    if [[ $expType == *baseline ]]; then

        ${ATLAS_PROD}/sw/atlasinstall_prod/atlasprod/experiment_metadata/condense_sdrf.pl -e $expAcc -f $expAcc/$expAcc-factors.xml -z -i -o $expAcc
        if [ $? -ne 0 ]; then
            echo "ERROR: Failed to generate $expAcc/${expAcc}.condensed-sdrf.tsv with Zooma mappings, trying without..."
            ${ATLAS_PROD}/sw/atlasinstall_prod/atlasprod/experiment_metadata/condense_sdrf.pl -e $expAcc -f $expAcc/$expAcc-factors.xml -i -o $expAcc
        fi
        if [ $? -ne 0 ]; then
            echo "ERROR: Failed to generate $expAcc/${expAcc}.condensed-sdrf.tsv"
            return 1
        fi
    else

        ${ATLAS_PROD}/sw/atlasinstall_prod/atlasprod/experiment_metadata/condense_sdrf.pl -e $expAcc -z -i -o $expAcc
        if [ $? -ne 0 ]; then
            echo "ERROR: Failed to generate $expAcc/${expAcc}.condensed-sdrf.tsv with Zooma mappings, trying without..."
            ${ATLAS_PROD}/sw/atlasinstall_prod/atlasprod/experiment_metadata/condense_sdrf.pl -e $expAcc -i -o $expAcc
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

    applyAllFixesForExperiment $expAcc
    if [ $? -ne 0 ]; then
	echo "ERROR: Applying fixes for experiment $e failed" >&2
	return 1
    fi

    rm -rf $expAcc/$expAcc-zoomifications-log.tsv
    popd
}

# ../bash_util/generic_routines.sh
applyAllFixesForExperiment() {
   exp=$1
   echo "Applying fixes for $exp ..."
    # Apply factor type fixes in idf file
    applyFixes $exp automatic_fixes_properties.txt idf.txt
    if [ $? -ne 0 ]; then
	echo "ERROR: Applying factor type fixes in idf file for $exp failed" >&2
	return 1
    fi

    # Commenting out SDRF bit as should not be needed any more.
    # Apply factor/sample characteristic type fixes to sdrf
    #applyFixes $exp automatic_fixes_properties.txt sdrf.txt
    #if [ $? -ne 0 ]; then
	#echo "ERROR: Applying sample characteristic/factor types fixes in sdrf file for $exp failed" >&2
	#return 1
    #fi
    # Apply sample characteristic/factor value fixes in sdrf file
    #applyFixes $exp automatic_fixes_values.txt sdrf.txt
    #if [ $? -ne 0 ]; then
	#echo "ERROR: Applying sample characteristic/factor value fixes in sdrf file for $exp failed" >&2
	#return 1
    #fi

    # Apply factor/sample characteristic type fixes to the condensed-sdrf file
    applyFixes $exp automatic_fixes_properties.txt condensed-sdrf.tsv
    if [ $? -ne 0 ]; then
	echo "ERROR: Applying sample characteristic/factor types fixes in sdrf file for $exp failed" >&2
	return 1
    fi
    # Apply sample characteristic/factor value fixes to the condensed-sdrf file
    applyFixes $exp automatic_fixes_values.txt condensed-sdrf.tsv
    if [ $? -ne 0 ]; then
	echo "ERROR: Applying sample characteristic/factor value fixes in sdrf file for $exp failed" >&2
	return 1
    fi
}

# ../bash_util/generic_routines.sh
# Applies fixes encoded in $fixesFile to $exp.$fileTypeToBeFixed.txt
applyFixes() {
    exp=$1
    fixesFile=$2
    fileTypeToBeFixed=$3
    atlasEnv=`atlas_env` 

    # Apply factor type fixes in ${fileTypeToBeFixed} file
    for l in $(cat $ATLAS_PROD/sw/atlasinstall_${atlasEnv}/atlasprod/experiment_metadata/$fixesFile | sed 's|[[:space:]]*$||g');
    do
	if [ ! -s "$exp/$exp.${fileTypeToBeFixed}" ]; then
	    echo "ERROR: $exp/$exp.${fileTypeToBeFixed} not found or is empty" >&2
	    return 1
	fi
	echo $l | grep -P '\t' > /dev/null
	if [ $? -ne 0 ]; then
	    echo  "WARNING: line: '$l' in automatic_fixes_properties.txt is missing a tab character - not applying the fix "
	fi
	correct=`echo $l | awk -F"\t" '{print $1}'`
	toBeReplaced=`echo $l | awk -F"\t" '{print $2}' | sed 's/[^-A-Za-z0-9_ ]/\\\&/g'`

	if [ "$fixesFile" == "automatic_fixes_properties.txt" ]; then
	    # in sdrf or condensed-sdrv fix factor/characteristic types only
	    #if [ "$fileTypeToBeFixed" == "sdrf.txt" ]; then
		#perl -pi -e "s|\[${toBeReplaced}\]|[${correct}]|g" $exp/$exp.${fileTypeToBeFixed}
	    if [ "$fileTypeToBeFixed" == "condensed-sdrf.tsv" ]; then
		# In condensed-sdrf, the factor/characteristic type is the penultimate column - so tabs on both sides
		perl -pi -e "s|\t${toBeReplaced}\t|\t${correct}\t|g" $exp/$exp.${fileTypeToBeFixed}
	    else
		# idf
		perl -pi -e "s|\t${toBeReplaced}\t|\t${correct}\t|g" $exp/$exp.${fileTypeToBeFixed}
		perl -pi -e "s|\t${toBeReplaced}$|\t${correct}|g" $exp/$exp.${fileTypeToBeFixed}
	    fi
	elif [ "$fixesFile" == "automatic_fixes_values.txt" ]; then
	    #if [ "$fileTypeToBeFixed" == "sdrf.txt" ]; then
		#perl -pi -e "s|\t${toBeReplaced}\t|\t${correct}\t|g" $exp/$exp.${fileTypeToBeFixed}
		#perl -pi -e "s|\t${toBeReplaced}$|\t${correct}|g" $exp/$exp.${fileTypeToBeFixed}
	    if [ "$fileTypeToBeFixed" == "condensed-sdrf.tsv" ]; then
		# In condensed-sdrf, the factor/characteristic value is the last column - so tab on the left and line ending on the right
		perl -pi -e "s|\t${toBeReplaced}$|\t${correct}|g" $exp/$exp.${fileTypeToBeFixed}
	    fi
	fi
    done
}

# Returns prod or test, depending on the Atlas environment in which the script calling it is running
# It is assuming that all atlasinstall_<env>s are under ${ATLAS_PROD}/sw (it will fail otherwise)
atlas_env() {
    #scriptDir=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
    #atlasInstallSubDir=$(echo $scriptDir | awk -F"/" '{print $8}')
    #echo $atlasInstallSubDir | awk -F"_" '{print $2}'
    echo 'prod'
}



# copy_experiment_from_analysis_to_atlas_exps
# and related functions

get_analysis_path_for_experiment_accession(){
	[ "$1" ] && find $ATLAS_PROD/analysis -maxdepth 4 -type d -name "$1" -print -quit
}

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


## get privacy status for any experiments
## -MTAB- experiments loaded by AE/Annotare uis checked by peach API
## -GEOD-/-ERAD-/-ENAD- are loaded as public from now on
peach_api_privacy_status(){
    expAcc=$1
    exp_import=`echo $expAcc | awk -F"-" '{print $2}'`

    if [ $exp_import == "MTAB" ]; then
        response=`curl -s "http://peach.ebi.ac.uk:8480/api/privacy.txt?acc=$expAcc"`
        if [ -z "$response" ]; then
            echo "WARNING: Got empty response from http://peach.ebi.ac.uk:8480/api/privacy.txt?acc=$expAcc" >&2
            exit 0
        fi
        privacyStatus=`echo $response | awk '{print $2}' | awk -F":" '{print $2}'`

    ## if not MTAB, ie. GEOD or ENAD or ERAD are all loaded as public
    else
        privacyStatus=`echo "public"`
    fi
    echo $privacyStatus
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
    		--include '*/' \
    		--exclude '*archive/**' \
    		--exclude '*condensed-sdrf*' \
    		--include '*.tsv' \
    		--include 'qc/**' \
    		--include '*.xml' \
    		--include '*.txt' \
    		--include '*.png' \
    		--include '*.bedGraph'\
    		--include '*.Rdata' \
    		--include '*.pdf' \
    		--include '*.tsv.gz' \
            --exclude 'logs/**' \
    	    --exclude 'lsf.yaml' \
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
