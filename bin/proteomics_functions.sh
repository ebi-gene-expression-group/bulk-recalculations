#!/bin/bash

check_tsv_data_file_not_empty_of_data() {
  data_file=$1
  data_col=$2
  if [ ! -f "$data_file" ]; then
    echo "WARNING: no data file $data_file available - this could be fine for some older datasets."
    return 0
  fi
  lines_first_data_col=$(cut -f $data_col $data_file | sort -u | wc -l);
  if [ "$?" -eq 0 ]; then
    if [ "$lines_first_data_col" -le 1 ]; then
      echo "ERROR: $data_file has no data (Column ${data_col} is empty) - this is a major error, seen in the past with older dataset."
      echo "Data will need to be recovered from archives or public/fallback instances."
      echo "If this is recent, you should see correct data available in .snapshot"
      exit 1
    fi
  fi
}

rename_files_baseline(){
  expAcc=$1

  targetAnalysisMethods=${expAcc}-analysis-methods.tsv
  if [ -s "analysis-methods.tsv" ]; then
    mv analysis-methods.tsv $targetAnalysisMethods
  elif [ ! -f "$targetAnalysisMethods" ]; then
    # Only show an error if the target file to produce is not there.
    echo "ERROR: analysis-methods.tsv file doesn't exist for $expAcc to produce $targetAnalysisMethods"
    exit 1
  fi

  data_files=( *MappedToGeneID*.txt )
  if [ -f "${data_files[0]}" ]; then
    echo "Found MappedToGeneID files, copying..."
    cp ${data_files[0]} ${expAcc}.tsv.undecorated.backup
  fi
}

rename_files_differential(){
  expAcc=$1
  if [ -s "E-PROT-analysis-methods.tsv" ]; then
    mv E-PROT-analysis-methods.tsv ${expAcc}-analysis-methods.tsv
  fi

  if [ ! -s "${expAcc}-analysis-methods.tsv" ]; then
    echo "ERROR: exiting as ${expAcc}-analysis-methods.tsv could not be generated or found"
    exit 1
  fi

  if [ -s "${expAcc}-analysis.txt" ]; then
    # this can probably be removed.
    echo "WARNING: found ${expAcc}-analysis.txt, please prefer ExpressionAtlas_FoldChange_groupcomparison.txt"
    cp ${expAcc}-analysis.txt ${expAcc}-analytics.tsv.undecorated
  fi
  if [ -s "ExpressionAtlas_FoldChange_groupcomparison.txt" ]; then
    cp ExpressionAtlas_FoldChange_groupcomparison.txt ${expAcc}-analytics.tsv.undecorated
  fi

  if [ ! -s "${expAcc}-analytics.tsv.undecorated" ]; then
    echo "ERROR: ${expAcc}-analytics.tsv.undecorated is absent or empty for $expAcc. Exiting."
    exit 1
  fi

  dos2unix ${expAcc}-analytics.tsv.undecorated
  ## remove gene name column
  sed -i -r 's/(\s+)?\S+//2' ${expAcc}-analytics.tsv.undecorated
 }

