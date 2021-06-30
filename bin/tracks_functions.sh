generate_differential_tracks() {
  expAcc=$1
  id=$2
  analyticsFile=$3
  gffFile=$4
  label=$5
  expPath=$6

  # tsv2bedGraph provided by irap-components, IRAP envs set there probably
  tsv2bedGraph --dup-pick-highest --out ${expPath}/${expAcc}.${id}.genes.pval --tsv ${analyticsFile} --gtf $gffFile --value_col ${id}.p-value --dup_value_col ${id}.log2foldchange --label "$label" --description "For more information please visit http://www.ebi.ac.uk/gxa/experiments/${expAcc}"
  if [ $? -ne 0 ]; then
    echo "ERROR: Failed to generate ${id}.genes.pval.bedGraph for $expAcc" >&2
    return 1
  fi
  tsv2bedGraph --dup-pick-highest --out ${expPath}/${expAcc}.${id}.genes.log2foldchange --tsv ${analyticsFile} --gtf $gffFile --value_col ${id}.log2foldchange --label "$label" --description "For more information please visit http://www.ebi.ac.uk/gxa/experiments/${expAcc}"
  if [ $? -ne 0 ]; then
    echo "ERROR: Failed to generate ${id}.genes.log2foldchange.bedGraph for $expAcc"  >&2
    return 1
  fi
}

# Generate fpkm tracks for baseline experiment $expAcc
generate_baseline_tracks() {
  expAcc=$1
  # assay id
  id=$2
  expressionsFile=$3
  gffFile=$4
  expPath=$5
  # assay label
  label=$6

  tpmORfpkm=$(basename $expressionsFile | awk -F"-" '{print $4}' | sed 's/.tsv//g')
  gxa_tsv2bedGraph.sh --out ${expPath}/${expAcc}.${id}.genes.expressions_${tpmORfpkm} --tsv $expressionsFile --gtf $gffFile --value_col $id --label "$label" --description "For more information please visit http://www.ebi.ac.uk/gxa/experiments/${expAcc}"
  if [ $? -ne 0 ]; then
    echo "ERROR: Failed to generate ${id}.genes.expressions.bedGraph_${tpmORfpkm} for $expAcc" >&2
    return 1
  fi
}
