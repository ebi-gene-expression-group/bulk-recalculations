get_contrast_colnum() {
    file=$1
    contrastId=$2
    pvalOrFoldChange=$3
    i=0
    # avoid spaces in Gene ID and Gene Name columns
    for col in $(head -1 $file | sed 's/ /_/g' | tr "\t" "\n"); do
      i=$[$i+1]
      if [ "$col" == "$contrastId.$pvalOrFoldChange" ]; then
        echo $i
      fi
    done
}

find_properties_file_gsea() {
    organism=$1
    property=$2

    if [ $property = "reactome" ]; then
      reactomeFile="$BIOENTITIES_PROPERTIES_PATH/reactome/${organism}.reactome.tsv.gsea.aux"
      if [ -s "$reactomeFile" ]; then
        echo "$reactomeFile"
      else
        >&2 echo "No annotation file found for organism $organism and property $property"
      fi
    else
      ensFile="$BIOENTITIES_PROPERTIES_PATH/ensembl/${organism}.ensgene.${property}.tsv"
      if [ -s "$ensFile" ]; then
        echo $ensFile
      else
        wbpsFile="$BIOENTITIES_PROPERTIES_PATH/wbps/${organism}.wbpsgene.${property}.tsv"
        if [ -s "$wbpsFile" ]; then
          echo $wbpsFile
        else
          >&2 echo "No annotation file found for organism $organism and property $property"
          exit 1
        fi
      fi
    fi
}
