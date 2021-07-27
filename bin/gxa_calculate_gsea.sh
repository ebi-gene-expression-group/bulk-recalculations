#!/usr/bin/env bash
# Calculate gsea for $expAcc and contrast in col nums $pvalColNum and $log2foldchnageColNum respectively, and output the result to ${expPath}/$analyticsFileRoot.<geneSetType>.gsea.tsv

if [ $# -ne 10 ]; then
        echo "Usage: $0 expAcc geneSetFile analyticsFile pvalColNum log2foldchangeColNum expPath contrastId plotTitle organism geneSetType"
        echo "e.g. $0 "
        exit 1
fi

scriptDir=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
# source ${scriptDir}/../bash_util/generic_routines.sh

expAcc=$1
geneSetFile=$2
analyticsFile=$3
pvalColNum=$4
log2foldchangeColNum=$5
expPath=$6
contrastId=$7
plotTitle=$8
organism=$9
geneSetType=${10}

outputFile=${expPath}/${expAcc}.${contrastId}.${geneSetType}.gsea

# Clean up the previous gsea analysis
rm -rf ${expPath}/${expAcc}.${organism}.${contrastId}.ensgene.${geneSetType}.tsv.aux*
rm -rf ${expPath}/${expAcc}.${organism}.${contrastId}.wbpsgene.${geneSetType}.tsv.aux*
rm -rf ${expPath}/${expAcc}.${organism}.${contrastId}.${geneSetType}.gsea*

# Set up environment to be able to run irap_GSE_piano
# done by conda
# source $IRAP_SINGLE_LIB/irap_environment.sh

geneSetFilePath=`ls $geneSetFile`
if [ $? -eq 0 ]; then

   # Prepare mapping files for decorating the gsea file with gene set accessions
   if [ "$geneSetType" == "reactome" ]; then
       mappingFile=$BIOENTITIES_PROPERTIES_PATH/$geneSetType/${organism}.${geneSetType}.tsv.decorate.aux
   elif [ "$geneSetType" == "interpro" ]; then
       mappingFile=$BIOENTITIES_PROPERTIES_PATH/$geneSetType/interproIDToTypeTerm.tsv.decorate.aux
   elif [ "$geneSetType" == "go" ]; then
       mappingFile=$BIOENTITIES_PROPERTIES_PATH/$geneSetType/goIDToTerm.tsv.decorate.aux
   fi
   # Make sure the mapping file has content.
   if [ ! -s "$mappingFile" ]; then
       echo "ERROR: Mapping file for GSEA is empty: $mappingFile"
       exit 1
   fi

   set -v
   irap_GSE_piano --tsv=$analyticsFile --pvalue-col=$pvalColNum --foldchange-col=$log2foldchangeColNum --title="$plotTitle" --pvalue=0.05 --gs_fdr=0.1 --method=fisher-exact --dup-use-best --plot-annot-only --top=10 --minsize 5 --maxsize 100 --descr $mappingFile --go=$geneSetFilePath --out=$outputFile 2>&1
   if [ $? -ne 0 ]; then
	echo "ERROR: Command: 'irap_GSE_piano --tsv=$analyticsFile --pvalue-col=$pvalColNum --foldchange-col=$log2foldchangeColNum --title=\"$plotTitle\" --pvalue=0.05 --gs_fdr=0.1 --method=fisher-exact --dup-use-best --plot-annot-only --top=10 --minsize 5 --maxsize 100 --descr $mappingFile --go=$geneSetFilePath --out=$outputFile' failed" >&2
	exit 1
   fi
   unset -v
   rm -rf ${outputFile}.Rdata

   # Decorate the gsea file with gene set accessions
   if [ -f $outputFile.tsv ]; then
       mv $outputFile.tsv $outputFile.tsv.aux0
       tail -n +2 $outputFile.tsv.aux0 | grep -vP '^\t' | sort -k1,1 -t$'\t'  > $outputFile.tsv.aux1
       head -1 $outputFile.tsv.aux0 | sed 's|^Name|Term\tAccession|' > $outputFile.tsv
       join -t $'\t' -1 1 -2 1 $mappingFile $outputFile.tsv.aux1 >> $outputFile.tsv
       if [ $? -ne 0 ]; then
	   echo "ERROR: Failed to decorate $organism $geneSetType file with gene set accessions"
	   exit 1
       fi
   fi
   # Clean up auxiliary files
   rm -rf $outputFile.tsv.aux*
   rm -rf ${expPath}/${expAcc}.${organism}.${contrastId}.ensgene.${geneSetType}.tsv.aux*
   rm -rf ${expPath}/${expAcc}.${organism}.${contrastId}.wbpsgene.${geneSetType}.tsv.aux*

else
   echo "WARNING: Failed to find $geneSetFile for $expAcc - cannot perform GSEA" >&2
   exit 0
fi
