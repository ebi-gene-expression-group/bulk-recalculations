#This is the fastest but is only possible if the condensed sdrf has been already written out
#We only condense the sdrf after successful analysis, and analysis needs organism name in e.g. GSEA, so we can't always use this
#This assumes the organism field there is unique (it fetches the first one that appears otherwise)
get_organism_from_condensed_sdrf(){
    cut -f 4,5,6 $1 | grep characteristic$'\t'organism$'\t' | head -n1 | cut -f 3 | to_ensembl_species_lowercase
}
