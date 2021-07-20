#!/bin/bash

datasets download genome taxon "Myxococcales" --assembly-level complete_genome,chromosome --exclude-gff3 --exclude-protein --exclude-rna --exclude-seq --filename myxococcales_meta.zip
dataformat tsv genome --package myxococcales_meta.zip --fields assminfo-level,assminfo-refseq-assm-accession,assmstats-total-sequence-len,organism-name,tax-id | sort -u >  myxococcales_info.tsv
gawk -F'\t' '{print $5}' myxococcales_info.tsv | tail -n +2 > taxids.list

while read id
do
        family_name="$(efetch -db taxonomy -id $id -format xml | xtract -pattern LineageEx -group Taxon -if Rank -equals family -tab '\n' -sep '\t' -element ScientificName,TaxId)"
        if [ -z "$family_name" ]
        then
                family_name="$(efetch -db taxonomy -id $id -format xml | xtract -pattern LineageEx -group Taxon -if Rank -equals suborder -tab '\n' -sep '\t' -element ScientificName,TaxId)"
        fi
        echo $family_name >> families.csv
done < taxids.list
./summarize_meta.R
