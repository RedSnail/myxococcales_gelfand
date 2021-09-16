#!/bin/bash

gawk -f handle_lstinfo.gawk $1_pan/LSTINFO-LSTINFO-${1:0:4}.lst | tail -n +2 | sort -k1,1 > $1_pan/name_mapping.tsv
gawk -F'\t' '{print $1 "\t" $10}' Bulk_pan/assembly* | tail -n +2 | sort -k1,1  > $1_pan/species_mapping.tsv
#paste $1_pan/name_mapping.tsv $1_pan/species_mapping.tsv > $1_pan/names_mapping.tsv
join -t$'\t' $1_pan/name_mapping.tsv $1_pan/species_mapping.tsv > $1_pan/names_mapping.tsv
# rm $1_pan/name_mapping.tsv
# rm $1_pan/species_mapping.tsv

while read assembly_ac new_name species 
do
	sed -i "s|${new_name}|${species}|g" $1_pan/*.treefile
done < $1_pan/names_mapping.tsv

