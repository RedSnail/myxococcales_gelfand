# source /home/demianchenko/anaconda3/bin/activate
# conda activate myxococcales

while IFS="" read -r fam || [ -n "$fam" ]
do
	mkdir ${fam}_pan &> /dev/null
	taxids=$(gawk -v family=${fam} -F'\t' '{if ($7 == family) {print $5}}' myxococcales_table.tsv | sort -u > /dev/stdout)
	echo $taxids
done < families_for_analysis.list
