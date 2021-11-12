#$ -S /bin/bash
#$ -pe smp 4
#$ -cwd
#PBS -l walltime=10:00:00,mem=50gb

send_word execution_started
source ~/anaconda3/bin/activate pan_no_prokka
source ~/PanACoTA/venv/bin/activate

while IFS="" read -r fam || [ -n "$fam" ]
do
	if [[ "$fam" == \#* ]]
	then
		continue
	fi
	echo $fam
	if [[ "$fam" == "WithOut" ]]
	then
	        taxids=$(tail -n +2 myxococcales_table.tsv | gawk -v family=${fam} -F'\t' '{print $5}' | sort -u | gawk 'BEGIN{ORS=","} {print }')
                acs=$(tail -n +2 myxococcales_table.tsv | gawk -v family=${fam} -F'\t' '{print $2}' | sort -u)
                # PanACoTA prepare -t $taxids -l complete,chromosome -o ~/myxococcales_gelfand/${fam}_pan -p 4 --min_dist 0 --max_dist 1
                # mv ${fam}_pan/LSTINFO*.txt ${fam}_pan/LSTINFO-${fam:0:4}.txt
                # PanACoTA annotate  --info ${fam}_pan/LSTINFO-$fam.txt -r ${fam}_pan --threads 4 -n ${fam:0:4}
	        # PanACoTA pangenome -l ~/myxococcales_gelfand/${fam}_pan/LSTINFO-LSTINFO-${fam:0:4}.lst -n ${fam:0:4} -d ~/myxococcales_gelfand/${fam}_pan/Proteins/ -o ~/myxococcales_gelfand/${fam}_pan/ -m proteinortho --threads 8 -v
                # PanACoTA corepers -p ~/myxococcales_gelfand/${fam}_pan/PanGenome-*.lst -o ~/myxococcales_gelfand/${fam}_pan/
                # PanACoTA align -c ~/myxococcales_gelfand/${fam}_pan/PersGenome_*.lst -l ~/myxococcales_gelfand/${fam}_pan/LSTINFO-LSTINFO-*.lst -d ~/myxococcales_gelfand/${fam}_pan/ -o ~/myxococcales_gelfand/${fam}_pan/ --threads 8 -n ${fam:0:4}
                # PanACoTA tree -a ~/myxococcales_gelfand/${fam}_pan/Phylo-*/*.aln -o ~/myxococcales_gelfand/${fam}_pan/ -b 1000
                bash rename_tree.sh ${fam}
                continue
        fi

	if [[ "$fam" == "Bulk" ]]
	then
		taxids=$(tail -n +2 myxococcales_table.tsv | gawk -v family=${fam} -F'\t' '{if ($7 != "Geobacteraceae") {print $5}}' | sort -u | gawk 'BEGIN{ORS=","} {print }')
		acs=$(tail -n +2 myxococcales_table.tsv | gawk -v family=${fam} -F'\t' '{if ($7 != "Geobacteraceae") {print $2}}' | sort -u)

	elif [[ "$fam" == *_* ]]
	then
		arr_fam=(${fam//_/ })
		fam=${arr_fam[1]}
		fams=${arr_fam[0]}
		taxids=$(gawk -v families=${fams} -F'\t' 'BEGIN { split(families, a, "|", seps) } {for (fam in a) {if ($7 == a[fam]) {print $5}}}' myxococcales_table.tsv | sort -u | gawk 'BEGIN{ORS=","} {print }')
		acs=$(gawk -v families=${fams} -F'\t' 'BEGIN { split(families, a, "|", seps) } {for (fam in a) {if ($7 == a[fam]) {print $2}}}' myxococcales_table.tsv | sort -u)

	else
		taxids=$(gawk -v family=${fam} -F'\t' '{if ($7 == family) {print $5}}' myxococcales_table.tsv | sort -u | gawk 'BEGIN{ORS=","} {print }')
		acs=$(gawk -v family=${fam} -F'\t' '{if ($7 == family) {print $2}}' myxococcales_table.tsv | sort -u)
	fi
	mkdir ${fam}_pan &> /dev/null
	mkdir ${fam}_pan/Proteins
	mkdir ${fam}_pan/Genes

	for ac in ${acs[@]}
	do
		cat ~/myxococcales_gelfand/WithOut_pan/LSTINFO-LSTINFO-With.lst | grep $ac
	done | \
	       gawk -F'\t' -v name=${fam:0:4} 'BEGIN {i = 1} {print name ".0921." sprintf("%05d", i) "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6; i += 1}' | \
	       (echo "$(head -n 1 ~/myxococcales_gelfand/WithOut_pan/LSTINFO-LSTINFO-With.lst)" && cat) > ~/myxococcales_gelfand/${fam}_pan/LSTINFO-LSTINFO-${fam:0:4}.lst


	for ac in ${acs[@]}
        do
                cat ~/myxococcales_gelfand/WithOut_pan/LSTINFO-LSTINFO-With.lst | grep $ac
        done | gawk -F'\t' -v name=${fam:0:4} 'BEGIN {i = 1} {print $1 "\t" name ".0921." sprintf("%05d", i); i += 1}' > ${fam}_pan/sname_mapping.tsv

	while read old new
	do
		cp WithOut_pan/Proteins/$old.prt ${fam}_pan/Proteins/$new.prt
		sed -i "s|${old}|${new}|g" ${fam}_pan/Proteins/$new.prt
		cp WithOut_pan/Genes/$old.gen ${fam}_pan/Genes/$new.gen
		sed -i "s|${old}|${new}|g" ${fam}_pan/Genes/$new.gen
	done < ${fam}_pan/sname_mapping.tsv
	PanACoTA pangenome -l ~/myxococcales_gelfand/${fam}_pan/LSTINFO-LSTINFO-${fam:0:4}.lst -n ${fam:0:4} -d ~/myxococcales_gelfand/${fam}_pan/Proteins/ -o ~/myxococcales_gelfand/${fam}_pan/ -m proteinortho --threads 8 -v
	PanACoTA corepers -p ~/myxococcales_gelfand/${fam}_pan/PanGenome-*.lst -o ~/myxococcales_gelfand/${fam}_pan/
	PanACoTA align -c ~/myxococcales_gelfand/${fam}_pan/PersGenome_*.lst -l ~/myxococcales_gelfand/${fam}_pan/LSTINFO-LSTINFO-*.lst -d ~/myxococcales_gelfand/${fam}_pan/ -o ~/myxococcales_gelfand/${fam}_pan/ --threads 8 -n ${fam:0:4}
	PanACoTA tree -a ~/myxococcales_gelfand/${fam}_pan/Phylo-*/*.aln -o ~/myxococcales_gelfand/${fam}_pan/ -b 1000

	# send_files ${fam}_pan/*.log
	# rm ${fam}_pan/*.log

	bash rename_tree.sh ${fam}
done < ~/myxococcales_gelfand/families_for_analysis.list

send_word annotation_done
