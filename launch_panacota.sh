#$ -S /bin/bash
#$ -pe smp 4
#$ -cwd
#PBS -l walltime=10:00:00,mem=50gb

send_word execution_started
source ~/anaconda3/bin/activate pan_no_prokka

while IFS="" read -r fam || [ -n "$fam" ]
do
	mkdir ${fam}_pan &> /dev/null
	taxids=$(gawk -v family=${fam} -F'\t' '{if ($7 == family) {print $5}}' myxococcales_table.tsv | sort -u | gawk 'BEGIN{ORS=","} {print }')
	send_word $fam
	PanACoTA prepare -t $taxids -l complete,chromosome -o /home/demianchenko/myxococcales_gelfand/${fam}_pan -p 4 --min_dist 1e-02 --max_dist 0.2
	PanACoTA annotate  --info ${fam}_pan/LSTINFO-*.txt -r ${fam}_pan --threads 4 -n ${fam:0:4}
	send_files ${fam}_pan/*.log
	# rm ${fam}_pan/*.log
done < ~/myxococcales_gelfand/families_for_analysis.list

send_word annotation_done
