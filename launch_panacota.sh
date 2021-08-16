#$ -S /bin/bash
#PBS -d .
#PBS -l walltime=10:00:00,mem=50gb

send_word execution_started
source ~/anaconda3/bin/activate
conda activate myxococcales

while IFS="" read -r fam || [ -n "$fam" ]
do
	mkdir ${fam}_pan &> /dev/null
	taxids=$(gawk -v family=${fam} -F'\t' '{if ($7 == family) {print $5}}' myxococcales_table.tsv | sort -u | gawk 'BEGIN{ORS=","} {print }')
	send_word $taxids
	PanACoTA prepare -t $taxids -l complete,chromosome -o /home/demianchenko/myxococcales_gelfand/${fam}_pan -p 4 --min_dist 1e-02 --max_dist 0.2
done < ~/myxococcales_gelfand/families_for_analysis.list

send_word downloading_done
