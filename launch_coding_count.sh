#$ -S /bin/bash
#$ -cwd
#PBS -d .
#PBS -l walltime=10:00:00,mem=50gb


rm overall_vs_code.tsv

for file in ./*_pan/gff3/*
do
  gawk -F'\t' -f measure_coding.awk $file >> overall_vs_code.tsv
done


send_word "done"
send_file overall_vs_code.tsv
