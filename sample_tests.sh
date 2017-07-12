#!/bin/bash

dt=$(date "+%m_%d_%M")
mkdir -p test_samples/$dt'_samps/genus'
mkdir -p test_samples/$dt'_samps/species'


for i in `seq 1 2`;do
	python3 rand_samp.py merged_assignment.txt [genus,species] test_samples/$dt'_samps'
done    

declare -a levels=("genus" "species")
nets_flder=$1




for lv in "${levels[@]}";do
	declare -a adjs_p=($nets_flder'/'$lv'/pears/'*_adj.tsv)
	declare -a lists_p=($nets_flder'/'$lv'/pears/'*_list.tsv)
	declare -a adjs_b=($nets_flder'/'$lv'/bins/'*_adj.tsv)
	declare -a lists_b=($nets_flder'/'$lv'/bins/'*_list.tsv)
	for fl in test_samples/$dt'_samps/'$lv'/'*.tsv;do
		for (( i = 0 ; i < ${#adjs_p[@]} ; i=$i+1 ));
		do
			echo ${lists_p[${i}]}
			echo ${adjs_p[${i}]}
			python3 sample_analysis.py $fl ${lists_p[${i}]} ${adjs_p[${i}]} $lv
		done
		for (( j = 0 ; j < ${#adjs_b[@]} ; j=$j+1 ));
		do
			echo ${lists_b[${j}]}
			echo ${adjs_b[${j}]}
			python3 sample_analysis.py $fl ${lists_b[${j}]} ${adjs_b[${j}]} $lv
		done
	done
done

