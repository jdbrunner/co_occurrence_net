#!/bin/bash

dt=$(date "+%m_%d_%M")
mkdir -p test_samples/$dt'_samps/genus/ranked'
mkdir -p test_samples/$dt'_samps/species/ranked'


for i in `seq 1 5`;do
	python3 rand_samp.py merged_assignment.txt [genus,species] test_samples/$dt'_samps'
done    

declare -a levels=("genus" "species")
nets_flder=$1

for lv in "${levels[@]}";do
	for fl in test_samples/$dt'_samps/'$lv'/'*.txt;do
		python3 sample_analysis.py $fl $nets_flder'/'$lv'/pears/clustered/'$lv'__list_clustered.tsv' $nets_flder'/'$lv'/pears/'$lv'__adj.tsv' $lv test_samples/$dt'_samps/'$lv'/ranked'
		python3 sample_analysis.py $fl $nets_flder'/'$lv'/bins/clustered/'$lv'__list_clustered.tsv' $nets_flder'/'$lv'/bins/'$lv'__adj.tsv' $lv test_samples/$dt'_samps/'$lv'/ranked'
	done
done

