#!/bin/bash

dt=$(date "+%m_%d_%H_%M")

declare -a levels=("genus" "species")

for lv in "${levels[@]}"
do
	mkdir -p ./$dt'_networks/'$lv'/bins'
	mkdir -p ./$dt'_networks/'$lv'/pears/validation_plots'
	python3 co_occurrence.py merged2.txt $lv $dt'_networks/'$lv True 300
	python3 cluster_net.py $dt'_networks/'$lv'/bins'
	python3 cluster_net.py $dt'_networks/'$lv'/pears'


	
	for np in $dt'_networks/'$lv'/pears/'*cor_adj.tsv; do
		python3 network_stats.py merged2.txt $lv True $np pears
	done
	for npt in $dt'_networks/'$lv'/pears/'*thr_adj.tsv; do
		python3 network_stats.py merged2.txt $lv True $npt pears_thr
	done
	
	python3 falsepm.py merged2.txt $dt'_networks/'$lv'/pears'
	
# 	python3 network_stats.py merged_assignment.txt $lv False True $dt'_networks/'$lv'/bins/' binned
done