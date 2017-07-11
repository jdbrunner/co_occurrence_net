#!/bin/bash

dt=$(date "+%m_%d_%M")

declare -a levels=("genus" "species")

for lv in "${levels[@]}"
do
	mkdir -p ./$dt'_networks/'$lv'/bins/clustered/cyto_input'
	mkdir -p ./$dt'_networks/'$lv'/pears/clustered/cyto_input'
	python3 co_occurrence.py merged_assignment.txt $lv $dt'_networks/'$lv True
	python3 co_occ_net_stats.py $dt'_networks/'$lv'/bins'
	python3 co_occ_net_stats.py $dt'_networks/'$lv'/pears'
	for fl in $dt'_networks/'$lv'/bins/clustered/'*.tsv; do
		python3 cleanup.py $fl
	done
	for fl in $dt'_networks/'$lv'/pears/clustered/'*.tsv; do
		python3 cleanup.py $fl
	done
done

