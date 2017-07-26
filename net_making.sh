#!/bin/bash

dt=$(date "+%m_%d_%H_%M")

declare -a levels=("genus" "species")

for lv in "${levels[@]}"
do
	mkdir -p ./$dt'_networks/'$lv'/bins'
	mkdir -p ./$dt'_networks/'$lv'/pears'
	python3 co_occurrence.py merged_assignment.txt $lv $dt'_networks/'$lv False True
	python3 cluster_net.py $dt'_networks/'$lv'/bins'
	python3 cluster_net.py $dt'_networks/'$lv'/pears'

	cd $dt'_networks/'$lv'/pears'
	npear=$(find . -name '*cor_adj.tsv')
	npearthr=$(find . -name '*thr_adj.tsv')
# 	cd ../'bins'
# 	nbinned=$()
	cd ../../..
	
# 	echo $npear
# 	echo $npearthr
			
	python3 network_stats.py merged_assignment.txt $lv False True $dt'_networks/'$lv'/pears/'$npear pears
	echo $npear
	python3 network_stats.py merged_assignment.txt $lv False True $dt'_networks/'$lv'/pears/'$npearthr pears_thr
	echo $nearthr
# 	python3 network_stats.py merged_assignment.txt $lv False True $dt'_networks/'$lv'/bins/' binned
done