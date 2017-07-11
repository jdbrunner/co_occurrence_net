#!/bin/bash

declare -a fld=(07_10_47_networks/genus/bins/*_adj.tsv)

for (( i = 0 ; i < ${#fld[@]} ; i=$i+1 ));
do
	echo $i
	echo ${fld[${i}]}
done
