#!/bin/bash

dt=$(date "+%m_%d")

declare -a levels=("genus" "species")

for lv in "${levels[@]}"
do
	mkdir -p ./$dt'_networks/'$lv'/clustered/cyto_input'
done

