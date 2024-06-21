#!/bin/bash

source ~/00envs/anaconda3/bin/activate /work/home/sdxgroup01/00envs/anaconda3/envs/sc_cellchat
wpath='~/Cellchat'
rdspath='~/rds_label_add'
IFS=$'\n'
files=$(find ~/rds_label_add -name "*.rds" -exec basename -s .rds {} \;) 
for i in ${files}
do
    opath="${wpath}/output/$i" && mkdir -p "${opath}"
    sbatch -N 1 -p YongDingHe --job-name="${i}"_cellchat --output="${opath}/${i}"_cellchat.out --error="${opath}/${i}"_cellchat.err --mem=100GB ${wpath}/code/cellchat.R -i ${rdspath}/"${i}".rds -o "${opath}"
done

