#!/bin/bash
#SBATCH -p XiaoYueHe
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=170GB
#Date:  2023-05-25
#Need to modify
# /work/home/sdxgroup01/Workspace/20240124rds/Spleen.rds

# Main

#parameter setting
cpu_num="30"
wpath=~/Workspace/Part2/pyscenic;cd ${wpath}

#input file prepare
tf_data="${wpath}/data/allTFs_hg38.txt"
motif_data="${wpath}/data/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
tss_data="${wpath}/data/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
rdspath="~/rds_label_add"

file=$1
tissue=$(basename "$file" .rds)
# opath="${wpath}/Liver/${tissue}$2"
opath="${wpath}/output20240305/${tissue}"
mkdir -p "${opath}" && cd "${opath}" 

# echo ${opath}
#01input file prepare
Rscript ${wpath}/code/pyscenic/01prepare_data.R -i "${file}" -o "${opath}"
python ${wpath}/code/01create_loom_file_from_scanpy.py -i "${opath}/${tissue}_exp.csv"

# # #02main
pyscenic grn "${opath}/${tissue}.loom" "${tf_data}" --num_workers "${cpu_num}" --output "${opath}/01_${tissue}_adjacency.tsv" --method grnboost2

pyscenic ctx "${opath}/01_${tissue}_adjacency.tsv" "${tss_data}" --annotations_fname "${motif_data}" --expression_mtx_fname "${opath}/${tissue}.loom" \
--mode "dask_multiprocessing" --output "${opath}/02_${tissue}_regulons.csv" --num_workers "${cpu_num}" --mask_dropouts

pyscenic aucell "${opath}/${tissue}.loom" "${opath}/02_${tissue}_regulons.csv" --output "${opath}/03_${tissue}_aucell.loom" --num_workers "${cpu_num}"

Rscript ${wpath}/code/02calcRSS_by_scenic.R -l "${opath}/03_${tissue}_aucell.loom" -m "${opath}/${tissue}_metadata.csv" -c cell_type

