 #!/bin/bash
IFS=$'\n'
rds_path="/work/home/sdxgroup01/Workspace/20240124rds"
file_list=$(find ${rds_path} -type f -name "*.rds")

Workspace="/work/home/sdxgroup01/Workspace/Part2/pyscenic/output20240305"
for file in ${file_list[@]}
do
    tissue=$(basename "$file" .rds)
    echo "#############"
    echo "$tissue"
    echo "$file"
    sbatch --job-name=${tissue}_pyscenic -p XiaoYueHe --nodes=1 --output="${Workspace}/${tissue}/${tissue}_pyscenic1.out" --error="${Workspace}/${tissue}/${tissue}_pyscenic1.err" /work/home/sdxgroup01/01user/zq/01scrna/code/pyscenic/submit_single_pyscenic.sh "$file" 
done

# find . -type d -exec sh -c 'if [ ! -e "{}"/03*loom ]; then basename "{}"; fi' \;
# find . -type d -exec sh -c 'if [ ! -e "{}"/*exp.csv ]; then basename "{}"; fi' \;
# tissue_values=("Liver")
# # tissue_values=("Liver" "Cerebral cortex" "Bone Marrow" "Tongue" "Circumvallate papilla" "Mammary gland" "Cerebellum" "Heart" "Cervical lymph node" "Esophage")

# for tissue in "${tissue_values[@]}"; do
#     file="${rds_path}/${tissue}.rds"
#     echo "############"
#     echo "$tissue"
#     echo "$file"
#     sbatch --job-name=${tissue}_pyscenic --mem=150GB -p XiaoYueHe --output="${Workspace}/${tissue}/${tissue}_pyscenic.out" --error="${Workspace}/${tissue}/${tissue}_pyscenic.err" /work/home/sdxgroup01/01user/zq/01scrna/code/pyscenic/submit_single_pyscenic.sh "$file"
#     # sbatch --job-name=${tissue}_pyscenic --nodes=1 --cpus-per-task=1 --mem=150GB -p XiaoYueHe --output="${Workspace}/${tissue}/${tissue}_pyscenic.out" --error="${Workspace}/${tissue}/${tissue}_pyscenic.err" /work/home/sdxgroup01/01user/zq/01scrna/code/pyscenic/submit_single_pyscenic.sh "$file"
#     echo "############"
# done
