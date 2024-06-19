##Reference: 
#https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorials

#Expression matrix (Cellranger)
path1 = "~/Bos_taurus_ARS_UCD1.2.gtf"  #change to your genome annotation file 
path2 = "~/Bos_taurus_ARS_UCD1.2.filtered.gtf"  #change to your corresponding name
path3 = "~/index/Bos_taurus"  #change to your specie name
path4 = "~/Bos_taurus_ARS_UCD1.2.fa"  #change to your reference genome
path5 = "~/sample_id"
#(1）filtering annotation information：choose the information you really interested in
cellranger mkgtf $path1 $path2\
                 --attribute=gene_biotype:protein_coding \
                 --attribute=gene_biotype:lncRNA \
                 --attribute=gene_biotype:antiBRse \
                 --attribute=gene_biotype:miRNA \
                 --attribute=gene_biotype:IG_LV_gene \
                 --attribute=gene_biotype:IG_V_gene \
                 --attribute=gene_biotype:IG_V_pseudogene \
                 --attribute=gene_biotype:IG_D_gene \
                 --attribute=gene_biotype:IG_J_gene \
                 --attribute=gene_biotype:IG_J_pseudogene \
                 --attribute=gene_biotype:IG_C_gene \
                 --attribute=gene_biotype:IG_C_pseudogene \
                 --attribute=gene_biotype:TR_V_gene \
                 --attribute=gene_biotype:TR_V_pseudogene \
                 --attribute=gene_biotype:TR_D_gene \
                 --attribute=gene_biotype:TR_J_gene \
                 --attribute=gene_biotype:TR_J_pseudogene \
                 --attribute=gene_biotype:TR_C_gene

#（2） Index of reference genome（0.5h）
cellranger mkref --genome=$path3 \
         --fasta=$path4 \
         --genes=$path2 \
         --memgb=50 \
--nthreads=16 \


#（3）counts matrix（4h）
cellranger count --id=$sample_id \ 
              --transcriptome= $path3 \
              --fastqs= $path5 \
              --sample=$sample_id \
              --localcores=8 \
              --localmem=64

#Important output folder: ~/sample_id/outs/filtered_feature_bc_matrix
#Turn to the seurat analysis

