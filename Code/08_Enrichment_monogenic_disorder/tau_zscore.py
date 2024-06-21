# import numpy as np
import tspex
# import seaborn as sns
import pandas as pd
import sys
import os
curPath = os.path.abspath(os.path.dirname("~/dat/cattle_scdata/Global atlas/All_rds/annotation_rds/CellType/CPM/"))
rootPath = os.path.split(curPath)[0]
sys.path.append(rootPath)


dat = pd.read_csv("~/dat/cattle_scdata/Global atlas/All_rds/annotation_rds/CellType/CPM/combine_all_genes_celltype_cpm_mean.csv", index_col=0)
#dat = dat.drop(dat.columns[[0]], axis = 1)
zscore = tspex.TissueSpecificity(dat, 'zscore', log=True)
tau = tspex.TissueSpecificity(dat, 'tau', log=True)
tissue_specificity1 = zscore.tissue_specificity
tissue_specificity2 = tau.tissue_specificity

output = os.path.join(curPath, "output")
if not os.path.exists(output):
    os.mkdir(output)
outfile1 = os.path.join(output, "zscore_specificity.csv")
outfile2 = os.path.join(output, "tau_specificity.csv")
tissue_specificity1.to_csv(outfile1)
tissue_specificity2.to_csv(outfile2)
print("success!")