import argparse
import os, sys

import loompy as lp;
import numpy as np;
import scanpy as sc;
def main():
    parser= argparse.ArgumentParser(description='make input for pySCENIC')
    parser.add_argument('-i', '--input', type=str, required=True, metavar='input_csv')
    args = parser.parse_args()

    x=sc.read_csv(args.input)
    row_attrs = {"Gene": np.array(x.var_names),}
    col_attrs = {"CellID": np.array(x.obs_names)}
    name = os.path.splitext(os.path.basename(args.input))[0].split('_')[0]
    lp.create(name+'.loom',x.X.transpose(),row_attrs,col_attrs)

if __name__ == '__main__':
    main()