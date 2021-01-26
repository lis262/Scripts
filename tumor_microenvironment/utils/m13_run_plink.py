'''
Created by Shangzhong.Li@pfizer.com on 2021/01/15, this is prepare input to run bianry gwas
'''
import argparse
import os,glob,gzip,tabix
import pandas as pd

parser = argparse.ArgumentParser(description='do case-control GWAS')

parser.add_argument('-c','--cancer',action='store',dest='cancer',help='cancer')
parser.add_argument('-r','--race',action='store',dest='race',help='population',default='eur')

# prepare the plink input
args = parser.parse_args()

cancer = args.cancer
race = args.race


clst_fn = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor/p06_cluster/cluster.txt'
clst_df = pd.read_csv(clst_fn,sep='\t',header=0)
dplt_sps = clst_df.query('cancer == "BRCA" and Stratification == "2"')['id'].tolist()
inft_sps = clst_df.query('cancer == "BRCA" and Stratification == "3"')['id'].tolist()

# prepare plink format
plink2 = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor/software/plink2'
raw_path = '/hpc/grid/wip_drm_targetsciences/projects/TCGA/SNP-Array/IMPUTE'
dose_fns = sorted(glob.glob(raw_path+'/*.dose.vcf.gz'))


