'''
Created by Shangzhong Li on 2020/08/12
this file prepares the inputfile for matrixQTL

* If you want to use the cytoreason scores, use line 45, if you want to use xcell scores, use line 46.
* need to change the name of the output, the foramt is cancer_population_scoreSource_scoreFormat.txt 
'''
import pandas as pd
import os
import argparse
import numpy as np
import glob

parser = argparse.ArgumentParser(description='prepare matrixQTL input files')
parser.add_argument('-c','--cancer',action='store',dest='cancer',help='cancer type')
parser.add_argument('-m','--mode',action='store',dest='mode',help='mode of cell abundance score',default='log')

args = parser.parse_args()
tumor_prefix = args.cancer
mode = args.mode


path = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor'
# expression file, in this case is the cell type abundance score
gwas_path = path + '/f03_gwas'
if not os.path.exists(gwas_path):
    os.mkdir(gwas_path)

tumor_path = gwas_path + '/' + tumor_prefix
if not os.path.exists(tumor_path):
    os.mkdir(tumor_path)

eur_sample_fn = gwas_path  + '/eur_samples4_gwas.txt'
# eur_sample_fn = gwas_path  + '/all_samples4_gwas.txt'
eur_filter_df = pd.read_csv(eur_sample_fn,sep='\t',header=0)

tumors = [tumor_prefix]
eur_tumor_df = eur_filter_df.query('ctype in @tumors')
eur_tumor_df[['IID','IID']].to_csv(
  tumor_path+'/eur_samples.txt',header=None, index=False,sep=' ')
eur_samples = eur_tumor_df['IID'].tolist()
eur_tumor_samples = ['-'.join(i.split('-')[:3]) for i in eur_samples]
print('number of european tumor samples',len(eur_tumor_samples))

# extract cytoreason subset
# cyto_fn = gwas_path + '/cytoreason_tcga.tsv'
cyto_fn = gwas_path + '/xcell.txt'
cyto_df = pd.read_csv(cyto_fn,sep='\t',header=0,index_col=0)
cyto_df = cyto_df[cyto_df.index.isin(eur_tumor_samples)]
cyto_samples = cyto_df.index.tolist()
print(len(cyto_samples),'samples have cytoreason info')

#-------  covariates file for matrix QTL
covari_fn = gwas_path + '/covariates.tsv'
covari_df = pd.read_csv(covari_fn,sep='\t',header=0,index_col=0)
covari_df = covari_df[covari_df.index.isin(eur_tumor_samples)]
covari_samples = covari_df.index.tolist()
print(len(covari_samples),'samples have covariate info')
analyze_samples = list(set(cyto_samples).intersection(covari_samples))
print(len(analyze_samples),'samples are shared between tcga and cytoreason')
covari_df = covari_df.transpose()
covari_df.index.name = 'covariate'
covari_df[analyze_samples].to_csv(tumor_path + '/covariates.txt',sep='\t')
t_cyto_df = cyto_df.transpose()
t_cyto_df.index.name = 'cell'
t_cyto_df = t_cyto_df[t_cyto_df.isna().sum(axis=1) < t_cyto_df.shape[1] * 0.1]
if mode == 'log':
    t_cyto_df = t_cyto_df.apply(np.log10)
elif mode == 'normal':
    import scipy.stats as st
    rank_cyto_df = (t_cyto_df.rank()-0.5) / t_cyto_df.rank().max()
    rank_cyto_df = rank_cyto_df.apply(st.norm.ppf)
    rank_cyto_df[analyze_samples].to_csv(tumor_path + '/cell_type_contri.txt',sep='\t')
    t_cyto_df = rank_cyto_df
t_cyto_df[analyze_samples].to_csv(tumor_path + '/cell_type_contri.txt',sep='\t')

# ---------------------- genotype file --------------------------
def gt_file4matrixQTL(gt_path,out_gt_fn):
    '''
    gt_path: folder that have all traw files
    out_gt_fn: out genotype file for matrixQTL
    '''
    in_gt_fns = sorted(glob.glob(in_gt_path + '/*.traw'))
    header = True
    if os.path.exists(out_gt_fn):
        os.remove(out_gt_fn)
    for in_gt in in_gt_fns:
        # prepare SNPfiles
        with open(in_gt) as in_f, open(out_gt_fn,'a') as out_f:
            head = in_f.readline().strip().split('\t')
            head = ['-'.join(c.split('-')[:3])[2:] if '-' in c else c for c in head]
            if header:
                out_f.write('\t'.join(['SNP'] + analyze_samples)+'\n')
                header = False
            indexes = [head.index(s) for s in analyze_samples]
            for line in in_f:
                item = line.strip().split('\t')
                gts = [float(item[i]) for i in indexes]
                out_f.write('\t'.join([item[1]]+[item[i] for i in indexes])+'\n')

in_gt_path = '/hpc/grid/wip_drm_targetsciences/projects/TCGA/SNP-Array/IMPUTE/plink_dosage/'
gt_fn = tumor_path + '/SNP.txt'
gt_file4matrixQTL(in_gt_path, gt_fn)

# ---------------------- run matrixQTL --------------------------
matrixQTL = '~/Code/Scripts/tumor_microenvironment/utils/m05_matrixQTL.R'

cmd = ('Rscript {r} {p} {c}_eur_xcell_log.txt').format(r=matrixQTL,p=tumor_path,c=tumor_prefix)
os.system(cmd)