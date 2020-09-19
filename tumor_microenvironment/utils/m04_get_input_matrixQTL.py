'''
Created by Shangzhong Li on 2020/08/12
this file prepares the inputfile for matrixQTL

* if you want to use european samples, use line 34, otherwise use line 35.
* If you want to use the cytoreason scores, use line 45, if you want to use xcell scores, use line 46.
* need to change the name of the output, the foramt is cancer_population_scoreSource_scoreFormat.txt 
'''
import pandas as pd
import os
import argparse
import numpy as np
import glob, gzip

parser = argparse.ArgumentParser(description='prepare matrixQTL input files')
parser.add_argument('-c','--cancer',action='store',dest='cancer',help='cancer type')
parser.add_argument('-m','--mode',action='store',dest='mode',help='mode of cell abundance score',default='log')
parser.add_argument('-r','--race',action='store',dest='race',help='population',default='eur')
parser.add_argument('-s','--score',action='store',dest='score',help='population',default='cyto')

args = parser.parse_args()
tumor_prefix = args.cancer
mode = args.mode
race = args.race
score = args.score # cell abundance score

path = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor'
# expression file, in this case is the cell type abundance score
gwas_path = path + '/f03_gwas'
if not os.path.exists(gwas_path):
    os.mkdir(gwas_path)

tumor_path = gwas_path + '/' + tumor_prefix
if not os.path.exists(tumor_path):
    os.mkdir(tumor_path)

if race == 'eur':
    eur_sample_fn = gwas_path  + '/f01_eur_samples4_gwas.txt'
elif race == 'all':
    eur_sample_fn = gwas_path  + '/f01_all_samples4_gwas.txt'
eur_filter_df = pd.read_csv(eur_sample_fn,sep='\t',header=0)

tumors = [tumor_prefix]
eur_tumor_df = eur_filter_df.query('ctype in @tumors')
eur_tumor_df[['IID','IID']].to_csv(
  tumor_path+'/eur_samples.txt',header=None, index=False,sep=' ')
eur_samples = eur_tumor_df['IID'].tolist()
eur_tumor_samples = ['-'.join(i.split('-')[:3]) for i in eur_samples]
print('number of european tumor samples',len(eur_tumor_samples))

# extract cytoreason subset
if score == 'cyto':
    cyto_fn = gwas_path + '/f02_cytoreason_tcga.tsv'
elif score == 'xcell':
    cyto_fn = gwas_path + '/f02_xcell.txt'
cyto_df = pd.read_csv(cyto_fn,sep='\t',header=0,index_col=0)
cyto_df = cyto_df[cyto_df.index.isin(eur_tumor_samples)]
cyto_samples = cyto_df.index.tolist()
print(len(cyto_samples),'samples have cytoreason info')

#-------  covariates file for matrix QTL
covari_fn = gwas_path + '/f02_covariates.tsv'
covari_df = pd.read_csv(covari_fn,sep='\t',header=0,index_col=0)
covari_df = covari_df[covari_df.index.isin(eur_tumor_samples)]
covari_samples = covari_df.index.tolist()
print(len(covari_samples),'samples have covariate info')
analyze_samples = list(set(cyto_samples).intersection(covari_samples))
print(len(analyze_samples),'samples are shared between tcga and cytoreason')
covari_df = covari_df.transpose()
covari_df.index.name = 'covariate'
if covari_df.loc['gender',:].tolist().count(1) in [0,covari_df.shape[1] - covari_df.loc['gender',:].isna().sum()]:
    covari_df = covari_df.drop('gender')
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
        with gzip.open(in_gt,'rt') as in_f, open(out_gt_fn,'a') as out_f:
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

qtl_fn = ('{c}_{pp}_{score}_log.txt').format(c=tumor_prefix,pp=race,score=score)
cmd = ('Rscript {r} {p} {out}').format(r=matrixQTL,p=tumor_path,out=qtl_fn)
os.system(cmd)
os.system('gzip {qtl}'.format(qtl=qtl_fn))


#------------------------ add rsid to the results
vari2rsid = '~/Code/Scripts/tumor_microenvironment/utils/m07_add_rsid2matrixQTL.py'
anno_path = '/hpc/grid/wip_drm_targetsciences/projects/TCGA/SNP-Array/IMPUTE/vep_anno_vcf'

cmd = ('python {code} -i {qtl} -a {anno}').format(code=vari2rsid,qtl=tumor_path + '/' + qtl_fn+'.gz',anno=anno_path)
os.system(cmd)
