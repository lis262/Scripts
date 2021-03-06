'''
Created by Shangzhong.Li@pfizer.com in 2020.
This file plot the boxplot of snp affect on the cell type abundance.
Input parameters:
snp: needs to be in the format chr:pos
'''
import argparse
import os
import pandas as pd
import tabix
import gzip
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
plt.style.use('ggplot')

parser = argparse.ArgumentParser(description='plot genotye count for a snp')


parser.add_argument('-s','--snp',action='store',dest='snp',help='snp id')
parser.add_argument('-c','--cell',action='store',dest='cell',help='cell type')
parser.add_argument('-t','--tumor',action='store',dest='tumor',help='tumor type')
parser.add_argument('-p','--race',action='store',dest='race',help='population',default='eur')
parser.add_argument('-o','--out',action='store',dest='out',help='output file')

args = parser.parse_args()
snp = args.snp
cell = args.cell
tumor = args.tumor
race = args.race
out_fn = args.out

path = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor'
work_path = path + '/plot'
tumor_path = work_path + '/' + tumor

if not os.path.exists(tumor_path):
    os.makedirs(tumor_path)
if race == 'eur':
    sample_fn = path  + '/f01_eur_samples4_qtl.txt'
elif race == 'all':
    sample_fn = path + '/f01_all_samples4_qtl.txt'


def get_samples(sample_df, tumor):
    # get tumor samples in the cacer in TCGA
    tumor_df = sample_df.query('ctype in @tumor')
    samples = tumor_df['IID'].tolist()
    tumor_samples = ['-'.join(i.split('-')[:3]) for i in samples]
    print('number of tumor samples in cacner',tumor+':',len(tumor_samples))
    return tumor_samples

def get_cyto_scores(cyto_df,tumor_samples):
    # extract cytoreason subset
    # cyto_fn = path + '/xcell.txt'
    cyto_df = cyto_df[cyto_df.index.isin(tumor_samples)]
    # remove cell types not in the cancer 
    t_cyto_df = cyto_df.transpose()
    t_cyto_df = t_cyto_df[t_cyto_df.isna().sum(axis=1) < t_cyto_df.shape[1] * 0.1]
    cyto_df = t_cyto_df.transpose()
    # get samples
    cyto_samples = cyto_df.index.tolist()
    print(len(cyto_samples),'samples have cytoreason info')
    # get cells
    cells = cyto_df.columns.tolist()
    return cyto_df,cyto_samples,cells

def get_genotype(snp,gt_path,cyto_samples):
    # 1. find the chromosome file with genotype
    chrom, pos= snp.split(':')[:2]
    chr_fn = gt_path + '/chr' + chrom + '.traw.gz'
    # 2. get the header
    with gzip.open(chr_fn,'rt') as f:
        head = f.readline().strip().split('\t')
        iids = ['-'.join(p[2:].split('-')[:3]) for p in head[6:]]
        indexes = [iids.index(s) for s in cyto_samples]
        tb = tabix.open(chr_fn)
        records = tb.query(chrom,int(pos), int(pos)+1)
        for record in records:
            record = record[6:]
            gt = [record[i] for i in indexes]
            gt_df = pd.DataFrame({'gt':gt},index=cyto_samples)
    return gt_df

def get_gt_with_cyto_score(snp,cell,gt_path,cyto_df,cyto_samples):
    # get genotype
    gt_df = get_genotype(snp,gt_path,cyto_samples)
    # get genotype and merge with cytoreason
    cell_gt_df = pd.merge(cyto_df[cell],gt_df,left_index=True,right_index=True)
    return cell_gt_df



gt_path = '/hpc/grid/wip_drm_targetsciences/projects/TCGA/SNP-Array/IMPUTE/plink'
# 1. get tumor sample
sample_df = pd.read_csv(sample_fn,sep='\t',header=0)
tumor_samples = get_samples(sample_df, tumor)
# 2. extract cytoreason score
cyto_fn = path + '/f02_cytoreason_tcga.tsv'
cyto_df = pd.read_csv(cyto_fn,sep='\t',header=0,index_col=0)
cyto_df,cyto_samples,cells = get_cyto_scores(cyto_df,tumor_samples)
# 3. get gt with cyto score
cell_gt_df = get_gt_with_cyto_score(snp,cell,gt_path,cyto_df,cyto_samples)

# 4. plot for a single plot
ax = sns.boxplot(y=cell,x="gt",data=cell_gt_df,boxprops=dict(alpha=.5))
ax = sns.swarmplot(y=cell,x="gt",data=cell_gt_df)
ax.set_title(snp)
plt.savefig(out_fn)

# ax = cell_gt_df.boxplot(column=cell,by='gt')
# ax.set_title('{snp} cytoreason {pop} patients in {tumor}'.format(snp=snp,pop=race,tumor=tumor))
# print(cell_gt_df['gt'].value_counts())
# plt.suptitle('')
# ax.set_ylabel(cell)
# plt.savefig(out_fn)

#------------ the following is for plotting multiple snps at once for a cancer