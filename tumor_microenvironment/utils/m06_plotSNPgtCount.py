import argparse
import pandas as pd
import tabix
import gzip
import matplotlib
import matplotlib.pyplot as plt
plt.style.use('ggplot')

parser = argparse.ArgumentParser(description='plot genotye count for a snp')


parser.add_argument('-s','--snp',action='store',dest='snp',help='snp id')
parser.add_argument('-c','--cell',action='store',dest='cell',help='cell type')
parser.add_argument('-t','--tumor',action='store',dest='tumor',help='tumor type')
parser.add_argument('-p','--race',action='store',dest='race',help='population',default='eur')
parser.add_argument('-o','--out',action='store',dest='out',help='output file',)

args = parser.parse_args()
snp = args.snp
cell = args.cell
tumor = args.tumor
race = args.race
out_fn = args.out

path = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor'
gwas_path = path + '/f03_gwas'
tumor_path = gwas_path + '/' + tumor

if race == 'eur':
    sample_fn = gwas_path  + '/f01_eur_samples4_gwas.txt'
elif race == 'all':
    sample_fn = gwas_path + '/f01_all_samples4_gwas.txt'

# get tumor samples in the cacer in TCGA
sample_df = pd.read_csv(sample_fn,sep='\t',header=0)
tumor_df = sample_df.query('ctype in @tumor')
tumor_df[['IID','IID']].to_csv(
  tumor_path+'/samples.txt',header=None, index=False,sep=' ')
samples = tumor_df['IID'].tolist()
tumor_samples = ['-'.join(i.split('-')[:3]) for i in samples]
print('number of tumor samples',len(tumor_samples))

# extract cytoreason subset
cyto_fn = gwas_path + '/f02_cytoreason_tcga.tsv'
# cyto_fn = gwas_path + '/xcell.txt'
cyto_df = pd.read_csv(cyto_fn,sep='\t',header=0,index_col=0)
cyto_df = cyto_df[cyto_df.index.isin(tumor_samples)]
cyto_samples = cyto_df.index.tolist()
print(len(cyto_samples),'samples have cytoreason info')

# get genotype and merge with cytoreason
# 1. find the chromosome file with genotype
chrom, pos= snp.split(':')[:2]
gt_path = '/hpc/grid/wip_drm_targetsciences/projects/TCGA/SNP-Array/IMPUTE/plink'
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
cell_gt_df = pd.merge(cyto_df[cell],gt_df,left_index=True,right_index=True)

ax = cell_gt_df.boxplot(column=cell,by='gt')
ax.set_title('{snp} cytoreason {pop} patients in {tumor}'.format(snp=snp,pop=race,tumor=tumor))
print(cell_gt_df['gt'].value_counts())
plt.suptitle('')
ax.set_ylabel(cell)
plt.savefig(out_fn)