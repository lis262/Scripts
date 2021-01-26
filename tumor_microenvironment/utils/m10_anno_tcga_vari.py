'''
Created by Shangzhong.Li@pfizer.com on 2020/10/02.
This annotate the TCGA imputated variants and add genes near 50k bps of the SNPs
'''

import pandas as pd
import tabix,glob,subprocess
from intervaltree import Interval, IntervalTree
import pandas as pd
import re,os
from natsort import natsorted



#--------------------- get interval of genes
def get_gene_interval(gff_fn):
   
    gff_df = pd.read_csv(gff_fn,sep='\t',comment='#',header=None)
    gene_df = gff_df[gff_df[2].values=='gene']
    gene_df = gene_df.reset_index(drop=True)
    cds_df = gff_df[gff_df[2].values == 'CDS']
    cds_df = cds_df.reset_index(drop=True)

    gene_df['gene'] = gene_df[8].map(lambda x:re.search('(?<=gene_name=).+?(?=;)',x).group(0))
    cds_df['gene'] = cds_df[8].map(lambda x:re.search('(?<=gene_name=).+?(?=;)',x).group(0))
    genes = cds_df['gene'].unique().tolist()
    gene_df = gene_df.query('gene in @genes')

    interval_dict = {}
    for idx, row in gene_df.iterrows():
        chrom, start, end = row[0], row[3], row[4]
        chrom = chrom[3:]
        if chrom in interval_dict:
            interval_dict[chrom][start-1:end] = row['gene']
        else:
            interval_dict[chrom] = IntervalTree()
    return interval_dict


def add_near50k_gene(row,interval_dict):
    chrom = row['chr']
    start = int(row['pos']) - 50000 -1
    end = int(row['pos']) + 50000
    itvs = interval_dict[chrom][start:end]
    near_genes = []
    for itv in itvs:
        near_genes.append(itv[2])
    if near_genes == []:
        near_genes = ['-']
    return ';'.join(near_genes)

# annotation path
# this add column [chr and position]
vep_path = '/hpc/grid/wip_drm_targetsciences/projects/TCGA/SNP-Array/IMPUTE/union_vcf/vep_anno_vcf'
anno_path = '/hpc/grid/wip_drm_targetsciences/projects/TCGA/SNP-Array/IMPUTE/vep_anno_vcf'
gff_fn = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/publicDB/gencode.v35lift37.annotation.gff3.gz'

interval_dict = get_gene_interval(gff_fn)
anno_vcfs = natsorted(glob.glob(vep_path + '/*.tsv.gz'))

for vcf in anno_vcfs:
    columns = ['Uploaded_variation', 'Consequence', 'NEAREST',
          'SYMBOL', 'Existing_variation', 'Gene', 'Feature','CANONICAL', 
          'CADD_PHRED']
    out_fn = anno_path + '/' + vcf.split('/')[-1][:-3]
    df = pd.read_csv(vcf, sep='\t',header=None,compression='gzip',comment='#',
                    names=columns)
    df['chr'] = df['Uploaded_variation'].map(lambda x:x.split(':')[0])
    df['pos'] = df['Uploaded_variation'].map(lambda x:x.split(':')[1])
    df['Near50kGene'] = df.apply(lambda row:add_near50k_gene(row,interval_dict),axis=1)
    df[['chr','pos']+columns+['Near50kGene']].to_csv(out_fn,sep='\t',index=False)
    os.system('bgzip -f {f}'.format(f=out_fn))
    subprocess.call('tabix -p bed -f -b 2 -e 2 -S 1 {f}.gz'.format(f=out_fn),shell=True)