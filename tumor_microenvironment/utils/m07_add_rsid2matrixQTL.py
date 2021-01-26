'''
Created by Shangzhong.Li@pfizer.com in 2020.
this file add rsid to the matrixQTL results
'''
import argparse
import pandas as pd
import tabix,glob
import subprocess,os
import time

start = time.time()

parser = argparse.ArgumentParser(description='add rsid')


parser.add_argument('-i','--input',action='store',dest='input',help='input matrixQTL file')
parser.add_argument('-a','--anno',action='store',dest='path',help='annotation path')


args = parser.parse_args()
qtl_fn = args.input
anno_path = args.path

# def get_rsid(snp):
#     chrom,pos = snp.split(':')[:2]
#     tb = tabix_dict['chr'+chrom]
#     records = tb.query(chrom,int(pos),int(pos)+1)
#     for record in records:
#         if record[6] != '-':
#             return record[6]
#         else:
#             return ':'.join(record[:2])

# path = os.path.dirname(qtl_fn)
# if path == '':
#     path = os.getcwd()


# tabix_dict = {}  # {ch:tabix_obj}
# anno_vcfs = glob.glob(anno_path + '/*.tsv.gz')
# for anno in anno_vcfs:
#     chrom = anno.split('/')[-1].split('.')[0]
#     tabix_dict[chrom] = tabix.open(anno)

snp_id_map_fn = anno_path + '/snp_id_map.txt'
snp_id_dict = {}
with open(snp_id_map_fn) as in_f:
    next(in_f)
    for line in in_f:
        snp, rsid = line.strip().split('\t')
        snp_id_dict[snp] = rsid

def get_rsid(snp):
    return snp_id_dict[snp]

header = True
inter_fn = path + '/inter.gz'
for qtl_df in pd.read_csv(qtl_fn,sep='\t',header=0,compression='gzip',chunksize=1000000):
    qtl_df['marker'] = qtl_df['SNP'].map(lambda x: get_rsid(x))
    qtl_df['REF'] = qtl_df['SNP'].map(lambda x: x.split(':')[2])
    qtl_df['ALT'] = qtl_df['SNP'].map(lambda x: x.split(':')[3])
    qtl_df['std'] = qtl_df['beta'].div(qtl_df['t-stat'])
    if header:
        qtl_df.to_csv(inter_fn,sep='\t',index=False,compression='gzip')
        header=False
    else:
        qtl_df.to_csv(inter_fn,sep='\t',index=False,compression='gzip',mode='a',header=False)

    end = time.time()
    print(f"run time is {end - start}")
os.rename(inter_fn, qtl_fn)