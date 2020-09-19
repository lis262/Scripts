import argparse
import pandas as pd
import tabix
import subprocess

parser = argparse.ArgumentParser(description='add rsid')


parser.add_argument('-i','--input',action='store',dest='input',help='input matrixQTL file')
parser.add_argument('-a','--anno',action='store',dest='path',help='annotation path')


args = parser.parse_args()
qtl_fn = args.input
anno_path = args.path

def get_rsid(snp,anno_path):
    chrom,pos = snp.split(':')[:2]
    anno_fn = anno_path + '/chr' + chrom + '.dose.tsv.gz'
    tb = tabix.open(anno_fn)
    records = tb.query(chrom,int(pos),int(pos)+1)
    for record in records:
        if record[6] != '-':
            return record[6]
        else:
            return ':'.join(record[:2])

qtl_df = pd.read_csv(qtl_fn,sep='\t',header=0,compression='gzip')
qtl_df['marker'] = qtl_df['SNP'].map(lambda x:get_rsid(x,anno_path))
qtl_df.to_csv(qtl_fn,sep='\t',index=False,compression='gzip')