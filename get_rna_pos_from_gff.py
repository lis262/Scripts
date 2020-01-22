# extract transcript start and end position from gff
import pandas as pd
import re
import argparse

parser = argparse.ArgumentParser(description='get gene/transcipt position from gff file')
parser.add_argument('-i','--gff',action='store',dest='input',help='path to input gff file')
parser.add_argument('-o','--out',action='store',dest='output',help='path to output txt file')

args = parser.parse_args()
in_gff = args.input
out = args.output


gff_df = pd.read_csv(in_gff,sep='\t',header=None, comment='#')

rna_df = gff_df[gff_df[2].values=='mRNA']
rna_df = rna_df.reset_index(drop=True)

rna_df['rna_id']=rna_df[8].map(lambda x:re.search('(?<=transcript:).+?(?=;)',x).group(0))
rna_df['gene_id']=rna_df[8].map(lambda x:re.search('(?<=gene:).+?(?=;)',x).group(0))
rna_df['gene_name']=rna_df[8].map(lambda x:re.search('(?<=Name=).+?(?=;)',x).group(0))
rna_df['gene_name'] = rna_df['gene_name'].map(lambda x: x.split('-')[0])

rna_df.rename(columns={0:'chrom',3:'start',4:'end',6:'strand'},inplace=True)

rna_df[['chrom','start','end','rna_id','gene_id','gene_name','strand']].to_csv(
out,sep='\t',index=False)
