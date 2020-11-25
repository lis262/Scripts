import glob
import tabix,gzip
import sys
import argparse



parser = argparse.ArgumentParser(description='extract UKBB variants')

parser.add_argument('-i','--input',action='store',dest='gene',help='input gene file')
parser.add_argument('-o','--output',action='store',dest='table',help='output table file')

args = parser.parse_args()

gene_file = args.gene
out_gene_fn = args.table

genes = []
with open(gene_file) as in_f:
    for line in in_f:
        genes.append(line.strip())


ukbb_path = '/lustre/workspace/projects/IB/P19_012_UKBB_200K_Vep/data'
tab_path = ukbb_path + '/vep_anno_tab_all'
# gene = 'OR4F5'

gene_pos_fn = tab_path + '/gene_region.txt'
tab_files = sorted(glob.glob(tab_path + '/*.tsv.gz'))

def get_gene_region(gene_pos_fn,gene):
    with open(gene_pos_fn) as in_f:
        for line in in_f:
            item = line.strip().split('\t')
            if item[-1] == gene:
                chrom = item[0]
                start = int(item[1])
                end = int(item[2])
                return [chrom, start, end]

with open(out_gene_fn,'w') as out:
    for g in genes:
        chrom, start, end = get_gene_region(gene_pos_fn,g)
        for tab in tab_files:
            ref_c,ref_s,ref_e = tab.split('.')[1].split('_')[-3:]
            ref_s = int(ref_s)
            ref_e = int(ref_e)
            if chrom == ref_c and (ref_s<=start<=ref_e or ref_s<=end<=ref_e):
                # read file, get header
                with gzip.open(tab, 'rt') as tab_h:
                    head = tab_h.readline()
                    out.write(head)
                tb = tabix.open(tab)
                records = tb.query(chrom, start, end + 1)
                for record in records:
                    if record[7] == g:
                        out.write('\t'.join(record) + '\n')