import glob,re
import argparse

'''
This file transfer vep vcf output format to tab format.
For multi-alt alleles, vep only use the 1st id, so here
we make id mapping to the correct alternative id.
'''

parser = argparse.ArgumentParser(description='change vep vcf output to tab')
parser.add_argument('-i','--vcf',action='store',dest='input',help='input vcf path')
parser.add_argument('-o','--out',action='store',dest='output',help='output txt file')

args = parser.parse_args()
vcf_path = args.input
tab_path = args.output




def vep_vcf2tab(vcf, tab):
    out_columns=['#Uploaded_variation','Location','Allele','Gene',
                 'Feature Feature_type','Consequence','cDNA_position',
                 'CDS_position','Protein_position','Amino_acids','Codons',
                 'Existing_variation','Extra']

    extra_fields=('IMPACT,STRAND,FLAGS,VARIANT_CLASS,HGVSc,HGVSp,LoF,LoF_info,'+ \
                            'CAROL,Condel,gnomad_OE,LoFtool').split(',')
    begin = True
    with open(vcf) as in_vcf, open(tab,'w') as out:
        # iterate each line in vcf format vep annotation
        for line in in_vcf:
            if line.startswith('##INFO=<ID=CSQ'): # extract CSQ column names
                CSQ_keys = re.search('(?<=Format: ).+?(?=\">)', line).group(0).split('|')
            elif line.startswith('#'): # other header
                if line.startswith('##FILTER') or line.startswith('##INFO'):
                    continue
                else:
                    out.write(line)
            else:
                if begin:
                    out.write('\t'.join(out_columns) + '\n')
                    begin = False
                # split each column of vep annotation in vcf format
                chrom,pos,vari_id,ref,alt,qual,filt,info,form=line.split('\t')
                vari_ids = vari_id.split(';')
                alts = alt.split(',')
                # extract CSQ content
                CSQ_items = re.search('(?<=CSQ=).+?(?=$)', info).group(0).split(',')
                # iterate each annotation item in CSQ
                for item in CSQ_items:
                    CSQ_dict = {}
                    CSQ_values = [i if i!='' else '-' for i in item.split('|')]
                    for k,v in zip(CSQ_keys, CSQ_values):
                        CSQ_dict[k] = v
                    allele_num = int(CSQ_dict['ALLELE_NUM']) -1
                    # generate tab output line
                    out_line = [vari_ids[allele_num], pos, CSQ_dict['Allele'],
                               CSQ_dict['Gene'],CSQ_dict['Feature'],CSQ_dict['Feature_type'],
                               CSQ_dict['Consequence'],CSQ_dict['cDNA_position'],
                               CSQ_dict['CDS_position'],CSQ_dict['Protein_position'],
                               CSQ_dict['Amino_acids'],CSQ_dict['Codons'],
                               CSQ_dict['Existing_variation']]
                    extra_field = '' # last column in tab format
                    for f in extra_fields:
                        if CSQ_dict[f] != '-':
                            extra_field = extra_field+f+'='+CSQ_dict[f]+';'
                    out_line.append(extra_field[:-1])
                    out.write('\t'.join(out_line) + '\n')


# path = '/lustre/workspace/projects/IB/P19_014_UKBB_300K_Vep/data'
vcfs = glob.glob(vcf_path + '/*.vcf')
for vcf in vcfs:
    tab = tab_path + '/' + vcf.split('/')[-1]
    if not os.path.exists(tab):
        vep_vcf2tab(vcf,tab)
