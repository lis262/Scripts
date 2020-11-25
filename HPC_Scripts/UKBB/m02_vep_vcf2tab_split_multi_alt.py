'''
Created by Shangzhong.Li@pfizer.com on 2020/05
This file transfer UKBB vep vcf output format to tab format.
For multi-alt alleles, vep only use the 1st id, so here
we make id mapping to the correct alternative id.
'''

import glob,re,os,gzip
import argparse
from collections import OrderedDict
from natsort import natsorted

parser = argparse.ArgumentParser(description='change vep vcf output to tab')
parser.add_argument('-i','--vcf',action='store',dest='input',help='input vcf path')
parser.add_argument('-o','--out',action='store',dest='output',help='output tab path')
parser.add_argument('-s','--start',action='store',dest='start',type=int,
                            help='start index of file',default=0)
parser.add_argument('-e','--end',action='store',dest='end',type=int,
                            help='end index of file',default=2000)

args = parser.parse_args()
vcf_path = args.input
tab_path = args.output
start = args.start
end = args.end

if not os.path.exists(tab_path):
    os.mkdir(tab_path)


def vep_vcf2tab(vcf, tab):
    # this one output all the transcripts instead of the one with worst prediction
    out_columns=['#Uploaded_variation','Location','REF','Allele','Gene',
                 'Feature', 'Feature_type','Consequence','cDNA_position',
                 'CDS_position','Protein_position','Amino_acids','Codons',
                 'Existing_variation','LoF','CADD_PHRED','Extra']

    extra_fields=('IMPACT,STRAND,FLAGS,VARIANT_CLASS,HGVSc,HGVSp,LoF,LoF_info,'+ \
                            'CAROL,Condel,gnomad_oe,LoFtool').split(',')
    begin = True
    with gzip.open(vcf,'rt') as in_vcf, open(tab,'w') as out:
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
                    out_line = [vari_ids[allele_num], pos, ref, CSQ_dict['Allele'],
                               CSQ_dict['Gene'],CSQ_dict['Feature'],CSQ_dict['Feature_type'],
                               CSQ_dict['Consequence'],CSQ_dict['cDNA_position'],
                               CSQ_dict['CDS_position'],CSQ_dict['Protein_position'],
                               CSQ_dict['Amino_acids'],CSQ_dict['Codons'],
                               CSQ_dict['Existing_variation'],
                               CSQ_dict['LoF'],CSQ_dict['CADD_PHRED']]
                    extra_field = '' # last column in tab format
                    for f in extra_fields:
                        if CSQ_dict[f] != '-':
                            extra_field = extra_field+f+'='+CSQ_dict[f]+';'
                    out_line.append(extra_field[:-1])
                    out.write('\t'.join(out_line) + '\n')


#============ The following part pick the most severe transcript =====
# define scores of variants effects

effect_dict = {
    'transcript_ablation':36,
    'splice_acceptor_variant':35,
    'splice_donor_variant':34,
    'stop_gained':33,
    'frameshift_variant':32,
    'stop_lost':31,
    'start_lost':30,
    'transcript_amplification':29,
    'inframe_insertion':28,
    'inframe_deletion':27,
    'missense_variant':26,
    'protein_altering_variant':25,
    'splice_region_variant':24,
    'incomplete_terminal_codon_variant':23,
    'start_retained_variant':22,
    'stop_retained_variant':21,
    'synonymous_variant':20,
    'coding_sequence_variant':19,
    'mature_miRNA_variant':18,
    '5_prime_UTR_variant':17,
    '3_prime_UTR_variant':16,
    'non_coding_transcript_exon_variant':15,
    'intron_variant':14,
    'NMD_transcript_variant':13,
    'non_coding_transcript_variant':12,
    'upstream_gene_variant':11,
    'downstream_gene_variant':10,
    'TFBS_ablation':9,
    'TFBS_amplification':8,
    'TF_binding_site_variant':7,
    'regulatory_region_ablation':6,
    'regulatory_region_amplification':5,
    'feature_elongation':4,
    'regulatory_region_variant':3,
    'feature_truncation':2,
    'intergenic_variant':1,
     '':0
}


def pick_most_severe_CSQ(csq_notes, csq_fields, effect_dict):
    """
    this function picks the most severe annotation for a variant
    {csq_field: value}
    * csq_notes: all csq information
    * csq_fields: column names in the csq fields
    * effect_dict: {effect: score}
    """
    gene_idx = csq_fields.index('Gene')
    
    csq_records = csq_notes.split(',')
    csq_dict = OrderedDict()
    
    effect_idx = csq_fields.index('Consequence')
    for csq in csq_records: # loop each transcript
        csq_values = csq.split('|')
        # for a transcript with &, choose the most severe term
        effect = csq_values[effect_idx]
        gene = csq_values[gene_idx]
        if gene not in csq_dict:
            csq_dict[gene] = {k:'' for k in csq_fields}
            
        if '&' in effect:
            # get severe scores for each SO term in the effect
            rna_effects = effect.split('&')
            scores = [effect_dict[i] for i in rna_effects]
            max_score = max(scores)
            max_idx = scores.index(max_score)
            max_effect = rna_effects[max_idx]
            csq_values[effect_idx] = max_effect # update the effect content
        else:
            max_effect = effect
        # update the server csq item
        if effect_dict[csq_values[effect_idx]] > effect_dict[csq_dict[gene]['Consequence']]:
            for k,v in zip(csq_fields, csq_values):
                csq_dict[gene][k] = v
    return csq_dict


def vep_vcf2tab_worst(vcf, tab):
    out_columns=['#Uploaded_variation','Location','REF','Allele','Gene',
                     'Feature', 'Feature_type','Consequence','cDNA_position',
                     'CDS_position','Protein_position','Amino_acids','Codons',
                     'Existing_variation','LoF','CADD_PHRED','Extra']

    extra_fields=('IMPACT,STRAND,FLAGS,VARIANT_CLASS,HGVSc,HGVSp,LoF,LoF_info,'+ \
                            'CAROL,Condel,gnomad_oe,LoFtool').split(',')

    import gzip,re
    begin = True
    with gzip.open(vcf,'rt') as in_vcf, open(tab,'w') as out:
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
                chrom,pos,vari_id,ref,alt,qual,filt,info,format=line.split('\t')[:9]
                vari_ids = vari_id.split(';')
                alts = alt.split(',')
                # extract CSQ content
                try:
                    csq_notes = re.search('(?<=CSQ=).+?(?=$)', info).group(0)
                except:
                    continue
                # split the multi allele csq notes
                csq_multi = [''] * len(vari_ids)
                for item in csq_notes.split(','):
                    CSQ_dict = {}
                    CSQ_values = [i if i!='' else '-' for i in item.split('|')]
                    for k,v in zip(CSQ_keys, CSQ_values):
                        CSQ_dict[k] = v
                    allele_num = int(CSQ_dict['ALLELE_NUM']) -1
                    csq_multi[allele_num] += item +','
                csq_multi = [i[:-1] for i in csq_multi]
                # output each vari in csq_multi
                for csq in csq_multi:
                    if csq == '':
                        continue
                    # get most severe variant
                    csq_dict = pick_most_severe_CSQ(csq, CSQ_keys, effect_dict)
                    
                    for k, gene_dict in csq_dict.items():
                        for k in gene_dict:
                            if gene_dict[k] == '':
                                gene_dict[k] = '-'
                        allele_num = int(gene_dict['ALLELE_NUM']) -1
                        out_line = [vari_ids[allele_num], pos, ref, gene_dict['Allele'],
                                  gene_dict['Gene'],gene_dict['Feature'],gene_dict['Feature_type'],
                                  gene_dict['Consequence'],gene_dict['cDNA_position'],
                                  gene_dict['CDS_position'],gene_dict['Protein_position'],
                                  gene_dict['Amino_acids'],gene_dict['Codons'],
                                  gene_dict['Existing_variation'], gene_dict['LoF'],
                                  gene_dict['CADD_PHRED']]
                        extra_field = '' # last column in tab format
                        for f in extra_fields:
                            if gene_dict[f] != '-':
                                extra_field = extra_field+f+'='+gene_dict[f]+';'
                        out_line.append(extra_field[:-1])
                        out.write('\t'.join(out_line) + '\n')    
    
                    
# path = '/lustre/workspace/projects/IB/P19_014_UKBB_300K_Vep/data'
vcfs = natsorted(glob.glob(vcf_path + '/*.vcf.gz'))
for vcf in vcfs[start:end]:
    tab = tab_path + '/' + vcf.split('/')[-1][:-6] + 'tsv'
    if not os.path.exists(tab):
        vep_vcf2tab_worst(vcf,tab)
        os.system('gzip -f {t}'.format(t=tab))
