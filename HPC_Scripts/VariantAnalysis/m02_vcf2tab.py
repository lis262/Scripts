'''
This file converts the vcf file from vcf format to table format
'''

import gzip,re
from collections import OrderedDict
import argparse
# vcf_file = vep_path + '/fsgs.chrY.vep.vcf.gz'
# tab_file = tab_path + '/fsgs.chrY.vep.tsv'

# import sys
# vcf_file = sys.argv[1]
# tab_file = sys.argv[2]

parser = argparse.ArgumentParser(description='transfer vcf to table format')

parser.add_argument('-i','--input',action='store',dest='vcf',help='input vcf file')
parser.add_argument('-o','--output',action='store',dest='table',help='output table file')

args = parser.parse_args()

vcf_file = args.vcf
tab_file = args.table



# define scores of variants effects
effect_dict = {
    'transcript_ablation':1.0,
    'splice_acceptor_variant':0.95,
    'splice_donor_variant':0.95,
    'stop_gained':0.95,
    'splice_region_variant':0.94,
    'frameshift_variant':0.94,
    'stop_lost':0.9,
    'start_lost':0.93,
    'transcript_amplification':0.6,
    'inframe_insertion':0.7,
    'inframe_deletion':0.7,
    'missense_variant':0.7,
    'protein_altering_variant':0.7,
    'incomplete_terminal_codon_variant':0.9,
    'start_retained_variant':0.65,
    'stop_retained_variant':0.65,
    'synonymous_variant':0.65,
    'coding_sequence_variant':0.65,
    'mature_miRNA_variant':0.65,
    '5_prime_UTR_variant':0.65,
    '3_prime_UTR_variant':0.65,
    'non_coding_transcript_exon_variant':0.64,
    'intron_variant':0.64,
    'NMD_transcript_variant':0.63,
    'non_coding_transcript_variant':0.63,
    'upstream_gene_variant':0.6,
    'downstream_gene_variant':0.6,
    'TFBS_ablation':0.6,
    'TFBS_amplification':0.6,
    'TF_binding_site_variant':0.6,
    'regulatory_region_ablation':0.6,
    'regulatory_region_amplification':0.6,
    'feature_elongation':0.6,
    'regulatory_region_variant':0.6,
    'feature_truncation':0.6,
    'intergenic_variant':0.6,
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
    csq_records = csq_notes.split(',')
    csq_dict = OrderedDict()
    csq_dict = {k:'' for k in csq_fields}
    effect_idx = csq_fields.index('Consequence')
    for csq in csq_records: # loop each transcript
        csq_values = csq.split('|')
        # for a transcript with &, choose the most severe term
        effect = csq_values[effect_idx]
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
        if effect_dict[csq_values[effect_idx]] > effect_dict[csq_dict['Consequence']]:
            for k,v in zip(csq_fields, csq_values):
                csq_dict[k] = v
    return csq_dict          


def update_prediction_score(csq_dict, key):
    """
    In vep, the sift polyphen, condel, coral merge the scores and description
    together, here we split them if the annotation has it
    """
    if key in csq_dict:
        content = csq_dict.pop(key,None)
        if ('(') in content:
            start = content.index('(')
            end = content.index(')')
            pred = (content[:start]).lower()
            score = content[start+1:end]
        else:
            score = ''
            pred = ''
        
        csq_dict[key+'_score'] = score
        csq_dict[key+'_pred'] = pred
    return csq_dict


def vcf2tab(vcf_file, tab_file):
    with gzip.open(vcf_file,'rt') as in_vcf, gzip.open(tab_file,'wb') as out:
        write_header = True # writ header or not
        for line in in_vcf:
            if line.startswith('##INFO=<ID=CSQ'): # extract CSQ column names
                try:
                    csq_fields = re.search('(?<=Format: ).+?(?=\">)', line).group(0).split('|')
                except:
                    raise 'header does not have CSQ field, please annotate with VEP'
            elif line.startswith('##'):
                continue
            else: # variants lines and column name lines are left and this step
                items = line.strip().split('\t')
                if line.startswith('#CHROM'): # column name line
                    columns = items
                    columns[0] = 'CHROM'
                    samples = columns[9:]  # get sample names in vcf file
                else:
                    # build a dictionary,key is sample name
                    # value is another dict, with format field as the key
                    samples_fmt_dict = OrderedDict()
                    samples_fmt_dict = {k:{} for k in samples}
                    # get each column value
                    chrom,pos,ID,ref,alt,qual,filt,info = items[0:8]
                    fmt = items[8].split(':')
                    # dictionary example {{sp1:{'GT':'0/1','DP':'53'}}}
                    for sp, gt in zip(samples, items[9:]):
                        for f,g in zip(fmt, gt.split(':')):
                            samples_fmt_dict[sp][f] = g
                    # build dictionary for info field
                    info_dict = OrderedDict()
                    for f in info.split(';'):
                        try:
                            k,v = f.split('=')
                            info_dict[k] = v
                        except:
                            info_dict[k] = True
                    # get samples genotype
                    genotypes = [samples_fmt_dict[s]['GT'] for s in samples]
                    depth = [samples_fmt_dict[s]['DP'] for s in samples]
                    try:
                        csq_notes = re.search('(?<=CSQ=).+?(?=$)', info).group(0)
                    except:
                        continue
                    csq_dict = pick_most_severe_CSQ(csq_notes, csq_fields, effect_dict)
                    csq_dict = update_prediction_score(csq_dict, 'SIFT')
                    csq_dict = update_prediction_score(csq_dict, 'PolyPhen')
                    csq_dict = update_prediction_score(csq_dict, 'CAROL')
                    csq_dict = update_prediction_score(csq_dict, 'Condel')
                    # remove some columns if exist
                    drop_keys = ['Allele','INTRON',
                                    'CDS_position','Protein_position',
                                     'ALLELE_NUM']
                    for drop in drop_keys:
                        _ = csq_dict.pop(drop,None)
                    # write to file
                    if write_header: # hasn't write header yet
                        keys = list(csq_dict.keys())

                        header = ['chr','pos','ref','alt'] + keys + \
                                ['gt_' + s for s in samples]
                        out_line = '\t'.join(header) + '\n'
                        out.write(out_line.encode('utf-8'))
                        write_header = False
                    # write data
                    csq_values = [v if v != '' else '-' for k,v in csq_dict.items()]
                    row = [chrom,pos,ref,alt] + csq_values + \
                            [g['GT'] for s,g in samples_fmt_dict.items()]
                    out_line = '\t'.join(row) + '\n'
                    out.write(out_line.encode('utf-8'))

vcf2tab(vcf_file, tab_file)