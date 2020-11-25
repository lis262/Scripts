'''
Created by Shangzhong.Li@pfizer.com on 2020/09/24.
this is to parse gnomad vcf file and add the allele frequency
and CADD score

'''
import gzip,re
import tabix
gnomad_fn = '/hpc/grid/wip_drm_targetsciences/users/shangl02/RefLib/Human/gnomad/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz'
tb = tabix.open(gnomad_fn)
plof_fn = '/hpc/grid/wip_drm_targetsciences/projects/gnomAD/gnomad_v2/gnomad.v2.1.1.all_lofs.txt.gz'
out_fn = '/hpc/grid/wip_drm_targetsciences/projects/gnomAD/gnomad_v2/gnomad.v2.1.1.all_lofs.af.txt.gz'
cadd_snp = '/hpc/grid/hgcb/workspace/projects/CADDv1_6/GRCh37/whole_genome_SNVs.tsv.gz'
cadd_indel = '/hpc/grid/hgcb/workspace/projects/CADDv1_6/GRCh37/gnomad.genomes.r2.1.1.indel.tsv.gz'

cadd_snp_tb = tabix.open(cadd_snp)
cadd_indel_tb = tabix.open(cadd_indel)

def get_freq(info, key):
    try:
        res = re.search('(?<=;{k}=).+?(?=;)'.format(k=key),info).group(0)
    except:
        res = '-'
    return res

def get_cadd_phred(tb_snp, tb_indel, chrom, pos, ref, alt):
    res = '-'
    if len(ref) == len(alt):
        records = tb_snp.query(chrom,pos-1,pos)
    else:
        records = tb_indel.query(chrom,pos-1,pos)
    for record in records:
        if record[2] == ref and record[3] == alt:
            res = record[5]
            break
    return res
        

add_columns = ['rsid','CADD_PHRED','AF','AF_afr','AF_amr', \
          'AF_asj','AF_eas','AF_fin','AF_nfe', \
                'AF_oth','AF_sas','AF_popmax','popmax']

with  gzip.open(out_fn, 'wb') as out_f, gzip.open(plof_fn,'rt') as in_f:
    head = in_f.readline().strip().split('\t') + add_columns
    out_head = '\t'.join(head) + '\n'
    out_f.write(out_head.encode('utf-8'))
    for line in in_f:
        items = line.strip().split('\t')
        chrom,pos,ref,alt = items[:4]
        pos = int(pos)
        records = tb.query(chrom,pos-1,pos)
        # define the inital values
        rsid = '-'
        afs = ['-'] * (len(add_columns) - 2)
        cadd_phred = '-'
        for record in records:
            if ref == record[3] and alt == record[4]:
                rsid = record[2]
                afs = [get_freq(record[7],k) for k in add_columns[2:]]
                break
            else:
                continue
        cadd_phred = get_cadd_phred(cadd_snp_tb, cadd_indel_tb, chrom, pos, ref, alt)
        out_line = '\t'.join(items + [rsid,cadd_phred] + afs) + '\n'
        # print(out_line)
        out_f.write(out_line.encode('utf-8'))