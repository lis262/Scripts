import gzip, re, os, glob
import argparse

parser = argparse.ArgumentParser(description='change vep vcf output to tab')
parser.add_argument('-i','--vcf',action='store',dest='input',help='input vcf file')
parser.add_argument('-o','--out',action='store',dest='output',help='output vcf file')

args = parser.parse_args()
in_vcf = args.input
out_vcf = args.output
# vcf = args.input
# tab = args.output

def search_af(anno, tag):
    try:
        res = re.search('(?<=;{t}=).+?(?=;)'.format(t=tag),anno).group(0)
    except:
        res = '-'
    return res

# fn = '/hpc/grid/wip_drm_targetsciences/projects/gnomAD/gnomad_v2/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz'
# out_fn = '/hpc/grid/wip_drm_targetsciences/projects/gnomAD/gnomad_v2/simple/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.gz'
pre = 'AF_'
with gzip.open(in_vcf,'rt') as in_f, gzip.open(out_vcf,'wb') as out:
    for line in in_f:
        if line.startswith('#'):
            out.write(line.encode('utf-8'))
        else:
            items = line.strip().split('\t')
            anno = items[7]
            AF = search_af(anno,'AF')
            AF_afr = search_af(anno,pre+'afr')
            AF_amr = search_af(anno,pre+'amr')
            AF_asj = search_af(anno,pre+'asj')
            AF_eas = search_af(anno,pre+'eas')
            AF_fin = search_af(anno,pre+'fin')
            AF_nfe = search_af(anno,pre+'nfe')
            AF_oth = search_af(anno,pre+'oth')
            items[7] = ';'.join(['AF='+AF,'AF_afr='+AF_afr,'AF_amr='+AF_amr,
                                'AF_asj='+AF_asj,'AF_eas='+AF_eas,'AF_fin='+AF_fin,
                                'AF_nfe='+AF_nfe,'AF_oth='+AF_oth])
            out_line = '\t'.join(items) + '\n'
            out.write(out_line.encode('utf-8'))

# transfer to bgzip format
cmd = ('gunzip -c {f} | bgzip > {o}').format(f=out_vcf,o=out_vcf[:3]+'.bgz')
os.system(cmd)
os.rename(out_vcf[:3]+'.bgz', out_vcf[:3]+'.gz')