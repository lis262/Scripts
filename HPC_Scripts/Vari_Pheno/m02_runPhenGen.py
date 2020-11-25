'''
Created on 2020/07/07 by Shangzhong.Li@pfizer.com
this file run phen_gene analysis

'''
import os
import glob
import argparse
import shutil

parser = argparse.ArgumentParser(description='run phen_gen')
# parser.add_argument('-i','--input',action='store',dest='input',help='input folder')
# parser.add_argument('-f','--phen',action='store',dest='phen_gen',help='phen_gen folder')
parser.add_argument('-n','--name',action='store',dest='name',help='out prefix name')
parser.add_argument('-p','--predictor',action='store',dest='pred',help='predictor, 1: coding, 0: genomic',default='1')
parser.add_argument('-d','--disease',action='store',dest='dise_mode',help='disease mode, 1: recessive, 0: dorminant',default='1')

args = parser.parse_args()

# folder = args.input
# phen_gen = args.phen_gen
name = args.name
pred = args.pred
dise_mode = args.dise_mode

folder = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/NIH_UDP/PheGen'
phen_gen = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/NIH_UDP/Software/Phen_Gene'
vcf = folder + '/m00_VCF/{name}.vcf.gz'.format(name=name)
ped = folder + '/m02_pedigree/{name}.ped'.format(name=name)
hpo = folder + '/m01_HPO/ncr_hpo/{name}.txt'.format(name=name)
out_path = folder + '/p01_phen_gen/' + name
if not os.path.exists(out_path):
	os.mkdir(out_path)

# copy vcf file
shutil.copy(vcf,phen_gen)
cmd = 'gunzip {f}'.format(f=phen_gen+'/'+vcf.split('/')[-1])
os.system(cmd)

# copy hpo file
out_hpo = phen_gen + '/' + hpo.split('/')[-1]
cmd = 'cut -f 3 {hpo} | sort | uniq > {out}'.format(hpo=hpo,out=out_hpo)
os.system(cmd)
# copy pedigree file
shutil.copy(ped, phen_gen)

# run phengene
p_vcf = vcf.split('/')[-1][:-3]
p_ped = ped.split('/')[-1]
p_hpo = hpo.split('/')[-1]

os.chdir(phen_gen)
cmd = ('perl phen-gen.pl input_phenotype={hpo} input_vcf={vcf} input_ped={ped} inheritance={dise} predictor={pred} stringency=1 discard_de_novo=1').format(hpo=p_hpo,vcf=p_vcf,ped=p_ped,dise=dise_mode,pred=pred)
print(cmd)
os.system(cmd)

# move results back to folder
os.remove(p_vcf)
os.remove(p_ped)
os.remove(p_hpo)

vcf_prefix = '.'.join(p_vcf.split('.')[:-1])
ped_prefix = '.'.join(p_ped.split('.')[:-1])
hpo_prefix = '.'.join(p_hpo.split('.')[:-1])

os.remove(p_vcf + '_genlog.txt')
os.remove(p_hpo + '_phenlog.txt')

score_fn = ped_prefix+'_'+vcf_prefix+'_genescores.txt'
os.rename(score_fn, out_path +'/'+name+'.scores.txt')
topgenes = ped_prefix+'_'+vcf_prefix+'_variants_for_topgenes.vcf'
os.rename(topgenes, out_path + '/'+name+'.topgenes.vcf')