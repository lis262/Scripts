'''
Created on 2020/08/24 by Shangzhong.Li@pfizer.com
this file run exomiser analysis
'''
import argparse
import shutil

sw = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/software/Exomiser/exomiser-rest-prioritiser-12.1.0.jar'
path = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/NIH_UDP/Exomiser'


parser = argparse.ArgumentParser(description='run exomiser')

parser.add_argument('-n','--name',action='store',dest='name',help='out prefix name')
parser.add_argument('-p','--predictor',action='store',dest='pred',help='predictor, 1: coding, 0: genomic',default='1')
parser.add_argument('-d','--disease',action='store',dest='dise_mode',help='disease mode, AR: recessive, AD: dorminant',default='AD')

args = parser.parse_args()
name = args.name
pred = args.pred
dise_mode = args.dise_mode


folder = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/NIH_UDP/Exomiser'
vcf = folder + '/m00_VCF/{name}.vcf.gz'.format(name=name)
ped = folder + '/m02_pedigree/{name}.ped'.format(name=name)
hpo = folder + '/m01_hpo/ncr_hpo/{name}.txt'.format(name=name)

hpo_ids = []
with open(hpo) as in_f:
    for line in in_f:
        item = line.strip().split('\t')
        hpo_ids.append(item[2])
hpo_ids = list(set(hpo_ids))

cmd = ('java -Xms2g -Xmx4g -jar {sw} --prioritiser hiphive -I {mode} -F 1 --hpo-ids {hpo} -v {vcf}').format(sw=sw,mode=dise_mode,hpo=','.join(hpo_ids),vcf=vcf)
print(cmd)
print(hpo_ids)
