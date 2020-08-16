# Created by Shangzhong Li on 2020/08/10
# After imputation from Michigan imputation server, we get the vcf file, to do the downstream analysis we want to transfer vcf back to plink format or other table format, if we recode using "A" mode, then we will need to provide recode allele to make sure the dosage is counted on the minor allele. This code is to prepare the recode allele file

# get the file for recode-allele
import glob, sys
import pandas as pd

# ipt_vcf_path = '/hpc/grid/wip_drm_targetsciences/projects/TCGA/SNP-Array/IMPUTE/'
ipt_vcf_path = sys.argv[1]
vcf_files = sorted(glob.glob(ipt_vcf_path + '/chr*.info.gz'))
for vcf in vcf_files:
    df = pd.read_csv(vcf,sep='\t',header=0,compression='gzip')
    df['count'] = df.apply(lambda row:row['REF(0)'] \
                if row['ALT_Frq'] != row['MAF'] else row['ALT(1)'],axis=1)
    df[['SNP','count']].to_csv(ipt_vcf_path+'/'+vcf.split('/')[-1][:-7]+\
                'ale.tsv',header=None,index=False,sep='\t')