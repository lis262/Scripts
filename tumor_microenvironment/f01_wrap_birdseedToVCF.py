'''
this file wraps function to run BIRDSEEDToVCF for our internal TCGA data
'''
import glob,os,sys
import pandas as pd
pd.set_option('display.max_columns', None)

start = int(sys.argv[1])
end = int(sys.argv[2])
#------------------------  define parameters
ref_fa = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/publicDB/hg19.fa'
snp_anno_fn = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor/GenomeWideSNP_6.na30.annot.hg19.csv.pickle'
path = '/hpc/grid/wip_drm_targetsciences/projects/TCGA/'
birdseed_code = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor/m01_BirdseedToVCF3.py'
header = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor/header.txt'
# derived parameters
snp_array_path = path + '/SNP-Array'
magetab_path = snp_array_path + '/MAGETAB'
sdrf = glob.glob(magetab_path + '/broad*/*.sdrf.txt')
meta_data = path + '/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81/TCGA-CDR-SupplementalTableS1.csv'


meta_df = pd.read_csv(meta_data,header=0)
samples = meta_df['bcr_patient_barcode'].tolist()

#------------------------ build sample and file name mapping dictionary get files that are not in the metadata file
not_in_seed = snp_array_path + '/array_with_no_meta.txt'
fn_sample_dict = {}
sample_fn_dict = {}
with open(not_in_seed,'w') as out_h:
    for sd in sdrf:
        try:
            df = pd.read_csv(sd,sep='\t',header=0,encoding='latin1')
            dic = df.set_index('Derived Array Data Matrix File.1')\
                      ['Comment [TCGA Barcode]'].to_dict()
            fn_sample_dict.update(dic)
            dic = df.set_index('Comment [TCGA Barcode]') \
                      ['Derived Array Data Matrix File.1'].to_dict()
            sample_fn_dict.update(dic)
            
            birdseed_samples = df['Comment [TCGA Barcode]'].tolist()
            birdseed_samples = ['-'.join(s.split('-')[:3]) for s in birdseed_samples]
            non_overlap_samples = set(birdseed_samples) - set(samples)
            for birdseed_fn in df['Derived Array Data Matrix File.1']:
                s = '-'.join(fn_sample_dict[birdseed_fn].split('-')[:3])
                if s in non_overlap_samples:
                    out_h.write('\t'.join([s, birdseed_fn,  \
                            sd.split('/')[-1]])+'\n')
        except:
            print(sd)
            assert False

#---------------- convert files -----------------------
birdseed_path = snp_array_path + '/BIRDSEED'
birdseed_files = glob.glob(birdseed_path + '/*/*.txt')

for birdseed in birdseed_files[start:end]:
    # 1. transfer birdseed to vcf file format
    out_sp = fn_sample_dict[birdseed.split('/')[-1]]
    out_vcf = os.path.dirname(birdseed) + '/' + out_sp + '.vcf'
    reheader_vcf = out_vcf[:-3] + 'rehead.vcf.gz'
    array_sp = birdseed.split('/')[-1].split('.')[0]
    if not os.path.exists(reheader_vcf):
        # 1. birdseed to vcf
        cmd = 'python {py} --birdseed {array} \
              --snp_annotation_file {anno} \
              --array_sample {a_s} \
              --output_vcf {vcf} \
              --vcf_sample {sp} \
              --fasta {ref_fa} --add_chr && \
              (grep ^\"#\" {vcf}; grep -v ^\"#\" {vcf} \
              | sort -k1,1 -k2,2n) | bgzip > {vcf}.gz && \
              tabix {vcf}.gz && rm {vcf}'.format(
                        py=birdseed_code,array=birdseed,anno=snp_anno_fn,
                        a_s=array_sp, vcf=out_vcf,
                        sp=out_sp.split('/')[-1],ref_fa=ref_fa)
        # os.system(cmd)
        # 2. reheader of vcf so the vcf file can be merged
        reheader_vcf = out_vcf[:-3] + 'rehead.vcf.gz'
        cmd = 'bcftools view -h {vcf} | cat {h} - | \
               bcftools reheader -h - {vcf} > {r_vcf} && \
               tabix {r_vcf} && \
               rm {vcf} && rm {vcf}.tbi'.format(
              h=header,vcf=out_vcf+'.gz',r_vcf=reheader_vcf)
        os.system(cmd)
        # print(cmd)


