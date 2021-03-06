"""
Created by Shangzhong.Li@pfizer on 2020/05
This generates commands to run vep on the UKBB vcf files,
output is in vcf format, use minimal and allele number to split the multi-alt allele.
also add rsid
"""
import os,glob
#===================== change these 3 parameters ========================
# path to the data
path = '/lustre/scratch/lis262/anno/lei'
# file that you want to save command to
command_file = path + '/cmd.sh'
# path to the container
container = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/container/vep_loftee_101.0_gnomad_pLI.sif'



#==================== don't change anything below this line ===============
prefix = '/media'
loftee_path = prefix + '/hpc/grid/wip_drm_targetsciences/projects/VEP_plugin/loftee_grch37'
pub_db_path = prefix + '/hpc/grid/wip_drm_targetsciences/users/shangzhong/publicDB'
vep_dir = prefix + '/hpc/grid/wip_drm_targetsciences/projects/VEP_plugin/vep_db_101_GRCH37'
ref_fa = prefix + '/hpc/grid/wip_drm_targetsciences/users/shangzhong/publicDB/hg19.fa'
CADD_SNP = prefix + '/hpc/grid/hgcb/workspace/projects/CADDv1_6/GRCh37/whole_genome_SNVs.tsv.gz'
CADD_INDEL = prefix + '/hpc/grid/hgcb/workspace/projects/CADDv1_6/GRCh37/gnomad.genomes.r2.1.1.indel.tsv.gz'
gnomADe = prefix + '/hpc/grid/wip_drm_targetsciences/users/shangl02/RefLib/Human/gnomad/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz'
gnomADg = prefix + '/hpc/grid/wip_drm_targetsciences/users/shangl02/RefLib/Human/gnomad/vcf/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz'

# vep_path = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/publicDB/vep_db_97_GRCH38'
vcf_files = glob.glob(path + '/*.vcf.gz')
out_path = path + '/vep_anno_vcf'
if not os.path.exists(out_path):
    os.mkdir(out_path)

with open(command_file,'w') as out_f:
    for vcf in vcf_files:
        vcf = prefix + vcf
        out_vcf = out_path + '/' + vcf.split('/')[-1] + '.gz'

        cmd = ('singularity run -B /:/media {container} \
                vep -i {in_f} \
               --plugin LoF,loftee_path:{loftee},human_ancestor_fa:{loftee}/human_ancestor_fa.gz,gerp_bigwig:{loftee}/gerp_conservation_scores.homo_sapiens.GRCh38.bw,conservation_file:{loftee}/phylocsf_gerp.sql \
               --dir_plugins {loftee} \
               --custom {ge},gnomADe,vcf,exact,0,AF_afr,AF_amr,AF_asj,AF_eas,AF_fin,AF_nfe,AF_oth \
               --custom {gg},gnomADg,vcf,exact,0,AF_afr,AF_amr,AF_asj,AF_eas,AF_fin,AF_nfe,AF_oth \
               --plugin Carol \
               --plugin Condel,/opt/.vep/Plugins/config/Condel/config,b \
               --plugin ExACpLI,/opt/gnomad.v2.1.1.oe_lof.by_gene.txt \
               --plugin LoFtool,/opt/.vep/Plugins/LoFtool_scores.txt \
               --plugin CADD,{CADD_SNP},{CADD_INDEL},\
               -o {out} --cache --force_overwrite --buffer_size 10000 \
               --species homo_sapiens --assembly GRCh37 \
               --dir {vep_dir} --offline --variant_class --fork 1 --hgvs -e \
               --fa {ref_fa} --minimal  --compress_output gzip \
               --allele_number --check_existing --vcf').format(
                in_f=vcf,
                container=container,
                ge=gnomADe,gg=gnomADg,
                loftee = loftee_path,
                CADD_SNP = CADD_SNP,
                CADD_INDEL = CADD_INDEL,
                vep_dir = vep_dir,
                out = prefix + '/' + out_path + '/' + vcf.split('/')[-1] + '.gz',
                ref_fa = ref_fa)
        pre_cmd = 'bsub -q short -o {log}.log.txt \"{cmd}\"'.format(
                  cmd=cmd,log=vcf.split('/')[-1])
        if not os.path.exists(out_vcf):
          out_f.write(pre_cmd + '\n')
    #     assert False
