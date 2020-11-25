'''
Created on 2020/07/04 by Shangzhong.Li@pfizer.com
This copies bam files and extract fastq file from bam files
'''
import os
import glob,sys

#--------------- get fq files by natural index ---------------------
# start = sys.argv[1]
# end = sys.argv[2]

# bam_path = '/hpc/grid/wip_drm_targetsciences/projects/p014_RepeatExpansion_ASPERA/BAM'
# out_path = '/lustre/scratch/lis262/RE/fq'
# bam_files = sorted(glob.glob(bam_path+'/*.bam'))

# for bam in bam_files[int(start):int(end)]:
#     prefix = out_path + '/' + bam.split('/')[-1][:-4]
#     fq1 = prefix + '_1.fq.gz'
#     fq2 = prefix + '_2.fq.gz'
#     cmd = 'bsub -q medium -o {log} \"gatk SamToFastq \
#          --INPUT {bam} --FASTQ {fq1} --SECOND_END_FASTQ \
#           {fq2}\"'.format(log=prefix+'.log',bam=bam,fq1=fq1,fq2=fq2)
#     os.system(cmd)
#     # print(bam)


#-------------- get fq files by samples barcode index ---------------
bam_path = '/hpc/grid/wip_drm_targetsciences/projects/p014_RepeatExpansion_ASPERA/BAM'
out_path = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/RE/batch19'
if not os.path.exists(out_path):
    os.mkdir(out_path)

file_index_fn = out_path + '/../file_index.txt'
with open(file_index_fn) as in_f:
    for line in in_f:
        idx = line.strip()
        try:  
            bam = glob.glob(bam_path + '/*_{index}.bam'.format(index=idx))[0]
        except:
            continue
        prefix = out_path + '/' + bam.split('/')[-1][:-4].split('_')[-1]
        fq1 = prefix + '_1.fq.gz'
        fq2 = prefix + '_2.fq.gz'  
        cmd = 'bsub -q medium -o {log} \"gatk SamToFastq \
             --INPUT {bam} --FASTQ {fq1} --SECOND_END_FASTQ \
              {fq2}\"'.format(log=prefix+'.log',bam=bam,fq1=fq1,fq2=fq2)
        os.system(cmd)
        # print(cmd)

