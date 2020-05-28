import os, glob

path = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor'
file_path = path + '/files'
file_lists = glob.glob(file_path + '/split*')



for f in file_lists:
	cmd = ('bsub -M 41457280 -o log_{log} -q medium "bcftools merge -l {l} -O z -o {p}/{log}.vcf.gz"').format(log=f.split('/')[-1],l=f,p=path)
	os.system(cmd)
	# print(cmd)

