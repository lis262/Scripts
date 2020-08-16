import os
cancer_types = ['ACC','BLCA','BRCA','CESC','CHOL','COAD',
                'DLBC','ESCA','GBM','HNSC','KICH','KIRC',
                'KIRP','LAML','LGG','LIHC','LUAD','LUSC',
                'MESO','OV','PAAD','PCPG','PRAD','READ',
                'SARC','SKCM','STAD','TGCT','THCA','THYM',
                'UCEC','UCS','UVM']

# cancer_types = ['CESC','DLBC','GBM','KIRC','LAML','LGG','OV','PCPG','PRAD','SARC','SKCM','STAD','TGCT','UCEC','UCS','UVM']
# cancer_types = ['BRCA']
path = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor'
for cancer in cancer_types:
    # 1. prepare matrixQTL input
    log = path + '/' + cancer + '.log'
    cmd = 'bsub -M 100457280 -q short -o {log} \"python ~/Code/Scripts/tumor_microenvironment/utils/m04_get_input_matrixQTL.py -c {c}\"'.format(log=log,c=cancer)
    # print(cmd)
    os.system(cmd)
