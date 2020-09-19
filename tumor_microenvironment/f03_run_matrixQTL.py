import os
cancer_types = ['ACC','BLCA','BRCA','CESC','CHOL','COAD',
                'DLBC','ESCA','GBM','HNSC','KICH','KIRC',
                'KIRP','LAML','LGG','LIHC','LUAD','LUSC',
                'MESO','OV','PAAD','PCPG','PRAD','READ',
                'SARC','SKCM','STAD','TGCT','THCA','THYM',
                'UCEC','UCS','UVM']

# cancer_types = ['DLBC','GBM','LAML','LGG','PCPG','SARC','TGCT','UCS','UVM']
# cancer_types = ['BRCA']
path = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor'
matrixQTL = '~/Code/Scripts/tumor_microenvironment/utils/m04_get_input_matrixQTL.py'
mode = 'log'
race = 'all'
score = 'cyto'
for cancer in cancer_types:
    # 1. run matrixQTL
    log = path + '/' + cancer + '.log'
    cmd = 'bsub -M 100457280 -q medium -o {log} \"python {code}  -c {c} \
        -m {m} -r {r} -s {s} \"'.format(log=log,code=matrixQTL,c=cancer,m=mode,r=race,s=score)
    # print(cmd)
    os.system(cmd)
