import glob,os
from functools import reduce


METAL = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/software/generic-metal/metal'
# path = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor/'
path = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor'
work_path = path + '/bakup_matrixQTL'

def overlap_cells(a,b):
    return list(set(a).intersection(set(b)))


cell_qtl_path = work_path + '/f04_cell_gwas'
meta_path = work_path + '/f06_meta_qtl'
if not os.path.exists(meta_path):
    os.mkdir(meta_path)

studies = ['KIRC','KIRP']
name = 'Kidney'
race = 'eur'
mode = 'norm'
score = 'cyto'  

# 1. get overlap cell types
cell_types = []
for study in studies:
    qtl_fns = glob.glob(cell_qtl_path+'/{c}/*_{c}_{r}_{s}_{m}.txt.gz'.format(
                    c=study,r=race,s=score,m=mode))
    qtl_fns = sorted(qtl_fns)
    cells = []
    for qtl in qtl_fns:
        idx = qtl.split('/')[-1].split('_').index(study)
        cells.append('_'.join(qtl.split('/')[-1].split('_')[:idx]))
    cell_types.append(cells)
cell_types = reduce(overlap_cells, cell_types)

# 2. run METAL for each cell type
analyze_path = meta_path + '/' + name
if not os.path.exists(analyze_path):
    os.mkdir(analyze_path)
for cell in cell_types:
    metal_files = [sorted(glob.glob(cell_qtl_path+'/{cancer}/{cell}_{cancer}_{r}_{s}_{m}*.gz'.format(
                  cell=cell,cancer=s,r=race,s=score,m=mode)))[0] for s in studies]
    config = analyze_path + '/config_{c}.txt'.format(c=cell)
    with open(config,'w') as out_f:
        config_dic = {'MARKER':'marker',
                      'ALLELE':'REF ALT',
                     'SEPARATOR':'TAB',
                     'SCHEME':'STDERR',
                     'EFFECT':'beta',
                      'STDERR':'std',
                      'PVALUE':'p-value',
                     'OUTFILE':analyze_path+'/'+cell+'_'+name+'_'+race+ ' .txt'}
        for k,v in config_dic.items():
            out_f.write(k+' '+v+'\n')
        for qtl in metal_files:
            out_f.write('PROCESS ' + qtl+'\n')
        out_f.write('ANALYZE HETEROGENEITY\nQUIT')
    # run in bsub
    cmd = 'bsub -q short -o {o}.bsub.log {metal} {config} {log}'.format(
        o=analyze_path+'/'+cell+'_'+name+'_'+race,
        metal=METAL,config=config,log=analyze_path+'/log.txt')
    os.system(cmd)
# for cell in cell_types:
#     meta_res_fn = analyze_path+'/'+cell+'_'+name+'_'+race+'1.txt'
#     meta_df = pd.read_csv(meta_res_fn,sep='\t',header=0)
#     meta_df = meta_df.sort_values('P-value')
#     meta_df.to_csv(meta_res_fn,sep='\t',index=False)