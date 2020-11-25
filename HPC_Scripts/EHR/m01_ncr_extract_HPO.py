'''
Created on July 6th 2020 by Shangzhong.Li@pfizer.com
used to extract HPO from clinical notes using NCR https://github.com/ccmbioinfo/NeuralCR.
'''
import ncrmodel
import argparse

parser = argparse.ArgumentParser(description='extract HPO from EHR')

parser.add_argument('-i','--input',action='store',dest='notes',help='input notes file')
parser.add_argument('-g','--model',action='store',dest='model',help='trained moel path')
parser.add_argument('-o','--output',action='store',dest='hpo',help='output file')


args = parser.parse_args()

EHR_file = args.notes
model_path = args.model
hpo_file = args.hpo

# path = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/NIH_UDP/ncr'
param_dir = model_path + '/model_params'
word_model_file = param_dir + '/pmc_model_new.bin'
model = ncrmodel.NCR.loadfromfile(param_dir, word_model_file)

with open(EHR_file) as in_f:
    text = in_f.readliens()
    text = ' '.join(text)
model.annotate_text(text, 0.8)


# -------------- run in command line
# python NeuralCR/annotate_text.py --params model_params --fasttext model_params/pmc_model_new.bin --input notes/UDP_11559.txt --output annotation/UDP_11559.txt
