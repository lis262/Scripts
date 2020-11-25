# merge overlapped bed regions. If overlapped, get the widest region.
from intervaltree import Interval, IntervalTree
from collections import defaultdict

fn = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/GRCh38_WES.bed'
# build dict {chr:IntervalTree(Interval(s, e))}
chr_tree = defaultdict(IntervalTree)
with open(fn) as in_f:
    n = 0
    for line in in_f:
        chrom, s, e = line.strip().split('\t')
        chr_tree[chrom].add(Interval(int(s),int(e)))       
# merge overlaps
for k in chr_tree:
    chr_tree[k].merge_overlaps()
# output to file
from natsort import natsorted
out_fn = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/publicDB/gatk_homo38_index/interval_WES.bed'
with open(out_fn, 'w') as out:
    for chrom in natsorted(chr_tree.keys()):
        for i in chr_tree[chrom]:
            out.write('\t'.join([chrom, str(i.begin), str(i.end)]) + '\n')