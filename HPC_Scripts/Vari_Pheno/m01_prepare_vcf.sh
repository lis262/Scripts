# merge vcf
chain=/hpc/grid/wip_drm_targetsciences/users/shangzhong/publicDB/hg38ToHg19.chain.gz
ref=/hpc/grid/wip_drm_targetsciences/users/shangzhong/publicDB/hg19.fa
# hg38 to hg19
CrossMap.py vcf $chain $1 $ref $2
# sort vcf and bgzip and tabix
cat $2 | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' | sed 's/^chr//' - | bgzip > $3

