import gzip,sys

in_vcf = sys.argv[1]
out_vcf = sys.argv[2]

prev = ''
with gzip.open(in_vcf,'rt') as in_f, gzip.open(out_vcf,'wb') as out:
	for line in in_f:
		if line.startswith('#'):
			out.write(line.encode('utf-8'))
		else:
			item = line.strip().split('\t')
			item[2] = '_'.join([item[0],item[1],item[3],item[4]])
			if item[2] == prev:
				continue
			else:
				prev = item[2]
			out_line = '\t'.join(item) + '\n'
			out.write(out_line.encode('utf-8'))

