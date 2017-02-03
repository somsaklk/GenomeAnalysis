#!/usr/bin/env python
#
# somsak.lik@biotec.or.th

import re 

match_table = 'OT_PacBIO7_gene_match.out'
f = open(match_table,'r')
handle = f.readlines()
f.close()
del handle[0]

w = open(match_table.replace('.out', '_sum.out'), 'w')

for line in handle:
	check = line.count("No match")
	gene = line.split('\t')[0]
	if check < 3:
		match_gene = re.findall(r'OT_(\d+)', line)[0]
		w.write(gene+'\t'+'OT_'+match_gene.zfill(5)+'\n')
w.close()
	
