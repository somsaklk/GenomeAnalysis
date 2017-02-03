#!/usr/bin/env python

table1 = 'sum_25AUG_OT_all-DL1.tsv'
table2 = 'sum_OT_PacBIO1-25AUG_OT_all.tsv'

f = open(table1, 'r')
genes1 = f.readlines()
f.close()

g = open(table2, 'r')
genes2 = g.readlines()
g.close()

header1 = "\t".join(genes2[0].split('\t')[:-2])
header2 = "\t".join(genes1[0].split('\t')[6:12])
header3 = "\t".join(genes2[0].split('\t')[-2:])

w = open('sum.out','w')
w.write(header1 + '\t' + header2 + '\t'+ header3)
del genes1[0]
del genes2[0]

dict_OT = {}

for line1 in genes1:
	features = line1.split('\t')
	ot_gene_old = features[0]
	if ot_gene_old != "No match":
		dict_OT[ot_gene_old] = features[6:12]

for line2 in genes2:
	#~ line2 = line2.replace('\n','')
	features = line2.split('\t')
	ot_gene_old = features[6]
	if ot_gene_old in dict_OT.keys():
		w.write("\t".join(features[:12]) + '\t' + "\t".join(dict_OT[ot_gene_old]) + '\t' + "\t".join(features[-2:]))
	else:
		w.write("\t".join(features[:12]) + '\t' + "\t".join(["No match"]*6) + '\t\t\n')
	
w.close()

			
	
