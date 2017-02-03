#!/usr/bin/env python
# Try to compare genes from 2 genomes based on chromosomal locations
# Use one as base and align the other on top of it.
import operator
from sys import stdout

gfile1 = 'OT_PacBIO7.out' # base gff
gfile2 = 'OT_PacBIO8.out' # ab initio gff
n = 20 # Threshold for position difference
outfile = 'output.tsv'
w = open(outfile, 'w')

f = open(gfile1, 'r')
g1 = f.readlines()
f.close()
del g1[0]

g = open(gfile2, 'r')
g2 = g.readlines()
g.close()
del g2[0]

contig = []

for line in g1:
	features = line.rstrip().split('\t')
	contig.append(features[1])
	
contig = set(contig)

gene_new = 0
new_genes = []

for eachcontig in contig:
	w.write(eachcontig+'\n')
	#~ dictstart1 = {}
	#~ dictend1 = {}
	#~ dictstart2 = {}
	#~ dictend2 = {}
	gene_ticks = []
	inter_ticks = []
	for line in g1:
		features = line.rstrip().split('\t')
		if features[1] == eachcontig:
			start_pos = int(features[3])
			end_pos = int(features[4])
			gene_ticks.append(start_pos)
			gene_ticks.append(end_pos)
			inter_ticks.append(start_pos-1)
			inter_ticks.append(end_pos+1)
	
	if len(gene_ticks) == 0:
		break
	
	# add first start and last end position for intergenic region
	inter_ticks.append(0)
	inter_ticks.append(max(gene_ticks)+int(0.5*max(gene_ticks)))
	gene_ticks = sorted(gene_ticks)
	inter_ticks = sorted(inter_ticks)
	
	gene_region = []
	for i in range(0,len(gene_ticks),2):
		gene_region.append((gene_ticks[i], gene_ticks[i+1]))
	inter_region = []
	for j in range(0,len(inter_ticks),2):
		inter_region.append((inter_ticks[j], inter_ticks[j+1]))
	
	regions = gene_region + inter_region
	regions = sorted(regions)

	#~ break
	for region in regions:
		#~ stdout.write('\r'+eachcontig)
		#~ stdout.flush() 
		#~ stdout.write('\r'+str(region))
		#~ stdout.flush()
		
		section1 = '\t' + str(region[0]) + '\t' + str(region[1])
		found1 = 0
		found_new = False
		while found1 == 0:
			for line1 in g1:
				features = line1.rstrip().split('\t')
				if features[1] == eachcontig:
					start_pos1 = int(features[3])
					end_pos1 = int(features[4])
					if abs(start_pos1 - region[0]) <= n and abs(end_pos1 - region[1]) <= n:
						section1 = "\t".join([features[0], features[3], features[4]])
						found1=1
			found1 = 1
			
		lines_to_section2 = []
		for line2 in g2:
			features = line2.rstrip().split('\t')
			if features[1] == eachcontig:
				start_pos2 = int(features[3])
				end_pos2 = int(features[4])
				print(region)
				if start_pos2 >= region[0] - n and end_pos2 <= region[1] - n:
					section2 = "\t".join([features[0], features[3], features[4]])
					lines_to_section2.append(section2)
		
		if len(lines_to_section2) > 0:
			for sec in lines_to_section2:
				w.write(section1 + '\t' + sec + '\n')
			# Check if there are new genes from g2
			if section1 == '\t' + str(region[0]) + '\t' + str(region[1]):
				if len(lines_to_section2) > 2:
					found_new = True
			else:
				if len(lines_to_section2) > 1:
					found_new = True
				
		else:
			w.write(section1 + '\t\t\t\n')
					
		if found_new:
			gene_new += 1
			new_genes.extend(lines_to_section2)

w.close()
print(gene_new)
for gene in new_genes:
	print(gene)
