#!/usr/bin/env python
# Extract all genes from GFF file derived from MAKER run output.
# Collect blast and est2genome match from GFF derived from MAKER run output.
# Also create FASTA file of genes according to gff feature file.
# somsak.lik@biotec.or.th
import re 

def pos_to_gene(pos, gene_dict):
	if pos in gene_dict.keys():
		gene_name = gene_dict[pos]
	else:
		gene_name = 'No match'
	return gene_name

def rvco(seq):
	bases = 'ATGCatgc'
	comp_bases = 'TACGtacg'
	basecode = '12345678'
	for x, y in zip(bases, basecode):
		seq = seq.replace(x,y)
	for y, x in zip(basecode, comp_bases):
		seq = seq.replace(y,x)
	return seq[::-1]

gffile = 'OT_read_PacBIO8.all.gff'
f = open(gffile, 'r')
gff = f.readlines()
f.close()
gene_nametag = 'Name'

# Collect gene sequences from sequence FASTA file.
fastafile = 'OT_read_PacBIO.fasta'
g = open(fastafile, 'r')
fasta = g.read().split('>')
g.close()
del fasta[0]
contig_seq = {}
for each in fasta:
	each = each.replace('\r','').split('\n')
	contig_id = each[0].split()[0]
	seq = "".join(each[1:])
	contig_seq[contig_id] = seq.upper()

outfile = 'OT_PacBIO8'
w = open(outfile+'.out', 'w')
w.write(outfile + '_id\tcontig\tstrand\tstart_position\tend_position\tsize\n')

outfasta = open(outfile+'_gene.fa', 'w')

list_chromosome = [] #list contigs which contain genes
list_genes = []
gene_start_pos = {}
est_match = {}
blastn_match = {}
blastx_match = {}

for line in gff:
	features = line.split('\t')
	if len(features) == 9:
		if features[2] == 'gene':
			gene_id = features[8].split(gene_nametag+'=')[1].split(';')[0].replace('\n','')
			contig = features[0]
			strand = features[6].replace('+','forward').replace('-','reverse')
			start_pos = features[3]
			end_pos = features[4]
			gene_size = int(end_pos) - int(start_pos) + 1
			gene_seq = contig_seq[contig][int(start_pos)-1:int(end_pos)]
            # If gene is on a reverse strand, reverse_complement it.
			if strand == 'reverse':
				gene_seq = rvco(gene_seq)

			w.write(gene_id+'\t'+contig+'\t'+strand+'\t'+start_pos+'\t'\
				+end_pos+'\t'+str(gene_size)+'\n')
			outfasta.write('>'+gene_id+' '+contig+' '+strand+' '+start_pos+'-'\
				+end_pos+'\n')
			for n in range(0, len(gene_seq), 60):
				outfasta.write(gene_seq[n:n+60]+'\n')
			list_chromosome.append(contig)
			list_genes.append(gene_id)
			gene_start_pos[start_pos] = gene_id
		elif features[2] == 'match_part':
			start_pos = features[3]
			match_name = line.split('Target=')[1].split(' ')[0]
			if features[1] == 'est2genome':
				est_match[start_pos] = match_name
			elif features[1] == 'blastn':
				blastn_match[start_pos] = match_name
			elif features[1] == 'blastx':
				blastx_match[start_pos] = match_name
w.close()
outfasta.close()

###########################################
# Count total gene and distribution in contigs
set_chromosome = set(list_chromosome)
set_genes = set(list_genes)
print "Unique genes = " + str(len(set_genes))
for each in set_chromosome:
	print each, list_chromosome.count(each)

###########################################
# Check against EST and Blast results with model EST provided when running MAKER
outgenematch = open(outfile+'_gene_match.out', 'w')
outgenematch.write(outfile+"_ID\test2genome\tblastn\tblastx\n")
for each_pos in gene_start_pos:
	outgenematch.write("\t".join([gene_start_pos[each_pos], pos_to_gene(each_pos, est_match),\
	pos_to_gene(each_pos, blastn_match), pos_to_gene(each_pos, blastx_match)])+'\n')
outgenematch.close()
