#!/usr/bin/env python
# Translate genbank refID xm_xxxxx into traditional gene symbol using
# gff file downloaded from NCBI

gffile = "DL1.gff"

f = open(gffile, 'r')
gff = f.readlines()
f.close()

id_symbol = {}
gb_id = {}

readingline = 1 

for line in gff:
	line = line.replace('\n','')
	features = line.split('\t')
	# Check reading line for debug
	#~ print readingline
	#~ readingline+=1
	
	if len(features) == 9:
		if features[2] == 'gene':
			primary_id = features[8].split("GeneID:")[1].split(';')[0]
			symbol_name = features[8].split("Name=")[1].split(';')[0]
			id_symbol[primary_id] = symbol_name
		elif features[2] == 'exon':
			primary_id = features[8].split("GeneID:")[1].split(',')[0]
			if "Genbank:" in features[8]:
				genbank_id = features[8].split("Genbank:")[1].split(';')[0]
			else:
				genbank_id = primary_id
			gb_id[genbank_id] = primary_id
	
outfile = gffile+"_name_table.tsv"
w = open(outfile, 'w')

for each in gb_id.keys():
	primary_id = gb_id[each]
	if primary_id in id_symbol.keys():
		symbol_name = id_symbol[primary_id]
	w.write(each+'\t'+symbol_name+'\n')

w.close()
