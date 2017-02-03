#!/usr/bin/env python
# summarise Blast result run on and downloaded from NCBI website
# exported as Hit table (csv)
import sys
import re
from sys import stdout
from Bio import Entrez

def entrez_to_symbol(genbankid, gene_prefix):
	Entrez.email = "somsak.lik@biotec.or.th"
	handle = Entrez.esearch(db="gene", term=genbankid)
	record = Entrez.read(handle)
	primary_id = ",".join(record["IdList"])
	summary = Entrez.efetch(db="gene", id=primary_id)
	summary_record = summary.read()
	name = re.search('"'+gene_prefix+'(.+?)"', summary_record)
	if name:
		name_found = name.group(1)
		return name_prefix+name_found
	else:
		return genbankid
		
######
# Create dictionary of ref Genbank ID translated to gene symbol ID

gffile = "DL1.gff"

f = open(gffile, 'r')
gff = f.readlines()
f.close()

id_symbol = {}
gb_id = {}

for line in gff:
	line = line.replace('\n','')
	features = line.split('\t')
	
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

######
blastfile = 'blastnweb_OT_PacBIO4_v_DL1.csv'
f = open(blastfile,'r')
handle = f.readlines()
f.close()

name_prefix = "HPODL"
gene_match = {}
for line in handle:
	cols = line.split(',')
	if len(cols) == 14:
		query = cols[2]
		evalue = float(cols[-2])
		match = cols[3]
		eset = 10000 # set initial e-value for placeholder, It can be set to arbitrarily high. It does not matter anyway.
		
		# show running record on screen
		stdout.write("\rRunning: %s" % query)
		stdout.flush()
		if not query in gene_match.keys():
			eset = evalue
			if match in gb_id.keys():
				primary_id = gb_id[match]
				if primary_id in id_symbol.keys():
					match_symbol = id_symbol[primary_id]
				else:
					match_symbol = entrez_to_symbol(match, name_prefix)
			else: 
				match_symbol = entrez_to_symbol(match, name_prefix)

			gene_match[query] = match_symbol+'\t'+str(evalue)
		
		else:
			eset = gene_match[query].split('\t')[1]
			eset = float(eset)
			if evalue < eset:
				if match in gb_id.keys():
					primary_id = gb_id[match]
					if primary_id in id_symbol.keys():
						match_symbol = id_symbol[primary_id]
					else:
						match_symbol = entrez_to_symbol(match, name_prefix)
				else: 
					match_symbol = entrez_to_symbol(match, name_prefix)
				
				gene_match[query] = match_symbol+'\t'+str(evalue)
	
print '\n\nTotal matches =', len(gene_match)

w = open(blastfile+'_sum.out', 'w')
for each in gene_match:
	w.write(each+'\t'+gene_match[each]+'\n')
w.close()
