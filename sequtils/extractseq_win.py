#!/usr/bin/python3
#Description: Extract a sequence from FASTA file, using start and end positions.
## This script uses only standard Python module. No Biopython dependencies required.
#Author: Somsak Likhitrattanapisal <somsak.lik@biotec.or.th>
#Last updated: 2017/01/15

import os

def rvco(seq):
	bases = 'ATGCatgc'
	comp_bases = 'TACGtacg'
	basecode = '12345678'
	for x, y in zip(bases, basecode):
		seq = seq.replace(x,y)
	for y, x in zip(basecode, comp_bases):
		seq = seq.replace(y,x)
	return seq[::-1]

fastalist = [f for f in os.listdir(r'.') if f.endswith('.fa') or f.endswith('.fasta') or f.endswith('.fna')]
fastalist = sorted(fastalist)

if len(fastalist) == 1:
	fastafile = fastalist[0]
	print("> %s is selected." % fastafile)
elif len(fastalist) == 0:
	print("No FASTA files detected in this directory")
	pass
else:
	c = 1
	for file in fastalist:
		print("(%d): %s" % (c, file))
		c+=1
	get_file = int(input("Please select file number: ")) - 1
	fastafile = fastalist[get_file]
	print("> %s is selected.\n" % fastafile)
	
f = open(fastafile, 'r')
handle = f.read().replace(r'\r','').split('>')
f.close()
del handle[0]

run = True
while run:
	seq_ids = []
	seq_descriptions = []
	seqs = []
	i = 1
	for record in handle:
		record = record.split('\n')
		seq_id = record[0].split()[0]
		seq_ids.append(seq_id)
		seq_descriptions.append(record[0])
		seq = "".join(record[1:])
		seqs.append(seq)
		if len(record[0]) > 70:
			print('('+str(i) + '): ' + record[0][:70] + '...\tlen='+str(len(seq)))
		else:
			print('('+str(i) + '): ' + record[0]+'\tlen='+str(len(seq)))
		i+=1	

	if len(seqs) == 1:
		selected = 0
	else:
		selected = int(input("Please select record number: ")) - 1
		
	if 0 <= selected < len(seqs):
		print('> The record "' + seq_ids[selected] + '" is selected.\n')
		seq = seqs[selected]
		try:
			start = int(input("Start position: "))
		except (SyntaxError, ValueError):
			start = 1
		try:
			end = int(input("End position: "))
		except (SyntaxError, ValueError):
			end = len(seqs[selected])
		if 0 < start < end <= len(seq):
			print("> Extract %s sequence region from %d to %d\n" % (seq_ids[selected], start, end)) 
			seq_extract = seq[start-1:end]
			try:
				rvco_opt = input("Reverse complement? (y/N): ")
				rvco_opt = rvco_opt.upper()
			except (SyntaxError, ValueError, NameError):
				rvco_opt = 'N'
			if rvco_opt == 'Y' or rvco_opt == 'YES':
				output_file = "".join(fastafile.split('.')[:-1])+'_rvco_extracted.fasta'
				seq_extract = rvco(seq_extract)
				g = open(output_file, 'w')
				g.write('>'+seq_descriptions[selected] + ', extracted region (reverse complement) ' + str(start) + '...' + str(end) + '\n')
				for n in range(0, len(seq_extract), 60):
					g.write(seq_extract[n:n+60]+'\n')
				g.close()
			else:
				output_file = "".join(fastafile.split('.')[:-1])+'_extracted.fasta'
				g = open(output_file, 'w')
				g.write('>'+seq_descriptions[selected] + ', extracted region ' + str(start) + '...' + str(end) + '\n')
				for n in range(0, len(seq_extract), 60):
					g.write(seq_extract[n:n+60]+'\n')
				g.close()		
			print('\n> OPERATION SUCCESSFUL!\n> The output is written to ' + '"' + output_file + '".')
			term_sig = input("(Press Enter to close)")
			run = False	
		else:
			print("Selected positions invalid. Try again\n")
			pass
	else:
		print("Record number invalid. Try again\n")
		pass
