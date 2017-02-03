#!/usr/bin/env python 

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

fasta_file = 'polished_assembly.fasta'
fasta_out = 'selected_contigs.fasta'

excluded = \
		['tig00000003',\
		'tig00000020',\
		'tig00000027',\
		'tig00000049',\
		'tig00000053',\
		'tig00000156',\
		'tig00000163']

genome_length = 0
selected_records = []
for seq_record in SeqIO.parse(fasta_file, "fasta"):
	if not seq_record.id in excluded:
		selected_records.append(seq_record)
		genome_length+=len(seq_record.seq)
		print(seq_record.id+'%s\tlen=%d' % (seq_record.id, len(seq_record.seq)))
		SeqIO.write(seq_record, seq_record.id + '.fasta', "fasta")

SeqIO.write(selected_records, fasta_out, "fasta") # SeqIO.write can parse SeqRecord object or list of SeqRecord objects into file.
print('genome total size = %d' % genome_length) 

