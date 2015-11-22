#! /usr/bin/python

# Written by Simon Watson and Matt Cotten, Wellcome Trust Sanger Institute
#
# Takes a pair of read files in FASTQ format, and BLASTs them against
# a BLAST database. Any chimeric reads are split, with the larger fragment
# kept in the file and the second framgment placed in a separate file (to
# keep the numbers in forward and reverse the same)
#
# "num_to_parse" changes the number of reads that are BLASTed together. Lower
# to reduce memory consumption.
#
# To use this script, you need the following:
# 1) FASTQ module from QUASR: https://github.com/simonjwatson/QUASR
# 2) A BLAST database to BLAST the reads against
# 3) The BLASTN binary


import sys, os
sys.path.append('/Users/sw10/Dropbox/Sanger/QUASR/QUASR_v6.09/modules/') # 1)
import fastq
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

def blast_reads(blast_string, reads, outfh, outExtra):
	blast_db = '/Users/sw10/Dropbox/Sanger/blastdb/ebola/Zaire_ebolavirus_KM034562' # 2)
	blast_binary = '/Applications/ncbi-blast-2.2.29+/bin/blastn' # 3)
	xml_outfile = '/tmp/test.xml'
	evalue = 0.01 
	cline = NcbiblastnCommandline(cmd=blast_binary, out=xml_outfile, outfmt=5, query="-", db=blast_db, evalue=evalue, max_target_seqs=1, num_threads=1)
	stdout, stderr = cline(blast_string)

	with open(xml_outfile, 'r') as blast_handle:
		blast_records = NCBIXML.parse(blast_handle)
		for blast_record in blast_records:
			name = blast_record.query
			for alignment in blast_record.alignments:
				count = 1
				hits = {}
				top_hitBool = True
				for hsp in alignment.hsps:
					seq = reads[name].sequence[hsp.query_start:hsp.query_end]
					qual = reads[name].quality[hsp.query_start:hsp.query_end]
					if hsp.sbjct_start > hsp.sbjct_end:
						tmp1 = [seq[i] for i in range(len(seq)-1,-1,-1)]
						seq = ''.join(tmp1)
						tmp2 = [qual[i] for i in range(len(qual)-1,-1,-1)]
						qual = ''.join(tmp2)
		
					header = '%s:%d' % (name, count)
					hits[header] = (seq, qual)
					count += 1
			for head, seq_tuple in sorted(hits, key=lambda head: len(hits[head][0]), reverse=True):
				if top_hitBool == True:
					outfh.write("@%s\n%s+\n%s\n" % (head, seq_tuple[0], seq_tuple[1]))
					top_hitBool = False
				else:
					outExtra.write("@%s\n%s+\n%s\n" % (head, seq_tuple[0], seq_tuple[1]))
				
	os.remove(xml_outfile)


if len(sys.argv) != 5:
	print('chimericReadSplitter.py <forward.fq> <reverse.fq> <outfile.fq> <num_to_parse>')
	sys.exit(0)

infile = sys.argv[1]
pairfile = sys.argv[2]
outprefix = sys.argv[3]
max_at_once = int(sys.argv[4])

outfileF1 = "%s.1.fq" % outprefix
outfileR1 = "%s.2.fq" % outprefix
outfileF2 = "%s.1b.fq" % outprefix
outfileR2 = "%s.2b.fq" % outprefix

counter = 0
readsF = {}
readsR = {}
blast_string = ''
with open(infile, 'r') as infh, open(pairfile, 'r') as revfh, open(outfileF1, 'w') as outf1, open(outfileR1, 'w') as outr1, open(outfileF2, 'w') as outf2, open(outfileR2, 'w') as outr2:
	reverse_reads = fastq.fastq_iterator(revfh)
	for headerF, sequenceF, qualityF in fastq.fastq_iterator(infh):
		headerR, sequenceR, qualityR = next(reverse_reads)
		blast_stringF += '>%s\n%s\n' % (headerF, sequenceF)
		blast_stringR += '>%s\n%s\n' % (headerR, sequenceR)
		readsF[header] = fastq.FastqRecord(headerF, sequenceF, qualityF)
		readsR[header] = fastq.FastqRecord(headerR, sequenceR, qualityR)
		counter += 1
		if counter == max_at_once:
			blast_reads(blast_stringF, readsF, outf1, outf2)
			blast_reads(blast_stringR, readsR, outr1, outr2)
			counter = 0
			reads = {}
			blast_string = ''
	blast_reads(blast_stringF, readsF, outf1, outf2)
	blast_reads(blast_stringR, readsR, outr1, outr2)
