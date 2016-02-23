# Sequence Read simulator for CS362 Project 3
# Authors: Charles Eyermann and Sam Hinh

import simulator, graph

#print "hi"

class Assembler:

	def __init__(self, reads, k):
		self.reads = reads
		self.k = k

	def __init__(self,reads):
		self.reads = reads
		

	def generate_kmers(self,k):
		kmer_list = []
		for read in self.reads:
			#print read
			for nucleotide in range(len(read)+1-k):
				kmer_list.append(read[nucleotide:nucleotide+k])
		#print len(kmer_list)
		return (read,kmer_list)


sequences = simulator.Simulator("sample.fasta.txt", 3, 10, 0.01)
reads = sequences.generate_reads()
#print len(reads)
a = Assembler(reads)
print len(a.generate_kmers(9)[0]), a.generate_kmers(9)[0]
print len(a.generate_kmers(9)[1])


