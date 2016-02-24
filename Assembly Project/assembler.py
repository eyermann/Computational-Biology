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
		big_list = []
		
		for read in self.reads:
			kmer_list = []
			for nucleotide in range(len(read)+1-k):
				kmer_list.append(read[nucleotide:nucleotide+k])
			kmer_info = (read, kmer_list)
			big_list.append(kmer_info)
		#print len(kmer_list)
		return big_list


sequences = simulator.Simulator("sample.fasta.txt", 3, 50, 0.01)
reads = sequences.generate_reads()
a = Assembler(reads)
print len(a.generate_kmers(3)[1][1])#, a.generate_kmers(9)
#print len(a.generate_kmers(9)[1])


