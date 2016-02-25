# Sequence Assembler for CS362 Project 3
# Authors: Charles Eyermann and Sam Hinh

from simulator import *
from graph import *

#print "hi"

class Assembler:

	def __init__(self,reads):
		self.reads = reads
		self.G = {}
		self.N = {}

	def generate_kmers(self,read,k):
		for nucleotide in range(len(read)+1-k):
			yield read[nucleotide:nucleotide+k]

	def build_DBG(self,k):	
		# for each read in set of all reads
		for read in self.reads:
			# for each kmer in set of kmers per read 
			for kmer in self.generate_kmers(read,k):
				# link each left and right k-1mer together
				lkmer = kmer[:-1]
				rkmer = kmer[1:]

				if lkmer not in self.G:
					self.G[lkmer] = Node(lkmer)
					#self.G[lkmer].neighbors.append(Node(rkmer))
				else:
					self.G[lkmer].neighbors.append(Node(rkmer))

				if rkmer not in self.G:
					self.G[rkmer] = Node(rkmer)
					#self.G[rkmer].neighbors.append(Node(lkmer))
				else:
					self.G[rkmer].neighbors.append(Node(lkmer))
		return self.G
				
	def dot_file_generator(self, dbg):
		output = ""
		output += "digraph {\n"
		
		for k in dbg.iterkeys():
			for v in dbg[k].get_neighbors():
				output += "%s -> %s;\n" % (k,v.name)
		output += "}"
		return output


if __name__ == "__main__":
	sequences = Simulator("sample.fasta.txt", 3, 20, 0.2)
	reads = sequences.generate_reads()
	a = Assembler(reads)
	dbg = a.build_DBG(13)
	dot = a.dot_file_generator(dbg)
	with open("out.dot", "w") as f:
		f.write(dot)
		f.close()



