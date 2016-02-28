# Sequence Assembler for CS362 Project 3
# Authors: Charles Eyermann and Sam Hinh

from simulator import *
from graph import *
import argparse

p = argparse.ArgumentParser(description='Global Affine Alignment Algorithm.\r\n Written by Charles Eyermann for CS362 Winter 2016')
p.add_argument('read_file', type=argparse.FileType('r'), help="This should be the file containing your sequence reads")
p.add_argument('kmer_length', type=int, help="Desired k-mer length (integer value).")
clargs = p.parse_args()

class Assembler:

	def __init__(self, reads, k):
		self.reads = []
		self.G = {}
		self.N = {}
		self.header = ""
		self.k = k

		for row in clargs.read_file:
			if row[0] != ">":
				self.reads.append(row)
				#print len(row)

	def generate_kmers(self,read,k):
		for nucleotide in range(len(read)+1-self.k):
			yield read[nucleotide:nucleotide+self.k]

	def build_DBG(self):	
		# for each read in set of all reads
		for read in self.reads:
			# for each kmer in set of kmers per read 
			for kmer in self.generate_kmers(read,self.k):
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
		
		for key in dbg.iterkeys():
			for v in dbg[key].get_neighbors():
				output += "%s -> %s;\n" % (key,v.name)
		output += "}"
		return output


if __name__ == "__main__":
	#sequences = Simulator("test.txt", 5, 8, 0.1)
	#reads = sequences.generate_reads()
	a = Assembler(clargs.read_file, clargs.kmer_length)
	dbg = a.build_DBG()
	dot = a.dot_file_generator(dbg)

	with open("out.dot", "w") as f:
		f.write(dot)
		f.close()



