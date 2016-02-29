# Sequence Read simulator for CS362 Project 3
# Authors: Charles Eyermann and Sam Hinh

import numpy as np 
import sys, os, collections, random, argparse



class Simulator:

	def __init__(self, input_file,coverage,read_length,error_rate):
		self.raw_data = ""
		self.header = ""
		self.coverage = coverage # coverage =  NL/G
		self.read_length = read_length
		self.error_rate = error_rate
		
		#with open(input_file, 'rU') as f:
		for row in clargs.fasta_sequence_file:
			if row[0] == ">":
				self.header += row.strip()
			else:
				self.raw_data += row.strip()
	
		self.G = len(self.raw_data)
		self.N = (self.coverage * self.G)/self.read_length
		self.NL = self.N * self.read_length
		self.reads = []
		self.mutation_count = 0

	def mutate(self, c):
		n = ['A','C','G','T']
		n.remove(c)
		return random.choice(n)

	def generate_random_nucleotide(self):
		n = ['A','C','G','T']
		return random.choice(n)

	def mangle_data(self, sequence):
		#change error rate to an int from 1-100 to make it easier for python's random number generator
		error_percent = 100*self.error_rate

		for i in range(len(sequence)): #for nucleotide in sequence: s[0] to s[19] i range 1 to 19
			# roll the dice
			error = random.randint(1,100)
			# check to see if mutation occured
			if error <= error_percent:
				#print "mutation occured at position i: ", i
				
				# edge case if mutation occurs in first index of string
				if i == 0:
					#print "<"+self.mutate(sequence[i])+">"+sequence[1:]
					sequence = self.mutate(sequence[i]) + sequence[1:]
					self.mutation_count += 1
				
				# edge case if mutation occurs in second index of string
				elif i == 1:
					#print sequence[0]+"<"+self.mutate(sequence[1])+">"+sequence[2:]
					sequence = sequence[0] + self.mutate(sequence[1]) + sequence[2:]
					self.mutation_count += 1

				# edge case if mutation occurs in second to last index of string
				elif i == (len(sequence) - 2):
					#print sequence[:-3]+"<"+self.mutate(sequence[i])+">"+sequence[-1:]
					sequence = sequence[:-3] + self.mutate(sequence[i]) + sequence[-1:]
					self.mutation_count += 1

				# edge case if mutation occurs in last index of string
				elif i == (len(sequence) - 1):
					#print sequence[:i]+"<"+self.mutate(sequence[i])+">"
					sequence[:i] + self.mutate(sequence[i])
					self.mutation_count += 1

				# handle non-edge cases
				else:
					if len(sequence) == self.read_length:
						#print sequence[:i]+"<"+self.mutate(sequence[i])+">"+sequence[i+1:]	
						sequence = sequence[:i] + self.mutate(sequence[i]) + sequence[i+1:]
						self.mutation_count += 1
						
		return sequence

	def generate_reads(self):
		index_range = self.G - self.read_length
		#print index_range
		for i in range(self.N):
			read_start_index = random.randint(0, index_range)
			read_end_index = read_start_index + self.read_length
			seq = self.raw_data[read_start_index:read_end_index]
			mutated_seq = self.mangle_data(seq)
			self.reads.append(mutated_seq)
		
		return self.reads


	def write_file(self, reads):
		with open("easy_output_reads.txt", "w") as f:
			h = ">" + self.header[1:30] + "| coverage: " + str(self.coverage) + "| read length: " + str(self.read_length) + "| error rate: " + str(self.error_rate) + "|\n"
			f.write(h)
			for row in self.reads:
				f.write(row + "\n")


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Sequence Read Simulator.\r\n Written by Charles Eyermann and Sam Hinh for CS362 Winter 2016')
	parser.add_argument("fasta_sequence_file", type=argparse.FileType('r'), help="The input file should be in FASTA format")
	parser.add_argument("coverage", type=int, help='How many times the entire genome should be sequenced')
	parser.add_argument("read_length", type=int, help='The length of each read generated')
	parser.add_argument("error_rate", type=float, help='Error rate between 0 and 1')
	clargs = parser.parse_args()
	#sequences = Simulator("easy_data_set.txt", 10, 50, 0.01)
	sequences = Simulator(clargs.fasta_sequence_file, clargs.coverage, clargs.read_length, clargs.error_rate)
	# print " " + "_"*78 + ""
	# print "|" + " "*78 + "|"
	# print "|" + " "*78 + "|"
	# print "|" + "_"*78 + "|"
	#print sequences.header
	reads = sequences.generate_reads()
	#print reads
	sequences.write_file(reads)


