# Sequence Read Simulator for CS362 Project 3
# Authors: Charles Eyermann and Sam Hinh

import sys, os, random, argparse

class Simulator:
	'''
	Class to encapsulate the Sequence Read Simulator. Currently the module
	will introduce errors at a specified rate. Future builds will include
	profiles to simulate the biases of specific sequencing machines.
	'''

	def __init__(self, input_file,coverage,read_length,error_rate):
		self.raw_data = ""
		self.header = ""
		self.coverage = coverage # coverage =  NL/G
		self.read_length = read_length
		self.error_rate = error_rate

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
		'''Helper function to produce meaningful mutation'''
		n = ['A','C','G','T']
		n.remove(c)
		return random.choice(n)


	def generate_random_nucleotide(self):
		''' Helper function to generate (semi)random nucleotides'''
		n = ['A','C','G','T']
		return random.choice(n)


	def mangle_data(self, sequence):
		'''Driver function for mutating sequences'''
		#change error rate to an int from 1-100 to make it easier for python's random number generator
		error_percent = 100*self.error_rate

		for i in range(len(sequence)): #for nucleotide in sequence: s[0] to s[19] i range 1 to 19
			# roll the dice
			error = random.randint(1,100)
			# check to see if mutation occured
			if error <= error_percent:

				# edge case if mutation occurs in first index of string
				if i == 0:
					sequence = self.mutate(sequence[i]) + sequence[1:]
					self.mutation_count += 1

				# edge case if mutation occurs in second index of string
				elif i == 1:
					sequence = sequence[0] + self.mutate(sequence[1]) + sequence[2:]
					self.mutation_count += 1

				# edge case if mutation occurs in second to last index of string
				elif i == (len(sequence) - 2):
					sequence = sequence[:-3] + self.mutate(sequence[i]) + sequence[-1:]
					self.mutation_count += 1

				# edge case if mutation occurs in last index of string
				elif i == (len(sequence) - 1):
					sequence[:i] + self.mutate(sequence[i])
					self.mutation_count += 1

				# handle non-edge cases
				else:
					if len(sequence) == self.read_length:
						sequence = sequence[:i] + self.mutate(sequence[i]) + sequence[i+1:]
						self.mutation_count += 1
		return sequence


	def generate_reads(self):
		'''Generate reads at (semi)random intervals along the source genome/sequence'''
		index_range = self.G - self.read_length
		for i in range(self.N):
			read_start_index = random.randint(0, index_range)
			read_end_index = read_start_index + self.read_length
			seq = self.raw_data[read_start_index:read_end_index]
			mutated_seq = self.mangle_data(seq)
			self.reads.append(mutated_seq)
		return self.reads


	def write_file(self, reads):
		'''
		Simply write the reads out to a file with a header giving
		information about the parameters used during read simulation.
		'''
		with open("output_reads.txt", "w") as f:
			h = ">" + self.header[1:30] + "| coverage: " + str(self.coverage) + "| read length: " + str(self.read_length) + "| error rate: " + str(self.error_rate) + "|\n"
			f.write(h)
			for row in self.reads:
				f.write(row + "\n")


if __name__ == '__main__':
	'''Main function to encapsulate argparse functionality'''
	#BEGIN ARGPARSE
	parser = argparse.ArgumentParser(description='Sequence Read Simulator.\r\n Written by Charles Eyermann and Sam Hinh for CS362 Winter 2016')
	parser.add_argument("fasta_sequence_file", type=argparse.FileType('r'), help="The input file should be in FASTA format")
	parser.add_argument("coverage", type=int, help='How many times the entire genome should be sequenced')
	parser.add_argument("read_length", type=int, help='The length of each read generated')
	parser.add_argument("error_rate", type=float, help='Error rate between 0 and 1')
	parser.add_argument('-v', '--verbose', help="Print out some more info about the program at runtime", action="store_true")
	clargs = parser.parse_args()
	#END ARGPARSE

	#BEGIN MAIN
	if clargs.verbose:
		print "Generating simulated read data with following parameters: "
		print "\tCoverage:", clargs.coverage
		print "\tRead length:", clargs.read_length
		print "\tError rate:", clargs.error_rate

	sequences = Simulator(clargs.fasta_sequence_file, clargs.coverage, clargs.read_length, clargs.error_rate)
	reads = sequences.generate_reads()
	if clargs.verbose:
		print "File written to: 'output_reads.txt'"
	sequences.write_file(reads)
