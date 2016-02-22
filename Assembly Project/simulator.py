# Sequence Read simulator for CS362 Project 3
# Authors: Charles Eyermann and Sam Hinh


import numpy as np
import sys, os, collections

class Simulator:

	def __init__(self, input_file):
		self.raw_data = ""
		self.header = ""

		with open(input_file, 'rU') as f:
			for row in f:
				if row[0] == ">":
					self.header += row
				else:
					self.raw_data += row
			




sequence = Simulator("sample.fasta.txt")
print sequence.header


