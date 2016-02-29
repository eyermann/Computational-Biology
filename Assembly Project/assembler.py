# Sequence Assembler for CS362 Project 3
# Authors: Charles Eyermann and Sam Hinh

from simulator import *
from graph import * 
import argparse


class Assembler:

	def __init__(self, reads, k):
		self.reads = []
		self.G = {}
		self.N = {}
		self.header = ""
		self.k = k

		for row in clargs.read_file:
			if row[0] != ">":
				self.reads.append(row.strip())
				# print len(row)

	def generate_kmers(self,read,k):
		for nucleotide in range(0,len(read)+1-self.k):
			yield read[nucleotide:nucleotide+self.k]

	def build_DBG(self):	
		# for each read in set of all reads
		for read in self.reads:
			#print "read: ", read
			# for each kmer in set of kmers per read 
			for kmer in self.generate_kmers(read,self.k):
				# link each left and right k-1mer together
				lkmer = kmer[:-1]
				rkmer = kmer[1:]

				if lkmer not in self.G:
					self.G[lkmer] = Node(lkmer)

					if rkmer not in self.G[lkmer].neighbors:
						self.G[lkmer].neighbors.append(Node(rkmer))
				else:
					self.G[lkmer].neighbors.append(Node(rkmer))

				if rkmer not in self.G:
					self.G[rkmer] = Node(rkmer)
				else:
					pass

				self.G[lkmer].outdegree += 1
				self.G[rkmer].indegree += 1

		root = [x for x in self.G.iterkeys() if self.G[x].is_head()][0]
		self.populate_parents_dfs(self.G,root)
		return self.G
				
	def dot_file_generator(self, dbg):
		output = ""
		output += "digraph {\n"
		
		for key in dbg.iterkeys():
			for v in dbg[key].get_neighbors():
				output += "%s -> %s;\n" % (key,v.name)
		output += "}"
		return output

	def eulerian_walk(self,dbg):
		walk = []
		#choose random node
		#start = random.choice(self.G.keys())
		start = self.G.iterkeys().next()
		print "Starting node: ", start
		def next(node):
			#print [x.name for x in self.G[node].get_neighbors()]
			#print len(self.G[node].get_neighbors())
			neighbors = self.G[node].get_neighbors()
			while len(neighbors) > 0:
				#print neighbors
				nxt = neighbors.pop()
				next(nxt.name)
			walk.append(node)
		next(start)
		walk = walk[::-1]
		return walk

	def find_leaves(self):
		leaves = []
		for x in self.G.iterkeys():
			if self.G[x].indegree == 1 and self.G[x].outdegree == 0:
				leaves.append(self.G[x])
		return leaves

	def get_cycle(self, dbg):
			path = []
			visited = set()
			
			def visit(node):
				if node in visited:
					return None
				visited.add(node)
				path.append(node)
				
				neighbors = [x for x in self.G[node.name].neighbors]
				for neighbor in neighbors:
					if neighbor in path or visit(neighbor):
						return path
				path.remove(node)
			for x in self.G.iterkeys():
				paths = visit(self.G[x])
			return paths
	
	# def get_parent(self,node):
	# 	#print "node: ", node.name
	# 	parents = []
	# 	for x in self.G.iterkeys():
	# 		#print self.G[x].name
	# 		#print "node:",x, "neighbors: ", [y.name for y in self.G[x].neighbors]
	# 		if node.name in [y.name for y in self.G[x].neighbors]:
	# 			#print y.name, y
	# 			parents.append(self.G[x])
	# 	return parents

	def dfs(self, dbg, start):
		visited, stack = [], [start]
		while stack:
			v = stack.pop()
			if v not in visited:
				visited.append(v)
				q = [item.name for item in dbg[v].neighbors if item not in visited]
				stack.extend(q)
		return visited

	def populate_parents_dfs(self, dbg, start):
		visited, stack = [], [start]
		while stack:
			v = stack.pop()
			if v not in visited:
				visited.append(v)
				q = [item.name for item in dbg[v].neighbors if item not in visited]
				for n in q:
					a.G[n].parents.append(v)
				stack.extend(q)
		return visited

	# def collapse(self): 
	# 		# while the indegree of a node is 1, collapse
	# 		def concat(node):
	# 			#if node.is_collapsible():
	# 			if self.get_parent(node)[0].indegree == 1:
	# 				#print len(self.get_parent(node))
	# 				if len(self.get_parent(node)) == 1:
	# 					parent = self.get_parent(node)[0]
	# 					cur = node.name
	# 					kmer_len = len(parent.name)
	# 					#print "kmer len: ", kmer_len
	# 					#print "current node: ", node.name
	# 					#print "'parent' node: ", parent.name
	# 					if self.get_parent(self.get_parent(node)[0]):
	# 						real_parent = self.get_parent(self.get_parent(node)[0])[0]
	# 						#print "real parent: ", real_parent.name
	# 						if parent.name[1:] == cur[:(kmer_len-1)]:
	# 							conc = parent.name[0] + cur
	# 							#print "concatenation: ", conc
	# 							newnode = self.G[conc] = Node(conc)
	# 							real_parent.neighbors.append(newnode)
	# 							del self.G[parent.name]
	# 							del self.G[node.name]
					
	# 							concat(newnode)

	# 		for leaf in self.find_leaves():
	# 			concat(leaf)	
	def collapse(self):
		# while the indegree of a node is 1, collapse
		def concat(node):
			if node.is_collapsible():
			#if self.get_parent(node)[0].indegree == 1:
				print len(node.get_parents())
				if len(node.get_parents()) == 1:

					parent = self.G[node.get_parents()[0]]
				 	cur = node
				 	kmer_len = len(parent.name)
				 	print kmer_len
				 	print "kmer len: ", kmer_len
				 	print "current node: ", cur.name
				 	print "'parent' node: ", parent.name

				 	if parent.get_parents()[0]:
						real_parent = self.G[parent.get_parents()[0]]
						print "real parent: ", real_parent.name
						if parent.name[1:] == cur.name[:(kmer_len-1)]:
							conc = parent.name[0] + cur.name
							print "concatenation: ", conc
							newnode = self.G[conc] = Node(conc)
							print "parent neighbors: ", parent.get_neighbors()[0].name
							print "real parent neighbors: ", real_parent.get_neighbors()[0].name
							real_parent.neighbors.append(newnode)
							print "neighbors: ", self.G[real_parent.name].neighbors
							#self.G[real_parent.name].neighbors.remove(parent)
							del self.G[parent.name]
							del self.G[node.name]
							#del real_parent
				
				 			concat(newnode)

		for leaf in self.find_leaves():
			concat(leaf)	






if __name__ == "__main__":
	#sequences = Simulator("test.txt", 5, 8, 0.1)
	#reads = sequences.generate_reads()
	p = argparse.ArgumentParser(description='Global Affine Alignment Algorithm.\r\n Written by Charles Eyermann for CS362 Winter 2016')
	p.add_argument('read_file', type=argparse.FileType('r'), help="This should be the file containing your sequence reads")
	p.add_argument('kmer_length', type=int, help="Desired k-mer length (integer value).")
	clargs = p.parse_args()
	
	a = Assembler(clargs.read_file, clargs.kmer_length)
	dbg = a.build_DBG()
	#print a.G['o_e'].parents
	#print "keys: ", [x for x in dbg.iterkeys()]
	#oot = [x for x in a.G.iterkeys() if a.G[x].is_head()][0]
	#print "root: ", root
	#print a.populate_parents_dfs(dbg, root)


	#dot = a.dot_file_generator(dbg)
	#print [x for x in a.G.iterkeys()]

	# for x in a.G.iterkeys():
	# 	print x, a.G[x].get_degree(), a.G[x].indegree, a.G[x].outdegree, a.G[x].is_head(), a.G[x].is_collapsible(), a.G[x].parents

	#print a.eulerian_walk(dbg)
	#to test eulerian walk output
	#superstr = a.eulerian_walk(dbg)
	#superstring = superstr[0] + ''.join(map(lambda x: x[-1], superstr[1:]))
	#print "Eulerian walk results: ", superstring

	a.collapse()
	#print a.dfs(dbg, root)
	#print [x for x in a.G.iterkeys()]
	dot = a.dot_file_generator(dbg)

	with open("out.dot", "w") as f:
		f.write(dot)
		f.close()

