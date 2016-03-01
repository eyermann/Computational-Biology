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

		roots = [x for x in self.G.iterkeys() if self.G[x].is_head()]
		for root in roots:
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

	def dot_file_generator_list(self, l):
		output = ""
		output += "graph {\n"
		
		for item in l:
			output += "%s;\n" % (item)
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

	def find_leaves(self,dbg):
		leaves = []
		for x in dbg.iterkeys():
			if dbg[x].indegree == 1 and dbg[x].outdegree == 0:
				leaves.append(dbg[x])
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

	def bfs(self, dbg, start):		
		visited, q = [], [start]
		while q:
			v = q.pop(0)
			if self.G[v] not in visited:
				visited.append(self.G[v])
				p = [item.name for item in self.G[v].neighbors if item not in visited]
				q.extend(p)
		return visited


	def dfs(self, dbg, start):
		visited, stack = [], [start]
		while stack:
			v = stack.pop()
			if self.G[v] not in visited:
				visited.append(self.G[v])
				q = [item.name for item in self.G[v].neighbors if item not in visited]
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
					#print self.G[v]
					a.G[n].parents.append(self.G[v])
				stack.extend(q)
		return visited

	def concat(self,node):
		contigs = []

		#print "is collapsible: ", node.is_collapsible(), 
		#print len(node.get_parents())
		#print node.get_parents()[0].get_parents()
		if node.get_parents():
			#print node.is_collapsible() and node.get_parents()[0].is_collapsible(),
			if node.is_collapsible() and node.get_parents()[0].is_collapsible():
			#if self.get_parent(node)[0].indegree == 1:
				#print "prev: ", node.get_parents()[0]
				if len(node.get_parents()) == 1:
					#print self.G[node.get_parents()[0]]

					parent = node.get_parents()[0]
					#print parent
					cur = node
					kmer_len = len(parent.name)

					if parent.get_parents():
						real_parent = parent.get_parents()[0]
						#print "real parent: ", real_parent.name, real_parent
						#print parent.name[1:] == cur.name[:(kmer_len-1)],
						if parent.name[1:] == cur.name[:(kmer_len-1)]:
							conc = parent.name[0] + cur.name
							#print "! Concatenation: ", conc
							newnode = self.G[conc] = Node(conc)
							real_parent.neighbors.append(newnode)
							newnode.parents.append(real_parent)
							newnode.indegree = len(newnode.parents)
							#print "real parent: ", real_parent.name, [x.name for x in real_parent.neighbors]#[0].name
							#self.G[real_parent.name].neighbors.remove(real_parent.neighbors[0])
							self.G[real_parent.name].neighbors = [x for x in self.G[real_parent.name].neighbors if x.name != parent.name]
							#print "real parent: ", real_parent.name, [x.name for x in real_parent.neighbors]
							del self.G[parent.name]
							del self.G[node.name]
							self.concat(newnode)

			elif node.is_collapsible():
				parent = node.get_parents()[0]
				kmer_len = len(parent.name)
				cur = node
				if parent.name[1:] == cur.name[:(kmer_len-1)]:
					conc = parent.name[0] + cur.name
					#print conc
					newnode = self.G[conc] = Node(conc)	
					#print newnode.name

					#del self.G[parent.name]
					#del self.G[node.name]
					#print type(conc)
					

					


if __name__ == "__main__":
	#sequences = Simulator("test.txt", 5, 8, 0.1)
	#reads = sequences.generate_reads()
	p = argparse.ArgumentParser(description='Global Affine Alignment Algorithm.\r\n Written by Charles Eyermann and Sam Hinh for CS362 Winter 2016')
	p.add_argument('read_file', type=argparse.FileType('r'), help="This should be the file containing your sequence reads")
	p.add_argument('kmer_length', type=int, help="Desired k-mer length (integer value).")
	clargs = p.parse_args()

	if clargs.kmer_length < 2:
		p.print_help()
		print "\n"
		print "FATAL: kmer length must be 2 or more! - suggested to use kmer length of 4 at minimum."
		sys.exit(1)
	
	a = Assembler(clargs.read_file, clargs.kmer_length)
	dbg = a.build_DBG()

	# to test eulerian walk output
	#superstr = a.eulerian_walk(dbg)
	#print len(superstr)
	#print superstr
	#superstring = superstr[0] + ''.join(map(lambda x: x[-1], superstr[1:]))
	#print "Eulerian walk results: ", superstring


	# Get all leaves in tree, collapse them!
	#print len(a.find_leaves())
	#for leaf in a.find_leaves(dbg):
	#	a.concat(leaf)

	# find all roots in tree
	roots = [x for x in dbg.iterkeys() if dbg[x].is_head()]
	print "number of roots found: ", len(roots)
	#print max(len(a.bfs(dbg,root)) for root in roots)

	# find longest depth-first path in graph
	max_root = None
	max_len = 0
	for root in roots:
		z = len(a.dfs(dbg,root))
		if z > max_len:
			max_root = root
			max_len = z
	print "root with longest path: ", max_root, " has length: ", max_len
	dfs = [x.name for x in a.dfs(dbg,max_root)]
	bfs = [x.name for x in a.bfs(dbg,max_root)]
	print "depth first search results: ", dfs
	print "breadth first search results: ", bfs
	
	#superstring = bfs[0] + ''.join(map(lambda x: x[-1], bfs[1:]))
	#print superstring

	print "roots found through bfs: ", [x for x in bfs if a.G[x].is_head()]
	print "leaves found through bfs: ", [x for x in bfs if a.G[x].is_leaf()]
	
	print "roots found through bfs: ", [x for x in dfs if a.G[x].is_head()]
	print "leaves found through dfs: ", [x for x in dfs if a.G[x].is_leaf()]
	
	# count = 0
	# for x in dbg.iterkeys():
	# 	if dbg[x].is_head():
	# 		count +=1
	# print count
	#maxlen = max([len(x) for x in a.G.iterkeys()])
	#contigs = [x for x in a.G.iterkeys() if len(x) == maxlen]

	dot = a.dot_file_generator(dbg)

	with open("out.dot", "w") as f:
		f.write(dot)
		f.close()

