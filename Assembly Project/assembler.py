# Sequence Assembler for CS362 Project 3
# Authors: Charles Eyermann and Sam Hinh

from simulator import *
from graph import * 
import argparse

class Assembler:

	def __init__(self, data, k = None):
		if k is None:
			self.reads = []
			self.G = {}
			self.keys = []
			self.header = ""
			self.G = self.dot_file_to_graph(data)


		else:
			self.reads = []
			self.G = {}
			self.N = {}
			self.header = ""
			self.k = k

			for row in clargs.read_file:
				if row[0] != ">":
					self.reads.append(row.strip())
					# print len(row)
			self.build_DBG()

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
						self.G[lkmer].neighbors.append(Node(rkmer).name)
				else:
					self.G[lkmer].neighbors.append(Node(rkmer).name)

				if rkmer not in self.G:
					self.G[rkmer] = Node(rkmer)
				else:
					pass

				self.G[lkmer].outdegree += 1
				self.G[rkmer].indegree += 1

		roots = [x for x in self.G.iterkeys() if self.G[x].is_head()]
		for root in roots:
			self.populate_parents_dfs(self.G,root)
		self.keys = self.G.keys()
		for item in self.G.iterkeys():
			self.G[item].unique_neighbors = len(set([x for x in self.G[item].neighbors]))
			self.G[item].unique_parents = len(set([x for x in self.G[item].get_parents()]))
		return self.G
				
	def dot_file_generator(self, dbg):
		output = ""
		output += "digraph {\n"
		
		for key in dbg.iterkeys():
			for v in dbg[key].get_neighbors():
				output += "%s -> %s;\n" % (key,v)
		output += "}"
		return output

	def weighted_dot_file_generator(self, dbg):
		wd = {}
		output = ""
		output += "digraph {\n"
		
		for key in dbg.iterkeys():
			for v in dbg[key].get_neighbors():
				k = "%s -> %s" % (key,v)
				if k not in wd:
					wd[k] = 1
				else:
					wd[k] += 1
				#output += "%s -> %s;\n" % (key,v.name)
		for z in wd.iterkeys():
			y = wd[z]
			output += "%s [ label=\"%s\" ];\n" % (z,y)
		output += "}"
		return output

	def weighted_graph_converter(self, dbg):
		wd = {}

		for key in dbg.iterkeys():
			for v in dbg[key].neighbors:
				k = "%s -> %s" % (key,v)
				if k not in wd:
					wd[k] = 1
				else:
					wd[k] += 1
		return wd

	def dot_file_generator_list(self, l):
		output = ""
		output += "graph {\n"
		
		for item in l:
			output += "%s;\n" % (item)
		output += "}"
		return output

	def dot_file_to_graph(self, dotfile):
		
		with open(dotfile, "rU") as f:
			for row in f:
				row = row.strip()[:-1]
				if row:
					l = row.split(" -> ") 
					lkmer,rkmer = l[0],l[-1]
					if lkmer not in self.G:
						self.G[lkmer] = Node(lkmer)

						if rkmer not in self.G[lkmer].neighbors:
							self.G[lkmer].neighbors.append(Node(rkmer).name)
					else:
						self.G[lkmer].neighbors.append(Node(rkmer).name)

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


	def eulerian_walk(self):
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
		return [x for x in a.G.iterkeys() if a.G[x].is_leaf()]

	def find_main_path(self,dbg):
		return

	def get_cycle(self):
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

	def eulerianess(self):
		return len([a.G[node].name for node in self.G.iterkeys() if a.G[node].get_degree() != 0])

	def populate_parents_dfs(self, dbg, start):
		visited, stack = [], [start]
		while stack:
			v = stack.pop()
			if v not in visited:
				visited.append(v)
				q = [item for item in dbg[v].neighbors if item not in visited]
				for n in q:
					#print self.G[v]
					self.G[n].parents.append(self.G[v].name)
				stack.extend(q)
		return visited

	def merge_nodes(self,n1,n2):
		#assume we are collapsing upwards i.e ATTGC + TTGCG = ATTGCG
		newname = n1.name[0] + n2.name
		#print newname
		newnode = Node(newname)
		newnode.parents = n1.parents
		newnode.neighbors = n2.neighbors
		newnode.indegree = n1.indegree
		newnode.outdegree = n2.outdegree
		a.G[newname] = newnode
		del a.G[n1.name]
		del a.G[n2.name]

	def collapse_driver(self):
		if self.keys:
			for k,v in self.G.items():
				print "\n"
				print "current list of nodes in dict: ", self.keys
				cur, nbor_list = k, [x for x in v.neighbors]
				print "node: '", cur, "' neighbors: ", nbor_list
				if nbor_list:
					try:

						# if all nodes in the neighbors list for the current node are the same and they overlap, 
						# we can collapse the two nodes 

						print a.G[nbor_list[0]].unique_neighbors
						#print len(a.G[cur].get_parents())
						next_node = nbor_list[0]
						if (nbor_list.count(next_node) == len(nbor_list)
						and (cur[1:] == next_node[:(len(cur)-1)])
						and (a.G[cur].unique_parents <= 1)
						and (a.G[next_node].unique_parents <= 1)
						and (a.G[cur].unique_neighbors <= 1)
						and (a.G[next_node].unique_neighbors <= 1)):


							print "Homogenous neighbor list AND",
							print "Should be identical (overlap): <" , cur[1:], "==", nbor_list[0][:(len(cur)-1)],">", 
							print ", so merge these nodes: ", cur, "and ", nbor_list[0]
							#print 

							n1,n2 = a.G[cur], a.G[nbor_list[0]]
							#print n1.name,n2.name
							#self.merge_nodes(a.G[cur],a.G[nbor_list[0]])
							newname = n1.name[0] + n2.name
							#print "LOLOLOL", newname
							
							newnode = Node(newname)
							newnode.parents = n1.parents
							newnode.neighbors = n2.neighbors
							newnode.unique_parents = len(set(n1.parents))
							newnode.unique_neighbors = len(set(n2.neighbors))
							#print "DATA DUMP: ",
							#newnode.data_dump()
							#newnode.indegree = n1.indegree
							#newnode.outdegree = n2.outdegree
							a.G[newname] = newnode
							self.keys.append(newname)	
							
							del a.G[n1.name]
							del a.G[n2.name]
							#print n1,n2
							self.keys.remove(n1.name)
							self.keys.remove(n2.name)
							
							self.collapse_driver()
						else:
							print "If False: ", nbor_list.count(nbor_list[0]) == len(nbor_list), "| AND",
							print "NOT identical: <" , cur[1:], nbor_list[0][:(len(cur)-1)],">",
							print "we won't merge these nodes: ", cur, "and ", nbor_list[0]
							continue

					except KeyError,e:

						print "ABORT ABORT: ", str(e)
						#self.collapse_driver()
						break
					
				else: 
					print ""
					#self.collapse_driver()
			


	def concat(self,node):
		contigs = []

		#print "is collapsible: ", node.is_collapsible(), 
		#print len(node.get_parents())
		#print node.get_parents()[0].get_parents()
		if type(node) == str:
			node = self.G[node]
			#print "this is also a node", node
		#else:
			#print "this is a node", node

		if node.get_parents():
			#print node.is_collapsible(), self.G[node.get_parents()[0]].is_collapsible()
			if node.is_collapsible() and self.G[node.get_parents()[0]].is_collapsible():
			#if self.get_parent(node)[0].indegree == 1:
				#print "prev: ", node.get_parents()[0]
				if len(node.get_parents()) == 1:
					#print self.G[node.get_parents()[0]]

					parent = node.get_parents()[0]
					#print parent
					cur = node
					kmer_len = len(parent)

					if self.G[parent].get_parents():
						real_parent = self.G[parent].get_parents()[0]
						#print "real parent: ", real_parent.name, real_parent
						#print parent.name[1:] == cur.name[:(kmer_len-1)],
						if parent[1:] == cur.name[:(kmer_len-1)]:
							conc = parent[0] + cur.name
							#print "! Concatenation: ", conc
							newnode = self.G[conc] = Node(conc)
							self.G[real_parent].neighbors.append(newnode.name)
							newnode.parents.append(real_parent)
							newnode.indegree = len(newnode.parents)
							newnode.name = conc
							#print "real parent: ", real_parent.name, [x.name for x in real_parent.neighbors]#[0].name
							#self.G[real_parent.name].neighbors.remove(real_parent.neighbors[0])
							self.G[real_parent].neighbors = [x for x in self.G[real_parent].neighbors if x != parent]
							#print "real parent: ", real_parent, [x for x in self.G[real_parent].neighbors]
							del self.G[parent]
							del self.G[node.name]
							self.concat(newnode.name)

			elif node.is_collapsible():
				parent = node.get_parents()[0]
				kmer_len = len(parent)
				cur = node
				if parent[1:] == cur.name[:(kmer_len-1)]:
					conc = parent[0] + cur.name
					#print conc
					newnode = self.G[conc] = Node(conc)	
					newnode.name = conc

					#del self.G[parent]
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
	# END ARGPARSE CODE

	# BEGIN MAIN
	a = Assembler(clargs.read_file, clargs.kmer_length)
	# a1 = Assembler("out.dot")
	# print a1.G
	print "Building DeBruijn graph... please wait"
	# just as a reminder, dbg is referencing self.G graph in assembler class. (a.G[key] == dbg[key])
	#a.build_DBG()
	print "Finished building DeBruijn graph!"
	for k,v in a.G.items():
		print k,v.unique_neighbors, v.unique_parents
	#print a.eulerianess()
	# to test eulerian walk output
	#superstr = a.eulerian_walk()
	#print len(superstr)
	#print superstr
	#superstring = superstr[0] + ''.join(map(lambda x: x[-1], superstr[1:]))
	#print "Eulerian walk results: ", superstring
	a.collapse_driver()
	# for x in a.G.keys():
	# 	cur, nbor_list = a.G[x].name, [x.name for x in a.G[x].neighbors]
	# 	print cur, nbor_list,
	# 	if nbor_list:
	# 		print nbor_list.count(nbor_list[0]) == len(nbor_list) 
	# 		a.merge_nodes(a.G[cur],a.G[nbor_list[0]])
	# 	else: print ""
		#print a.G[x+1]

	# Get all leaves in tree, collapse them!
	#print len([x for x in dbg.iterkeys()])
	#print "LEAVES: ", a.find_leaves()
	#for leaf in a.find_leaves():
		#print leaf
	#	a.concat(leaf)
	#print len([x for x in dbg.iterkeys()])

	# find all roots in tree - tested this and works in all edge cases!!!!
	#roots = [x for x in dbg.iterkeys() if dbg[x].is_head()]

	#print "number of roots found: ", len(roots)
	#print max(len(a.bfs(dbg,root)) for root in roots)

	# find longest depth-first path in graph

	# max_root = None
	# max_len = 0
	# for root in roots:
	# 	z = len(a.dfs(a.G,root))
	# 	if z > max_len:
	# 		max_root = root
	# 		max_len = z

	# print "root with longest path: ", max_root, " has length: ", max_len
	# dfsl = dfs(a.G, max_root)
	# for i in range(len(dfsl)-1):
	# 	if a.G[dfsl[i]].is_collapsible() and a.G[dfsl[i+1]].is_collapsible() and a.G[dfsl[i+1]] in a.G[dfsl[i]].neighbors:
	# 		merge_nodes(a.G[dfsl[i]],a.G[dfsl[i+1]])

	#dfs = [x.name for x in a.dfs(dbg,max_root)]
	#print dfs
	#bfs = [x.name for x in a.bfs(dbg,max_root)]
	#print "depth first search results: ", dfs
	#print "breadth first search results: ", bfs

	#print "roots found through bfs: ", [x for x in bfs if a.G[x].is_head()]
	#print "leaves found through bfs: ", [x for x in bfs if a.G[x].is_leaf()]
	#print "roots found through bfs: ", [x for x in dfs if a.G[x].is_head()]
	#print "leaves found through dfs: ", [x for x in dfs if a.G[x].is_leaf()]
	
	# count = 0
	# for x in dbg.iterkeys():
	# 	if dbg[x].is_head():
	# 		count +=1
	# print count
	#maxlen = max([len(x) for x in a.G.iterkeys()])
	#contigs = [x for x in a.G.iterkeys() if len(x) == maxlen]

	dot = a.weighted_dot_file_generator(a.G)
	#print dot

	with open("out.dot", "w") as f:
		f.write(dot)
		f.close()

