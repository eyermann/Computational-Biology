# Sequence Assembler for CS362 Project 3
# Authors: Charles Eyermann and Sam Hinh

from simulator import *
from node import * 
import argparse, copy
from collections import Counter

class Assembler:

	def __init__(self, data, k = None):
		'''
		Pseudo-overloaded constructor to allow graphs to be built from scratch
		with reads files, but also allows partial graph reconstruction from dot
		files. This feature is a work in progress.
		'''

		if k is None:
			self.reads = []
			self.G = {}
			self.subgraphs = []
			self.contigs = []
			self.keys = []
			self.header = ""
			self.G = self.dot_file_to_graph(data)
			self.balanced_nodes = 0
			self.semi_balanced_nodes = 0
			self.unbalanced_nodes = 0
			self.branching_nodes = 0

		else:
			self.reads = []
			self.G = {}
			self.subgraphs = []
			self.contigs = []
			self.header = ""
			self.k = k
			self.balanced_nodes = 0
			self.semi_balanced_nodes = 0
			self.unbalanced_nodes = 0
			self.branching_nodes = 0

			for row in clargs.read_file:
				if row[0] == ">":
					self.header += row.strip()
				else:
					self.reads.append(row.strip())


			if clargs.verbose:
				print "Building de Bruijn graph with", str(clargs.kmer_length) + "-mers... please wait!"

			self.build_DBG()

			if clargs.verbose:
				print "Finished building de Bruijn graph!\n"


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

		self.populate_parents()
		self.keys = self.G.keys()
		for item in self.G.iterkeys():
			self.G[item].unique_neighbors = len(set([x for x in self.G[item].neighbors]))
			self.G[item].unique_parents = len(set([x for x in self.G[item].parents]))
			if self.G[item].is_branching():
				self.branching_nodes += 1
			if self.G[item].get_balance():
				self.balanced_nodes += 1
			elif self.G[item].get_semi_balance():
				self.semi_balanced_nodes += 1
			else:
				self.unbalanced_nodes += 1
		
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

	def prettyprint(self):
		print "This is the graph!"
		print "| Key  || Neighbors"
		print "_"*78
		for k,v in self.G.iteritems():
			print "|", k, "||", [n for n in v.neighbors]

	def eulerian_walk(self,start):
		walk = []
		gc = copy.deepcopy(a.G)
		#choose random node
		#start = random.choice(self.G.keys())
		#start = self.G.iterkeys().next()
		#print "Starting node: ", start
		def next(node):
			#print [x.name for x in self.G[node].get_neighbors()]
			#print len(self.G[node].get_neighbors())
			neighbors = gc[node].get_neighbors()
			while len(neighbors) > 0:
				#print neighbors
				nxt = neighbors.pop()
				next(nxt)
			walk.append(node)
		next(start)
		walk = walk[::-1]
		return walk

	def find_leaves(self):
		return set([x for x in self.G.iterkeys() if self.G[x].is_leaf()])

	def get_eulerianess(self):
		return ((self.unbalanced_nodes == 0 and self.semi_balanced_nodes == 0) or 
				(self.unbalanced_nodes == 0 and self.semi_balanced_nodes == 2))

	def get_cycle(self, dbg):
		path = []
		pathset = set()
		visited = set()
		
		def visit(node):
			if node in visited:
				return None
			visited.add(node)
			path.append(node)
			pathset.add(node)
			neighbors = self.G[node].neighbors
			for neighbor in neighbors:
				try: 
					if neighbor in pathset or visit(neighbor):
							return path
				except RuntimeError as e:
					"recursion depth exceeded"
					pass
			pathset.remove(node)
			path.remove(node)

		for key in self.G.iterkeys():
			paths = visit(self.G[key].name)
		return pathset, path

	def bfs(self, dbg, start):		
		visited, q = set(), [start]
		while q:
			v = q.pop(0)
			if self.G[v] not in visited:
				visited.add(self.G[v])
				p = [item for item in self.G[v].neighbors if item not in visited]
				q.extend(p)
		return visited

	def is_connected(self):
		return len([x for x in a.G.iterkeys() if a.G[x].is_head()]) == 1

	def dfs(self, dbg, start):
		visited, stack = set(), [start]
		while stack:
			v = stack.pop()
			if self.G[v] not in visited:
				visited.add(self.G[v])
				q = [item for item in self.G[v].neighbors if item not in visited]
				stack.extend(q)
		return visited

	def populate_parents_dfs(self, dbg, start):
		visited, stack = [], [start]
		while stack:
			v = stack.pop()
			if v not in visited:
				visited.append(v)
				q = [item for item in dbg[v].neighbors if item not in visited]
				for n in q:
					self.G[n].parents.append(self.G[v].name)
				stack.extend(q)
		return visited

	def populate_parents(self):
		for key in self.G:
			for neighbor in self.G[key].neighbors:
				self.G[neighbor].parents.append(key)

	def merge_nodes(self,n1,n2):
		#assume we are collapsing upwards i.e ATTGC + TTGCG = ATTGCG
		newname = n1.name[0] + n2.name
		newnode = Node(newname)
		newnode.parents = n1.parents
		newnode.neighbors = n2.neighbors
		newnode.indegree = n1.indegree
		newnode.outdegree = n2.outdegree
		self.G[newname] = newnode
		del self.G[n1.name]
		del self.G[n2.name]

	def concat(self,node):
		contigs = []
		if self.G[node.name].get_parents():
			parents = self.G[node.name].get_parents()
			try:
				if node.is_collapsible() and self.G[parents[0]].is_collapsible():
					if len(node.get_parents()) == 1:
						parent = self.G[node.get_parents()[0]]
						cur = node
						kmer_len = len(parent.name)

						if parent.get_parents():
							real_parent = parent.get_parents()[0]
							if parent.name[1:] == cur.name[:(kmer_len-1)]:
								conc = parent.name[0] + cur.name
								newnode = self.G[conc] = Node(conc)
								self.G[real_parent].neighbors.append(newnode.name)
								newnode.parents.append(real_parent)
								newnode.indegree = len(newnode.parents)
								newnode.unique_parents = len(set(self.G[real_parent].parents))
								newnode.unique_neighbors = 0								
								newnode.name = conc
								self.G[real_parent].neighbors = [x for x in self.G[real_parent].neighbors if x != parent]
								del self.G[parent.name]
								del self.G[node.name]
								self.concat(newnode)

				elif node.is_collapsible():
					parent = node.get_parents()[0]
					kmer_len = len(parent)
					cur = node

					if parent[1:] == cur.name[:(kmer_len-1)]:
						conc = parent[0] + cur.name
						newnode = self.G[conc] = Node(conc)	
						newnode.name = conc
			except KeyError,e:
				print "keyerror"
				pass

			
	def collapse_driver(self):
			if self.keys:
				for k,v in self.G.items():
					#print "current list of values in dict: ", [x for x in self.G.keys()]
					#print "current list of nodes in dict: ", self.keys
					cur, nbor_list = k, [x for x in v.neighbors]
					#print "node: '", cur, "' neighbors: ", nbor_list
					if nbor_list:
						try:

							# if all nodes in the neighbors list for the current node are the same and they overlap, 
							# we can collapse the two nodes 

							#print self.G[nbor_list[0]].unique_neighbors
							#print len(self.G[cur].get_parents())

							next_node = nbor_list[0]
							if (nbor_list.count(next_node) == len(nbor_list)
							and (cur[1:] == next_node[:(len(cur)-1)])
							and (self.G[cur].unique_parents <= 1)
							and (self.G[next_node].unique_parents <= 1)
							and (self.G[cur].unique_neighbors <= 1)
							and (self.G[next_node].unique_neighbors <= 1)):

								#print "Homogenous neighbor list AND",
								#print "Should be identical (overlap): <" , cur[1:], "==", nbor_list[0][:(len(cur)-1)],">", 
								#print ", so merge these nodes: ", cur, "and ", nbor_list[0]

								n1,n2 = self.G[cur], self.G[next_node]
								newname = n1.name[0] + n2.name

								newnode = Node(newname)
								newnode.parents = n1.parents
								newnode.neighbors = n2.neighbors
								newnode.unique_parents = len(set(n1.parents))
								newnode.unique_neighbors = len(set(n2.neighbors))

								for parent in a.G[n1.name].parents:
									a.G[parent].neighbors = [newnode.name if x == n1.name else x for x in a.G[parent].neighbors]

								print "DATA DUMP: ",
								newnode.data_dump()
							
								self.G[newname] = newnode
								self.keys.append(newname)
														
								#if n1.name in self.keys:
								self.keys.remove(n1.name)
								#if n2.name in self.keys:
								self.keys.remove(n2.name)
								del self.G[n1.name]
								del self.G[n2.name]
								#print len(self.keys)
								self.keys = self.G.keys()
								self.collapse_driver()


							else:
								#print "If False: ", nbor_list.count(nbor_list[0]) == len(nbor_list), "| AND",
								#print "NOT identical: <" , cur[1:], nbor_list[0][:(len(cur)-1)],">",
								#print "we won't merge these nodes: ", cur, "and ", nbor_list[0]						
								continue

						except KeyError,e:
							print "ABORT ABORT: ", str(e)						
							break	
					else: 
						pass	

			return self.G

	def trim_cycles(self):
		if clargs.verbose:
			print "Graph size pre-trim: ", len(self.G)
		if self.get_cycle(self.G)[0]:
			while self.get_cycle(self.G)[0]:
				lst_set =  self.get_cycle(self.G)[0]
				lst1 = [x for x in self.get_cycle(self.G)[1]]
				for item in lst_set:
					del self.G[item]
				out = lst1[0]
				for i in range(1, len(lst1)):
					out += lst1[i]
		if clargs.verbose:
			print "Graph size post-trim: ", len(self.G)


	def get_contigs(self):
		roots = [x for x in self.G.iterkeys() if self.G[x].is_head()]
		passed_over = 0
		if clargs.verbose:
			print "Generating contigs: ", len(roots), "possible"
			counter = len(roots)
		for root in roots:
			try:
				superstr = self.eulerian_walk(root)
				results = superstr[0]
				for i in range(1,len(superstr)):
					if superstr[i] != superstr[i-1]:
						results += superstr[i][-1]
				self.contigs.append(results)
				counter -= 1

				if clargs.verbose:
					#print "CONTIG: ", results
					#print "\n"
					print counter
			except RuntimeError as e:
				if clargs.verbose:
					#print "cycle found - passing over this contig!\n"
					passed_over += 1
					counter -= 1
					print counter
				continue
		if clargs.verbose:
			print "Finished generating contigs"

		return self.contigs, passed_over

	def flush_contigs(self,contigs):
		out = ""
		with open("our_contigs.txt", "w") as f:
			header = "> Contigs generated: "+str(len(self.contigs))+" | Contigs skipped due to cycles: "+str(contigs[1])+" | Source Header: '"+self.header[1:]+"'"
			out += header+"\n"
			for contig in self.contigs:
				out += contig[::-1]+"\n"
			f.write(out)
		if clargs.verbose:
			print "Wrote contigs file to 'our_contigs.txt'"

	def trimmer(self,node):
		length_threshold = 2*self.k #cutoff for branch being considered junk
		counter = 0
		start = self.G[node]
		cur = self.G[node]
		path = []

		while cur.unique_neighbors != 0:
			counter += 1
			path.append(cur)
			# print len(cur.neighbors)
			#b = min(set(cur.neighbors) key=cur.neighbors.count)
			least_common = Counter(cur.neighbors).most_common()[-1][0]
			cur = a.G[least_common]

		for node in path:
			del a.G[node.name]



	def trim_branches(self):
		length_threshold = 2*self.k #cutoff for branch being considered junk
		branch_nodes = [x for x in self.G.iterkeys() if a.G[x].is_branching() and a.G[x].unique_neighbors > 1]
		merger_nodes = [x for x in self.G.iterkeys() if a.G[x].is_branching() and a.G[x].unique_parents > 1]
		print len(branch_nodes), len(merger_nodes)
		to_trim = []
		if clargs.verbose:
			counter = len(branch_nodes)
		for branch_root in branch_nodes:
			# if the branch is below the 2k threshold as cited in Pevzner et al. 2001, mark it for trimming
			if clargs.verbose:
				if (counter % 25) == 0:
					print counter, "Paths left"
				counter -= 1
			try:
				self.trimmer(branch_root)
			except KeyError,e:
				#print "already removed this node"
				pass



if __name__ == "__main__":
	# BEGIN ARGPARSE CODE
	p = argparse.ArgumentParser(description='de Bruijn Graph-Driven Sequence Assembler.\r\n Written by Charles Eyermann and Sam Hinh for CS362 Winter 2016')
	p.add_argument('read_file', type=argparse.FileType('r'), help="This should be the file containing your sequence reads")
	p.add_argument('kmer_length', type=int, help="Desired k-mer length (integer value).")
	p.add_argument('-v', '--verbose', help="Print out some more info about the program at runtime", action="store_true")
	clargs = p.parse_args()

	if clargs.kmer_length < 2:
		p.print_help()
		print "\n"
		print "FATAL: kmer length must be 2 or more! - suggested to use kmer length of 4 at minimum."
		sys.exit(1)
	# END ARGPARSE CODE

	# BEGIN MAIN CODE
	if clargs.verbose:
		print "\nVERBOSE MODE ON\n"

	a = Assembler(clargs.read_file, clargs.kmer_length)

	if clargs.verbose:

		print "Initial graph statistics:"
		print "\tUnbalanced nodes: ", a.unbalanced_nodes
		print "\tBalanced nodes: ", a.balanced_nodes
		print "\tSemi-balanced nodes: ", a.semi_balanced_nodes
		print "\tTotal nodes: ", a.unbalanced_nodes + a.balanced_nodes + a.semi_balanced_nodes
		print "\tBranching nodes: ", a.branching_nodes
		print "\tGraph is Eulerian: ", a.get_eulerianess()
		print "\n"

	print len(a.G)
	a.trim_branches()
	print len(a.G)
	#a.trim_cycles()
	#c = a.get_contigs()
	#a.flush_contigs(c)
	#print [x for x in a.G.keys() if a.G[x].is_branching()]
	

	#if not a.is_connected():
	# roots = [x for x in a.G.iterkeys() if a.G[x].is_head()]
	# leaves = [x for x in a.G.iterkeys() if a.G[x].is_leaf()]
	# print "There exist", len(roots), "roots"
	# print "There exist", len(leaves), "leaves"
	# 	for root in roots:
	# 		sub = a.dfs(a.G,root)
	# 		a.subgraphs.append(sub)
	#print a.get_cycle()
	# print len(a.subgraphs)
	#print a.get_cycle(a.G)[1]

	# to test eulerian walk output
	#superstr = a.eulerian_walk()
	#print len(superstr)
	#print superstr
	#superstring = superstr[0] + ''.join(map(lambda x: x[-1], superstr[1:]))
	#print "Eulerian walk results: ", superstring

	#for k,v in a.G.items():
	#	print k, v.unique_neighbors, v.unique_parents
	#a.collapse_driver()
	# for k,v in a.G.items():
	# 	print k,v.unique_neighbors, v.unique_parents


	# Get all leaves in tree, collapse them!
	#print len([x for x in dbg.iterkeys()])
	#print "LEAVES: ", [x for x in a.find_leaves()]
	# for leaf in a.find_leaves():
	# 	print leaf
	# 	a.concat(a.G[leaf])
	#print len([x for x in dbg.iterkeys()])

	# find all roots in tree - tested this and works in all edge cases!!!!
	#roots = [x for x in a.G.iterkeys() if a.G[x].is_head()]
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

	# rtd = [x for x in roots if x != max_root]
	# for r in rtd:
	# 	cur = a.G[r]
	# 	nxt = a.G[cur.neighbors[0]]
	# 	print cur.name,nxt.name
	# 	while nxt.unique_parents < 2:
	# 		try:
	# 			del a.G[r]
	# 			cur = nxt
	# 		except KeyError,e:
	# 			break


	#dfsl = a.dfs(a.G, max_root)
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

	with open("pretrim.dot", "w") as f:
		f.write(dot)
		f.close()

