# Sequence Assembler for CS362 Project 3
# Authors: Charles Eyermann and Sam Hinh

from simulator import *
from node import *
import argparse, copy, csv, timeit
from collections import Counter


class Subgraph:
	"""a subgraph class to allow for separating main graph into
	its connected components"""

	def __init__(self, graph):
		self.g = graph
		self.balanced_nodes = 0
		self.semi_balanced_nodes = 0
		self.unbalanced_nodes = 0
		self.branching_nodes = 0
		self.total_nodes = len(self.g.keys())

		for item in self.g.iterkeys():
			self.g[item].unique_neighbors = len(set([x for x in self.g[item].neighbors]))
			self.g[item].unique_parents = len(set([x for x in self.g[item].parents]))
			if self.g[item].is_branching():
				self.branching_nodes += 1
			if self.g[item].get_balance():
				self.balanced_nodes += 1
			elif self.g[item].get_semi_balance():
				self.semi_balanced_nodes += 1
			else:
				self.unbalanced_nodes += 1

		if clargs.verbose:
			self.data_dump()

	def data_dump(self):
		print("|total_nodes: ", self.total_nodes, end=' ')
		print("|balanced_nodes: ", self.balanced_nodes, end=' ')
		print("|semi_balanced_nodes: ", self.semi_balanced_nodes, end=' ')
		print("|unbalanced_nodes: ", self.unbalanced_nodes)
		print("|branching_nodes: ", self.branching_nodes, )
		print("|Eulerianess: ", self.get_eulerianess())

	def get_eulerianess(self):
		"""returns a boolean whether or not the main graph is eulerian or not"""
		return ((self.unbalanced_nodes == 0 and self.semi_balanced_nodes == 0) or
				(self.unbalanced_nodes == 0 and self.semi_balanced_nodes == 2))


class Assembler:

	def __init__(self, data, k=None):
		"""
		Pseudo-overloaded constructor to allow graphs to be built from scratch
		with reads files, but also allows partial graph reconstruction from dot
		files. This feature is a work in progress.
		"""

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
				print("Building de Bruijn graph with", str(clargs.kmer_length) + "-mers... please wait!")

			self.build_DBG()

			if clargs.verbose:
				print("Finished building de Bruijn graph!\n")

	def generate_kmers(self, read, k):
		"""Generate all possible kmers for a given sequence"""
		for nucleotide in range(0, len(read) + 1 - self.k):
			yield read[nucleotide:nucleotide + self.k]

	def build_DBG(self):
		"""Primary function to build de Bruijn graph, which is stored within the
		assembler object"""
		for read in self.reads:
			for kmer in self.generate_kmers(read, self.k):
				# split each kmer into its left and right "k-1 mer"
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

		if clargs.connected:
			self.get_connected_components()

		return self.G

	def dot_file_generator(self, dbg):
		"""build a DOT file from the de Bruijn graph"""
		output = ""
		output += "digraph {\n"

		for key in dbg.iterkeys():
			for v in dbg[key].get_neighbors():
				output += "%s -> %s;\n" % (key, v)
		output += "}"
		return output

	def weighted_dot_file_generator(self, dbg):
		"""build a weighted DOT file from the de Bruijn graph"""
		wd = {}
		output = ""
		output += "digraph {\n"

		for key in dbg.iterkeys():
			for v in dbg[key].get_neighbors():
				k = "%s -> %s" % (key, v)
				if k not in wd:
					wd[k] = 1
				else:
					wd[k] += 1
		for z in wd.iterkeys():
			y = wd[z]
			output += "%s [ label=\"%s\" ];\n" % (z, y)
		output += "}"
		return output

	def weighted_graph_converter(self, dbg):
		"""build a DOT file from the de Bruijn graph"""
		wd = {}
		for key in dbg.iterkeys():
			for v in dbg[key].neighbors:
				k = "%s -> %s" % (key, v)
				if k not in wd:
					wd[k] = 1
				else:
					wd[k] += 1
		return wd

	def dot_file_to_graph(self, dotfile):
		"""allows our assembler to take in a weighted DOT file instead of a reads file"""
		with open(dotfile, "rU") as f:
			for row in f:
				row = row.strip()[:-1]
				if row:
					l = row.split(" -> ")
					lkmer, rkmer = l[0], l[-1]
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

	def prettyprint(self):
		'''prints out the graph to stdout'''
		print("This is the graph!")
		print("| Key  || Neighbors")
		print("_" * 78)
		for k, v in self.G.iteritems():
			print("|", k, "||", [n for n in v.neighbors])

	def eulerian_walk(self, start):
		walk = []
		gc = copy.deepcopy(self.G)

		# choose random node
		# start = random.choice(self.G.keys())
		# start = self.G.iterkeys().next()
		# print "Starting node: ", start
		def next(node):
			# print [x.name for x in self.G[node].get_neighbors()]
			# print len(self.G[node].get_neighbors())
			try:
				neighbors = gc[node].get_neighbors()
				while len(neighbors) > 0:
					# print neighbors
					nxt = neighbors.pop()
					next(nxt)
				walk.append(node)
			except KeyError as e:
				# skip if node already removed
				pass

		next(start)
		# walk = walk[::-1]
		return walk

	def find_leaves(self):
		"""return a set of node names that are considered leaves"""
		return set([x for x in self.G.iterkeys() if self.G[x].is_leaf()])

	def get_eulerianess(self):
		"""returns a boolean whether or not the main graph is eulerian or not"""
		return ((self.unbalanced_nodes == 0 and self.semi_balanced_nodes == 0) or
				(self.unbalanced_nodes == 0 and self.semi_balanced_nodes == 2))

	def get_cycle(self, dbg):
		"""finds and returns the first cycle in our graph"""
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
					pass

			pathset.remove(node)
			path.remove(node)

		for key in self.G.iterkeys():
			paths = visit(self.G[key].name)
		return pathset, path

	def bfs(self, dbg, start):
		"""standard breadth-first search implemented for our graph"""
		visited, q = set(), [start]
		while q:
			v = q.pop(0)
			if self.G[v] not in visited:
				visited.add(self.G[v])
				p = [item for item in self.G[v].neighbors if item not in visited]
				q.extend(p)
		return visited

	def is_connected(self):
		"""returns a boolean describing whether the graph is connected or not"""
		return len([x for x in a.G.iterkeys() if a.G[x].is_head()]) == 1

	def dfs(self, dbg, start):
		"""standard depth-first search implemented for our graph"""
		visited, stack = set(), [start]
		while stack:
			v = stack.pop()
			if self.G[v] not in visited:
				visited.add(self.G[v])
				q = [item for item in self.G[v].neighbors if item not in visited]
				stack.extend(q)
		return visited

	def populate_parents_dfs(self, dbg, start):
		'''a DFS-based traversal of our graph that populates parental information
		for each node.'''
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
		"""faster implementation of our DFS-based algorithm"""
		for key in self.G:
			for neighbor in self.G[key].neighbors:
				self.G[neighbor].parents.append(key)

	def merge_nodes(self, n1, n2):
		"""stand alone node-merging function that assumes we are
		collapsing upwards (i.e ATTGC + TTGCG = ATTGCG)"""
		newname = n1.name[0] + n2.name
		newnode = Node(newname)
		newnode.parents = n1.parents
		newnode.neighbors = n2.neighbors
		newnode.indegree = n1.indegree
		newnode.outdegree = n2.outdegree
		self.G[newname] = newnode
		del self.G[n1.name]
		del self.G[n2.name]

	def concat_driver(self):
		"""driver function for concat"""
		for leaf in self.find_leaves():
			self.concat(self.G[leaf])

	def concat(self, node):
		"""first iteration of our graph simplification algorithm
		performs a bottom up concatenation of leaf nodes"""
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
							if parent.name[1:] == cur.name[:(kmer_len - 1)]:
								conc = parent.name[0] + cur.name
								newnode = self.G[conc] = Node(conc)
								self.G[real_parent].neighbors.append(newnode.name)
								newnode.parents.append(real_parent)
								newnode.indegree = len(newnode.parents)
								newnode.unique_parents = len(set(self.G[real_parent].parents))
								newnode.unique_neighbors = 0
								newnode.name = conc
								self.G[real_parent].neighbors = [x for x in self.G[real_parent].neighbors if
																 x != parent]
								del self.G[parent.name]
								del self.G[node.name]
								self.concat(newnode)

				elif node.is_collapsible():
					parent = node.get_parents()[0]
					kmer_len = len(parent)
					cur = node

					if parent[1:] == cur.name[:(kmer_len - 1)]:
						conc = parent[0] + cur.name
						newnode = self.G[conc] = Node(conc)
						newnode.name = conc
			except KeyError as e:
				if clargs.verbose:
					print("keyerror")
				pass

	def collapse_driver(self):
		"""second iteration of our concat function. Can collapse internal
		nodes within graph"""
		if self.keys:
			for k, v in self.G.items():
				cur, nbor_list = k, [x for x in v.neighbors]
				if nbor_list:
					try:

						# if all nodes in the neighbors list for the current node are the same and they overlap,
						# we can collapse the two nodes

						next_node = nbor_list[0]
						if (nbor_list.count(next_node) == len(nbor_list)
								and (cur[1:] == next_node[:(len(cur) - 1)])
								and (self.G[cur].unique_parents <= 1)
								and (self.G[next_node].unique_parents <= 1)
								and (self.G[cur].unique_neighbors <= 1)
								and (self.G[next_node].unique_neighbors <= 1)):

							# print "Homogenous neighbor list AND",
							# print "Should be identical (overlap): <" , cur[1:], "==", nbor_list[0][:(len(cur)-1)],">",
							# print ", so merge these nodes: ", cur, "and ", nbor_list[0]

							n1, n2 = self.G[cur], self.G[next_node]
							newname = n1.name[0] + n2.name

							newnode = Node(newname)
							newnode.parents = n1.parents
							newnode.neighbors = n2.neighbors
							newnode.unique_parents = len(set(n1.parents))
							newnode.unique_neighbors = len(set(n2.neighbors))

							# update parents!
							for parent in a.G[n1.name].parents:
								a.G[parent].neighbors = [newnode.name if x == n1.name else x for x in
														 a.G[parent].neighbors]

							# if clargs.verbose:
							# 	print "DATA DUMP: ",
							# 	newnode.data_dump()

							self.G[newname] = newnode
							self.keys.append(newname)

							self.keys.remove(n1.name)
							self.keys.remove(n2.name)

							del self.G[n1.name]
							del self.G[n2.name]
							self.keys = self.G.keys()
							self.collapse_driver()


						else:
							# print "If False: ", nbor_list.count(nbor_list[0]) == len(nbor_list), "| AND",
							# print "NOT identical: <" , cur[1:], nbor_list[0][:(len(cur)-1)],">",
							# print "we won't merge these nodes: ", cur, "and ", nbor_list[0]
							continue

					except KeyError as e:
						# if clargs.verbose:
						#	print "ABORT ABORT: ", str(e)
						break
				else:
					pass

		return self.G

	def trim_cycles(self):
		"""a function that can delete cycles within the graph"""

		if clargs.verbose:
			print("Graph size pre-trim: ", len(self.G))
		try:
			if self.get_cycle(self.G)[0]:
				while self.get_cycle(self.G)[0]:
					lst_set = self.get_cycle(self.G)[0]
					lst1 = [x for x in self.get_cycle(self.G)[1]]
					for item in lst_set:
						del self.G[item]
					out = lst1[0]
					for i in range(1, len(lst1)):
						out += lst1[i]

		except KeyError as e:
			pass
		if clargs.verbose:
			print("Graph size post-trim: ", len(self.G))

	def get_contigs(self):
		"""function that finds contigs by finding the heads of each
		disconnected graph in the forest and performing an eulerian walk
		through it"""
		roots = [x for x in self.G.iterkeys() if self.G[x].is_head()]
		passed_over = 0
		if clargs.verbose:
			print("Generating contigs: ", len(roots), "possible")
		counter = len(roots)
		for root in roots:
			try:
				superstr = self.eulerian_walk(root)
				results = superstr[0]
				for i in range(1, len(superstr)):
					if superstr[i] != superstr[i - 1]:
						results += superstr[i][-1]
				counter -= 1
				if len(results) > (2 * self.k):
					self.contigs.append(results)

					if clargs.verbose:
						print("CONTIG: ", results)
						print("\n")
			except RuntimeError as e:
				if clargs.verbose:
					passed_over += 1
					counter -= 1
				continue
		if clargs.verbose:
			print("Finished generating contigs")

		return self.contigs, passed_over

	def flush_contigs(self, contigs):
		"""function to write found contigs to file"""
		out = ""
		with open("contigs.txt", "w") as f:
			# header = "> Contigs generated: "+str(len(self.contigs))+" | Contigs skipped due to cycles: "+str(contigs[1])+" | Source Header: '"+self.header[1:]+"'"
			# out += header+"\n"
			for contig in self.contigs:
				out += contig[::-1] + "\n"
			f.write(out)
		if clargs.verbose:
			print("Wrote contigs file to 'contigs.txt'")

	def trimmer(self, node):
		"""implementation of a trimming heuristic for tips"""
		length_threshold = 2 * self.k  # cutoff for branch being considered junk
		counter = 0
		start = self.G[node]
		cur = self.G[node]
		path = []

		while cur.unique_neighbors != 0:
			counter += 1
			path.append(cur)
			least_common = Counter(cur.neighbors).most_common()[-1][0]
			cur = a.G[least_common]

		for node in path:
			del a.G[node.name]

	def trim_branches(self):
		"""driver function for tip-trimming heuristic similar to the Velvet assembler"""
		length_threshold = 2 * self.k  # cutoff for branch being considered junk
		branch_nodes = [x for x in self.G.iterkeys() if a.G[x].is_branching() and a.G[x].unique_neighbors > 1]

		if clargs.verbose:
			counter = len(branch_nodes)
		for branch_root in branch_nodes:
			if clargs.verbose:
				if (counter % 25) == 0:
					print(counter, "Paths left")
				counter -= 1
			try:
				self.trimmer(branch_root)
			except KeyError as e:
				pass

	def get_connected_components(self):
		gc = copy.deepcopy(self.G)

		head_list = [x for x in gc.iterkeys() if gc[x].is_head()]
		if clargs.verbose:
			print("Potential number of individual connected components: ", len(head_list))
		for head in head_list:
			try:
				path = self.eulerian_walk(head)
			except RuntimeError as e:
				# passing on cycles
				continue

			d = {}
			for n in path:
				try:
					d[n] = gc[n]
					del gc[n]
				except KeyError as e:
					# ignore paths with parts that have been removed already
					pass
			sub = Subgraph(d)
			self.subgraphs.append(sub)

		return self.subgraphs


if __name__ == "__main__":
	# BEGIN ARGPARSE CODE
	p = argparse.ArgumentParser(
		description='de Bruijn Graph-Driven Sequence Assembler.\r\n Written by Charles Eyermann and Sam Hinh for CS362 Winter 2016')
	p.add_argument('read_file', type=argparse.FileType('r'),
				   help="This should be the file containing your sequence reads")
	p.add_argument('kmer_length', type=int, help="Desired k-mer length (integer value).")
	p.add_argument('-v', '--verbose', help="Print out some more info about the program at runtime", action="store_true")
	p.add_argument('-l', '--logging', help="Write run information to log file", action="store_true")
	p.add_argument('-w', '--weighted_dot',
				   help="Output dot file with weighted edges instead of unweighted edges. Cleaner end result.",
				   action="store_true")
	p.add_argument('-c', '--connected',
				   help="First steps for solving eulerian superpath problem. Huge hit to performance with large graphs. Enabling this feature will find all unique eulerian paths and create a list of subgraph objects.",
				   action="store_true")
	clargs = p.parse_args()

	if clargs.kmer_length < 3:
		p.print_help()
		print("\n")
		print("FATAL: kmer length must be 3 or more! - suggested to use kmer length of 4 at minimum.")
		sys.exit(1)

	if clargs.read_file:
		acceptable_start_chars = ["A", "T", "C", "G", ">", "a", "t", "c", "g"]
		if clargs.read_file.next()[0] not in acceptable_start_chars:
			print("\n")
			print("FATAL: unrecognized input file format! - please check your input file.")
			sys.exit(1)
	# END ARGPARSE CODE

	# BEGIN MAIN CODE
	if clargs.verbose:
		print("\nVERBOSE MODE ON\n")

	start = timeit.default_timer()

	a = Assembler(clargs.read_file, clargs.kmer_length)

	stop = timeit.default_timer()
	runtime = stop - start

	if clargs.verbose:
		print("Initial graph statistics:")
		print("\tInput file: ", clargs.read_file.name)
		print("\tUnbalanced nodes: ", a.unbalanced_nodes)
		print("\tBalanced nodes: ", a.balanced_nodes)
		print("\tSemi-balanced nodes: ", a.semi_balanced_nodes)
		print("\tTotal nodes: ", a.unbalanced_nodes + a.balanced_nodes + a.semi_balanced_nodes)
		print("\tBranching nodes: ", a.branching_nodes)
		print("\tGraph is Eulerian: ", a.get_eulerianess())
		print("\tNumber of disconnected subgraphs: ", len(a.subgraphs))
		print("\tGraph generation time: ", runtime)
		print("\n")

	if clargs.logging:
		with open('log.csv', 'ab') as f:
			writer = csv.writer(f)
			row = (clargs.read_file.name, str(clargs.kmer_length),
				   str(a.unbalanced_nodes), str(a.balanced_nodes),
				   str(a.semi_balanced_nodes), str(len(a.G)),
				   str(a.branching_nodes), str(a.get_eulerianess()), runtime)
			writer.writerow(row)

	# print len(a.G)
	# a.collapse_driver()
	# print len(a.G)
	# a.trim_branches()
	# print len(a.G)
	# a.trim_cycles()
	# print len(a.G)
	c = a.get_contigs()
	a.flush_contigs(c)

	if clargs.weighted_dot:
		dot = a.weighted_dot_file_generator(a.G)
	else:
		dot = a.dot_file_generator(a.G)

	with open("out.dot", "w") as f:
		f.write(dot)
		f.close()
