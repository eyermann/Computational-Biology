# Standalone node class for de Bruijn graph Implemented by Charles Eyermann
# and Sam Hinh for Computational Biology (CS362) Winter 2016

class Node:
	def __init__(self,name):
		self.name = name
		self.neighbors = []
		self.unique_neighbors = 0
		self.parents = []
		self.unique_parents = 0
		self.indegree = 0
		self.outdegree = 0


	def __eq__(self, other):
		this = self.name
		that = other.name
		if this == that:
			return True
		else:
			return False

	def __hash__(self):
		return hash(self.name)

	def data_dump(self):
		print "|name: ", self.name,
		print "|unique neighbors: ", self.unique_neighbors,
		print "|unique parents: ", self.unique_parents,
		print "|indegree: ", self.indegree,
		print "|outdegree: ", self.outdegree, "|"

	def get_name(self):
		return self.name

	def get_neighbors(self):
		return self.neighbors

	def get_parents(self):
		return self.parents

	def get_balance(self):
		if self.indegree == self.outdegree:
			return True
		else:
			return False

	def get_semi_balance(self):
		if abs(self.indegree - self.outdegree) == 1:
			return True
		else:
			return False

	def get_degree(self):
		return self.indegree - self.outdegree

	def is_head(self):
		if self.get_degree() < 0 and self.indegree == 0:
			return True
		else:
			return False

	def is_branching(self):
		return (self.unique_neighbors*self.unique_parents) > 1

	def is_leaf(self):
		if self.outdegree == 0:
			return True
		else:
			return False

	def is_collapsible(self):
		if self.indegree == 1 and (self.outdegree == 0 or 1):
			return True
		else:
			return False
