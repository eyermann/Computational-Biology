# Graph Class for CS362 Project 3
# Authors: Charles Eyermann and Sam Hinh

class Node:
	def __init__(self,name):
		self.name = name
		self.neighbors = []
		self.parents = []
		self.indegree = 0
		self.outdegree = 0

	# def __eq__(self, other):
	# 	this = self.name
	# 	that = other.name
	# 	if this == that:
	# 		return True
	# 	else:
	# 		return False

	def __hash__(self):
		return hash(self.name)

	def get_name(self):
		return self.name

	def get_neighbors(self):
		return self.neighbors

	def get_parents(self):
		return self.parents

	def get_balance(self):
		if self.in_nodes == self.out_nodes:
			return True
		else:
			return False

	def get_semi_balance(self):
		if abs(self.in_nodes - self.out_nodes) == 1:
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

	def is_collapsible(self):
		if self.indegree == 1 and (self.outdegree == 0 or 1):
			return True
		else:
			return False

class Edge:
	def __init__(self, node1, node2):
		self.node1 = node1
		self.node2 = node2

	def __init__(self, node1, node2, weight):
		self.node1 = node1
		self.node2 = node2
		self.weight = weight

	# Assume that edge with direction node1->node2 represented as (node1,node2).
	# To reverse edge, init pairing as (node2,node1).
	def __init__(self, node1, node2, weight, direction):
		self.node1 = node1
		self.node2 = node2
		self.weight = weight

class Graph:

	def __init__(self, graph_data={}):
		self.graph_data = graph_data

	def nodes(self):
		return list(self.graph_data.keys())

	def edges(self):
		return self.get_edges()

	def get_edges(self):
		edge_list = []
		for node in self.graph_data:
			for adj in self.graph_data[node]:
				if {adj, node} not in edge_list:
					edge_list.append({node, adj})
		return edge_list

	def get_node_edges(self,node):
		if node in self.graph_data:
			return self.graph_data[node]

	def add_edge(self, edge):
		(node1,node2) = tuple(edge)
		if node1 in self.graph_data:
			self.graph_data[node1].append(node2)
		else:
			self.graph_data[node1] = [node2]

	def remove_edge(self,edge):
		(node1,node2) = tuple(edge)
		if node1 in self.graph_data:
			self.graph_data[node1].remove(node2)

	def add_node(self, node):
		if node not in self.graph_data:
			self.graph_data[node] = []

	def remove_node(self, node):
		if node in self.graph_data:
			self.graph_data.remove(node)


if __name__ == "__main__":
	pass