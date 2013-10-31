#!/usr/bin/python

import networkx as nx
import numpy as np
class Graph:
	def __init__(self, fname, sep, direct=False):
		self.G = nx.Graph()
		self.fname = fname
		self.sep = sep
		self.direct = direct
		if direct == True:	
			self.G = nx.DiGraph()

		graphfile = open(self.fname)
		for line in graphfile:
			el = line.strip("\n")
			if el[0] == "#":
				continue
			edge = el.split(self.sep)
			fe = int(edge[0])
			te = int(edge[1])
			self.G.add_edge(fe , te)
		
class weightedGraph:
	def __init__(self, fname, sep):
		self.G = nx.Graph()
		self.fname = fname
		self.sep = sep
	def readWeightedGraph(self):
		rowid = 0
		graphfile = open(self.fname)
		for line in graphfile:
			rowid += 1
			line = line.strip("\n")
			str = line.split(self.sep)
			u = int(str[0])
			v = int(str[1])
#			if u == v:
#				continue
			w = int(str[2])
			if self.G.has_edge(u, v):
				self.G[u][v]['weight'] += w
			else:
				self.G.add_edge(u, v, weight=w)
		return self.G

			
			
		
		
		
