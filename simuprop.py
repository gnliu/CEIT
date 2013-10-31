'''
Created on 2013-1-4

@author: guannliu@gmail.com
'''
from __future__ import division
import random
from readGraph import Graph, weightedGraph
import networkx as nx
import numpy as np
import util

R = 1000
class Propagation():

    def __init__(self, prop, initsize, immset, G, gn):
		
		self.prop = prop
		self.initsize = initsize
		self.G = G
		self.immset = immset
		self.initset = []
		self.gn = gn
        
    def initialInfection(self):
        while True:
            tempinf = random.choice(self.G.nodes())
            if tempinf not in self.immset and tempinf not in self.initset:
                self.initset.append(tempinf)
            if len(self.initset) == self.initsize:
                break
    def cascade(self):
		r = 0
		tinfsize = 0
		tinfloss = 0
		while r < R:
			nextr = []
			nextr.extend(self.initset)
			evinf = nextr
			infloss = 0
			while nextr:
				thisr = set(nextr)
				nextr = []
				thisinf = []
				for s in thisr:
					for nbr in self.G.neighbors(s):
						if self.gn == 'case':
							ep = 1 - 0.99*np.exp(-0.01*self.G[s][nbr]['weight'])
						else:
							ep = self.prop
						if random.random() < ep and \
							nbr not in self.immset and nbr not in evinf:
							nextr.append(nbr)
							evinf.append(nbr)
							infloss += util.lossFunc(self.G, nbr, self.gn)
			tinfsize += len(evinf)-len(self.initset)
			tinfloss += infloss
			r += 1
		return tinfsize / R, tinfloss / R
 

if __name__ == "__main__":
	
	G = weightedGraph("data/case.txt","\t").readWeightedGraph()
	# calculate the communication freq. for each node, then the loss is a function
	# with the frequency
	for n in G.nodes():
		G.node[n]['comfreq'] = 0
		for nn in G.neighbors(n):
			G.node[n]['comfreq'] += G[n][nn]['weight']
	imset = map(lambda x:int(x.strip()), file("immunization_set/case_5_immset.txt").readlines())
	pro = Propagation(0.05, 10, imset, G, "case")
	pro.initialInfection()
	infsize, infloss = pro.cascade()
	print infsize, infloss
	