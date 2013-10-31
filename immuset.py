'''
Created on 2013-1-4

@author: guannliu@gmail.com
'''
from __future__ import division
from readGraph import Graph, weightedGraph
import sys
import os
import random
import pickle
import numpy as np
import networkx as nx
import util
from simuprop import Propagation

# Bread-first to traverse a node   
def BFSNode(G,v,c):
    source = [v]
    G.node[v]['visited'] = 1
    cc = []
    cc.append(v)
    while len(source) != 0:
        s = source[0]
        G.node[s]['scc'] = c
        for nbr in G.neighbors(s):
            if G.node[nbr]['visited'] == 0:
                G.node[nbr]['visited'] = 1
                source.append(nbr)
                cc.append(nbr)
        source.remove(s)
    return cc

# Traverse the graph using breadth-first-search
def BFSGraph(G):
    graph_list = []
    cc = []
    for v in G.nodes():
        G.node[v]['visited'] = 0
    for v in G.nodes():
        c = 0
        if G.node[v]['visited'] == 1:
            continue
        else:
            cc.append(BFSNode(G,v,c))
            c += 1
    cc.sort(key=len,reverse=True)
    for cci in cc:
        graph_list.append(G.subgraph(cci).copy())
    return graph_list

'''
The exponential function to calculate frequency related probability
for each edge
'''
def weightProb(freq):
    return 1 - 0.99*np.exp(-0.01*freq)

# Bond percoalation on the network with infection probability p
def infSample(p,G):
    SG = G.copy()
    for u,v in SG.edges():
        if graph_name == 'case':
			p = weightProb(G[u][v]['weight'])
        if p <= random.random():
            SG.remove_edge(u,v)
    return SG

def expectedInfection(imSet,pG,sability):
    gn = graph_name
    nodesnum = G.number_of_nodes()
    pIG = pG.copy() 
    for i in imSet:
        if i in pIG:
            pIG.remove_node(i)
    cc = BFSGraph(pIG)
    x = 0  #counter
    print 'max,',[len(cci) for cci in cc][0:100]
    for i in range(0,1):
        for v in cc[i].nodes():
            reducedProb = {}
            for node in pIG.nodes():
                reducedProb[node] = 0
            x += 1
            if v in imSet:
                continue
            else:
                #temporarily add v to the immunization set.
                iPG = cc[i].copy()
                iPG.remove_node(v)
                ccr = BFSGraph(iPG)
                for ecc in ccr:
                    for savedn in ecc.nodes():
                        reducedProb[savedn] = (len(cc[i]) - len(ecc)) * \
                            util.lossFunc(G, savedn, gn) / nodesnum
            sa = sum(reducedProb.values()) / util.immcostFunc(G, v, gn)
            sability[v] += sa

# calculate the incremental savability for immunizing each candidate
def incsavability(node,p,imSet,R):
    r = 1
    incsa = 0
    gn = graph_name
    while r <= R:
        reducedProb = {}
        for v in G.nodes():
            reducedProb[v] = 0
        sampleNets = infSample(p, G)
        #immunize the candidate node to derive the incremental benefit
        imNets = sampleNets.copy()
        imNets.remove_nodes_from(imSet)
        for v in imNets.nodes():
            imNets.node[v]['visited'] = 0
        cc = imNets.subgraph(BFSNode(imNets, node, 0)).copy()
        # get the component node lies in
        theCpn = cc.copy()
        theCpn.remove_node(node)
        ncc = BFSGraph(theCpn)
        for ecc in ncc:
            for savedn in ecc.nodes():
                reducedProb[savedn] = (len(cc) - len(ecc)) * \
                            util.lossFunc(G, savedn, gn) / len(G)
        sa = sum(reducedProb.values()) / util.immcostFunc(G, node, gn)
        incsa += sa
        r += 1
    return incsa / R

def imNodes(size,p,R,G,graph_name):
    """
    Parameters:
    size: immunize size
    p: probability to infect neighbors
    R: simulation repeated R times for each vertex
    G: the graph
    """
    imSet = []
    sability = {}
    r = 1
    for eachnode in G.nodes():
        sability[eachnode] = 0
	# first round, sort the savability of all the nodes.
    while r <= R:
        sampleNets = infSample(p, G)
        expectedInfection(imSet, sampleNets, sability)
        stb = {}
        for candim in sability:
            stb[candim] = sability[candim] / r
        sortsab = sorted(stb.items(), key=lambda d:d[1], reverse=True)
        r += 1
    for candim in sability:
        stb[candim] = sability[candim] / r
    sortedlist = sorted(stb.items(), key=lambda d:d[1], reverse=True)
    topcand, satop = sortedlist[0]
    imSet.append(topcand)
    del sortedlist[0]
    print 'The 1 immunization target: %d'%topcand
    
    imsize = 2
    tmpflag = 0
    validity = {}
    for k,v in sortedlist:
        validity[k] = False
    while True:
        candtar, sa = sortedlist[0]
        nextcand, nsa = sortedlist[1]
        if validity[candtar]:
            imSet.append(candtar)
            del sortedlist[0]
            # print "Valid!!The %dth immunization target: %s"%(imsize, candtar)
            imsize += 1
            for k,v in sortedlist:
                validity[k] = False
            tmpflag = 0
        else:
            inc = incsavability(candtar, p, imSet, R)
            # print 'Incremental Benefit of %s is %.3f.......'%(candtar,inc)
            validity[candtar] = True
            if inc >= nsa:
                imSet.append(candtar)
                del sortedlist[0]
                print "The %d immunization target: %s"%(imsize, candtar)
                imsize += 1
                for k,v in sortedlist:
                    validity[k] = False
                tmpflag = 0
            else:
                # put the candtar to the right order of incremental benefit
                del sortedlist[0]
                for i in range(1, len(sortedlist)):
                    ithcand, ithsa = sortedlist[i]
                    if ithsa <= inc:
                        sortedlist.insert(i, (candtar, inc))
                        break
        tmpflag += 1
        if imsize > size:
            break        

    # record the obtained immunization set.
    PATH = "immunization_set"
    if not os.path.exists(PATH):
        os.mkdir(PATH)
    imm_file = open("%s/%s_%d_immset.txt"%(PATH,graph_name,size),"w")
    for t in imSet:
	    imm_file.write("%d\n"%t)
    imm_file.close()
	
    return imSet
if __name__=="__main__":
	if len(sys.argv) <= 3:
		print 'USAGE: python immuset.py [graph_name] [immu_size] [diffusion_p]'
		sys.exit(1)
	graph_name = sys.argv[1]
	graph_path = "data/%s.txt"%graph_name

	if graph_name == 'gdmb':
		G = Graph(graph_path,"\t", False).G
		# load each blogger's article number, treat as immunization cost
		af = file("data/blogger-article.txt","r").readlines()
		b = map(lambda x:int(x.split("\t")[0]), af)
		a = map(lambda x:int(x.strip().split("\t")[1]), af)
		ba = dict(zip(b,a))
		# put the number of articles on the node
		for u in G.nodes():
			G.node[u]['num'] = ba[u]
	elif graph_name == 'case':
		# load the case graph with weight on each edge
		G = weightedGraph(graph_path,"\t").readWeightedGraph()
		# calculate the communication freq. for each node, then the loss is a function
		# with the frequency
		for n in G.nodes():
			G.node[n]['comfreq'] = 0
			for nn in G.neighbors(n):
				G.node[n]['comfreq'] += G[n][nn]['weight']
	else:
		G = Graph(graph_path,"\t", False).G
	
	print 'Graph loaded successfully! Nodes: %d, Edges: %d'%(G.number_of_nodes(), G.number_of_edges())
	# print sum([G.degree(n) for n in G.nodes()])/G.number_of_nodes()
	# print sum(nx.clustering(G).values())/G.number_of_nodes()
	imm_size = int(sys.argv[2])
	diffu_prop = float(sys.argv[3])
	
	imset = imNodes(imm_size,diffu_prop,5,G,graph_name)
	print 'Immunizaiton targets generated! Find them in immunization_set/%s_%d_immset.txt'%(graph_name,imm_size)
	print '============Now infection simulations with the immunization targets==========='
	# import propagation and test the diffusion
	pro = Propagation(diffu_prop, imm_size, imset, G, graph_name)
	pro.initialInfection()
	infsize, infloss = pro.cascade()
	print 'Total Infection Sizes: %.3f'%infsize
	print 'Total Infection Loss: %.3f'%infloss
	# record the results
	wf = open("experiment_results.txt","a")
	wf.write("The immunization results for %s with %d CEIT immunization targets:\n"%(graph_name, imm_size))
	wf.write("=====================Total Infection Sizes: %.3f====================\n"%infsize)
	wf.write("===================Total Infection Loss: %.3f=====================\n\n\n"%infloss)
	wf.close()

