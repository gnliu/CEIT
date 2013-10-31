import networkx as nx
import numpy as np

rf = file("settings.txt").readlines()
rfs = map(lambda x:x.strip(), rf)
alpha = float(rfs[0].split("=")[1])
beta = float(rfs[1].split("=")[1])
alpha_prime = float(rfs[2].split("=")[1])
beta_prime = float(rfs[3].split("=")[1])

# cost functions for each immunization candidate
def immcostFunc(G,node,gn):
	if gn == 'gdmb':
		return 1 + alpha*(1-np.exp(-beta*G.node[node]['num']))
	else:
		return 1 + alpha*(1-np.exp(-beta*G.degree(node)))

# loss function
def lossFunc(G,node,gn):
	if gn == 'case':
		return 1 + alpha_prime*(1-np.exp(-beta_prime*G.node[node]['comfreq']))
	else:
		return 1 + alpha_prime*(1-np.exp(-beta_prime*G.degree(node)))