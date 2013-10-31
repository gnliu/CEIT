***********************************
COST-EFFECTIVE IMMUNIZATION TARGETS
***********************************

Guannan Liu
liugn.10@sem.tsinghua.edu.cn

(C) Copyright 2013, Guannan Liu (liugn.10@sem.tsinghua.edu.cn)

---------------------------------------------------------------------

This is a Python implementation of the proposed algorithm CEIT.
The details of the algorithm is described in the manuscript submitted
to Decision Support Systems (DECSUP-D-13-00082). The code can be used
for review purpose only.

---------------------------------------------------------------------

A. ENVIORONMENT
The code can be run on any platforms with Python 2.7, and open-source
packages NetworkX, Numpy should be pre-installed.

---------------------------------------------------------------------

B. SETTINGS
The parameters settings can be found in settings.txt
alpha and beta represent the parameters for the cost function,
alpha_prime and beta_prime represent the parameters settings for the
loss function.

To replicate the experiment results in the experiment sections, it is
recommended to use the same parameter settings (alpha, beta, 
alpha_prime and beta_prime) in the manuscript.

---------------------------------------------------------------------

C. DATA
The graph data and data concerning cost/loss calculations in the 
experiments are under the directory of 'data/', including:
  a. ccmp.txt: the network structure of CCMP dataset
  b. eumail.txt: the network structure of EUMail dataset
  c. gdmb.txt: the network structure of GDMB dataset
  d. blogger-article.txt: each blogger's number of articles for GDMB
  e. case.txt: the network structure for the dataset used in Case Study 
  Section, and the file is in the format [u v communication frequency]
  
----------------------------------------------------------------------

D. CODE
The packages contain several python scripts.
 |---readGraph.py : load graph file in the program
 |---util.py      : define cost and loss functions and set parameters
 |---immuset.py   : implement the CEIT algorithm and output the results
 |---simuprop.py  : diffusion simulations under independent cascade model
 

The code can be run by:
  python immuset.py [graph_name] [immu_size] [diffusion_p]
    |-- graph_name can be "ccmp", "gdmb", "eumail" and "case"
    |-- immu_size can be set in the range [1,100]
    |-- diffusion_p is recommended to be set in the range (0.02,0.1)
  The immunization set would be generated under the directory immunization_set.
  The results of the diffusion simulations will be recorded in the file 
  experiment_results.txt
  
Note: In implementing CEIT for the case study, the diffusion 
probability is different as a function of each edge's communication 
frequency, thus the program is slightly different, and the probability setting can
be any numbers.
