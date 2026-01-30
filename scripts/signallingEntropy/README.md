# Signalling entropy estimation with integration of correction method on the interaction networks
The given code computes signalling entropy by integrating a protein interaction network with an expression dataset.

## DoIntegPIN_JMK.R
Includes the functions necessary for integrating a PIN with a gene expression set and compute the interactions probabilities (p.m).

`DoIntegPIN`
Takes the PIN, as adjacency matrix adj.m, and the expression matrix, exp.m, in input, both with genes annotated to ENTREZ ID.
It returns a list with: 
* a, adjacency matrix representing the connected network upon integration of PPI and expression data
* e: corresponding gene expression data, same number of genes as a

`DoIntegPIN.full`
Same as DoIntegPIN, the difference is that the whole network is kept (instead of extracting the fully connected component)
* a: corresponds to the full network, not only the largest connected component
### CompStochMatrix
From a gene expression sample, represented by a vector exp.v, and a network in the form of an adjacency matrix, adj.m it computes and returns the probability matrix (p.m). 

## CompSR_AS.R 
Includes the functions necessary for the entropy calculation (adapted from Teschendorff 2014):
`CompS`
Calculates the local entropy of a node taking as input the stochastic vector p.v.
`CompSN```
Returns the normalised local entropy of a node from the input stochastic vector p.v.
	
`CompMaxSR`
Computes the maximum entropy rate of a network with adjacency matrix adj.m (assumed connected).

`CompSR.ne`
Computes the non-equilibrium entropies for the full networks from the probability matrix p.m. 
It is based on the normalised local entropies.

`CompSR`
Takes a probability matrix, p.m, and the maximum, maxSR, of the adjacency matrix and returns the normalised (by maxSR) global entropies
In this case the non-normalised local entropies are calculated.

## perturbPIN.R
Contains the functions to generate the four types of network perturbations. An adjacency matrix, a percentage value (in the range 0-1) and a seed are required by each function. The modified network, in the form of an adjacency matrix, is returned, together with the interactions that have been changed. 
`reducePIN`: eliminates interactions from the network.
`enlargePIN`: introduces new interactions
`rewirePIN`: keeps one node of the interactions and links it to a different one  
`flipPIN`: removes a certain number of existing interactions and adds new ones


## correctionMethods.R
Comprises the functions to apply the three types of correction we used in the study, including String filtering, topological and semantic.

`correctSTRING`: takes in input a probability matrix (p.m), a threshold thresh and a list with protein interactions and the corresponding combined score (0-1 value) downloaded from STRING  
It returns the adjacency matrix (adj.m) corresponding to the filtered network by the combined score from STRING.

`correctTopology`: takes a probability matrix (p.m), a threshold thresh and the desired method among "jaccard", "dice", "invlogweighted"
It returns the adjacency matrix (adj.m) corresponding to the filtered network by the chosen topological method	

`correctSemantic`: takes a probability matrix (p.m), a threshold thresh and the respective computed semantic neighborhood matrix for the chosen correction method "Resnik", "Lin", "Rel", "Jiang", "Wang"
It returns the adjacency matrix (adj.m) corresponding to the filtered network by the chosen semantic method
