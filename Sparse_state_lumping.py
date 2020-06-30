
from state_lumping_network import StateNetwork
import matplotlib.pyplot as plt
import numpy as np


## Generate State Network from Pathways

pathNetwork = StateNetwork()
pathNetwork.generateStateNetworkFromPaths("output/CPE_Network_Pathways.net",
                                          "output/CPE_states_training_order_2.net",
                                          outputValidationFilename="output/CPE_states_validation_order_2.net",
                                          markovOrder=2, splitWeight=True, validationProb=0.3,
                                          #minPathLength=3, maxPathLength=4, 
                                          seed=3)


pathNetwork = StateNetwork()
pathNetwork.generateStateNetworkFromPaths("output/CPE_Network_Pathways.net",
                                          "output/CPE_order_m3.net",
                                          outputValidationFilename="output/null.net",
                                          markovOrder=3, splitWeight=True, validationProb=0.0,
                                          #minPathLength=3, maxPathLength=4, 
                                          seed=3)

pathNetwork = StateNetwork()
pathNetwork.readFromFile("output/CPE_order_m3.net")

##  Feature matrix

#The feature matrix for a physical node is simply the rows of the state transition 
#matrix for the state nodes belonging to that physical node. To simplify, there is 
#a getFeatureMatrix method that removes all all-zero rows and columns in the feature 
#matrix and provides a mapping back to the original state network. It takes the physical 
#node id as input parameter and returns a tuple (X, T), where X is the feature matrix 
#(np.array) of size (numNonDanglingStateNodes, numLinkedNodes) and T is a dictionary 
#transforming row index in the feature matrix to state node id.

X, rowToStateId = pathNetwork.getFeatureMatrix(101)
print("Feature matrix for the central physical node: \n{}\n rowToStateId: {}".format(X, rowToStateId))



##  Measure pairwise similarity

#Now we can compare rows pairwise and cluster the most similar rows together. 
#The Jensen-Shannon distance is unfortunately not implemented in scikit-learn 
#(though it exist in a pull request), so let's create it.

import numpy as np
from sklearn.metrics import pairwise_distances

def plogp(x):
    x = np.asarray(x)
    return x * np.log2(x, where = x>0)

def entropy(x):
    return -np.sum(plogp(x))

def jensen_shannon_distance(x1, x2):
    x1 = np.asarray(x1)
    x2 = np.asarray(x2)
    mix = (x1 + x2) / 2
    return np.sqrt(entropy(mix) - (entropy(x1) + entropy(x2)) / 2)

def jensen_shannon_distances(X):
    return pairwise_distances(X, metric=jensen_shannon_distance)

print(jensen_shannon_distance(X[0], X[1]))
print(jensen_shannon_distance(X[1], X[2]))




##  Cluster with scikit-learn

#Now we can use general scikit-learn clustering algorithm that takes a custom 
#pairwise distance function as a metric, like Agglomerative Clustering. We can 
#also use for example cosine instead with similar results. Many cluster algorithms 
#take as input the number of clusters you want. For the example feature matrix, it's two.

from sklearn import cluster

model = cluster.AgglomerativeClustering(
    linkage="complete",
    affinity=jensen_shannon_distances,
#     affinity="cosine",
    n_clusters=2
)

labels = model.fit_predict(X)
print("Cluster labels in feature matrix space: {}\nCluster labels in state node space: {}".format(
    labels,
    {rowToStateId[i]:clusterId for i,clusterId in enumerate(labels)}
))



##  Lump entire network
#Now we are ready to run this clustering code on the entire network. For convenience, 
#StateNetwork provides a method clusterStateNodes that takes an argument clusterFeatureMatrix 
#where you can send in a custom clustering function. This function gets a feature 
#matrix as input argument and expects an array of cluster labels back.

def getFeatureClusterFunction(clusterRate=0.25):
    def calcClusters(X):
        numStates, numFeatures = X.shape
        if numStates < 2 or numFeatures < 2:
            # Don't cluster if too small
            return list(range(numStates))

        # Can be an adaptive number of clusters based on entropy reduction
        n_clusters = max(1, int(clusterRate * numStates))
        model = cluster.AgglomerativeClustering(
            linkage="complete",
            affinity=jensen_shannon_distances,
#             affinity="cosine",
            n_clusters=n_clusters
        )

        labels = model.fit_predict(X)
        return labels
    return calcClusters

pathNetwork.clusterStateNodes(clusterFeatureMatrix=getFeatureClusterFunction())



##  Did we lose any information?
#The state network has two methods calcEntropyRate() and calcLumpedEntropyRate() 
#to calculate the average number of bits required to encode the random walk on 
#each physical node.

h1 = pathNetwork.calcEntropyRate()
h2 = pathNetwork.calcLumpedEntropyRate()
print("Entropy rate before: {}, after: {}".format(h1, h2))

from pathlib import Path
pathNetwork.writeLumpedStateNetwork("output/CPE_Network_lumped.net")
print(Path('output/CPE_Network_lumped.net').read_text())





