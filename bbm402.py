import time
import networkx as nx
import random
import matplotlib
import array
import math
import sys
import collections
from random import choices
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from networkx.drawing.nx_agraph import graphviz_layout

edge_pE = {}# A dictionary to map edges with tauE values and pE values
S = []
S2 = []
count = array.array('i', (0 for i in range(0, 6)))
A2 = [0, 1, 2, 4, 6, 12]
A3 = [0, 0, 1, 0, 4, 12]
A4 = [0, 0, 0, 1, 1, 3]
A5 = [0, 0, 0, 0, 1, 6]
A6 = [0, 0, 0, 0, 0, 1]
N1 = 0
Cb = array.array('f', (0 for i in range(0, 6)))
# computes N1
def findN1(G):
    global N1
    vertices = list(G.nodes)
    for v in vertices:
        if G.degree(v) > 3:
            N1 += combination(G.degree(v), 3)

def randomNeighbour(G, this, that, dict_of_graph):
    randNeighbour = random.choices([i for i in list(dict_of_graph[this]) if i not in [that]])[0]
    return randNeighbour

# computes W and tauE
def findW(adamic_adar_scores, k, G):
    W = 0.0
    global edge_pE

    # Calculate tauE, lambdaE and W, L:
    for u, v, score in adamic_adar_scores:
        if (score > k and score <= 1):
            tauE = (G.degree[u]-1) * (G.degree[v] - 1)
            W += tauE
            edge_pE[(u, v)] = tauE

    # Calculate pE:
    for edge in edge_pE:
        edge_pE[edge] = edge_pE[edge] / W
    return W



# algorithm 1
def sample(G, edge_pE, dict_of_graph):

    randomEdge = random.choices(list(edge_pE.keys()), weights=list(edge_pE.values()))[0]
    # Select random uniform neighbour u' such that u' != v
    uPrime = randomNeighbour(G, randomEdge[0], randomEdge[1], dict_of_graph)
    # Select random uniform neighbour v' such that v' != u
    vPrime = randomNeighbour(G, randomEdge[1], randomEdge[0], dict_of_graph)
    return ((randomEdge[0], randomEdge[1]), (randomEdge[0], uPrime), (randomEdge[1], vPrime))

# check if subgraph is triangle
def triangle(nodes):
    # if u' = v' it is a triangle
    if nodes[2] == nodes[3]:
        return True

def determineSubgraphMotif(nodes, G):
    if not triangle(nodes):
        u = nodes[0]
        v = nodes[1]
        uPrime = nodes[2]
        vPrime = nodes[3]
        # 0 -> u, v; 1 -> u, u'; 2 -> v, v'; 3 -> u, v'; 4 -> v, u'; 5 -> u', v'
        exists = {0: True, 1: True, 2: True, 3: False, 4: False, 5: False}
        if G.has_edge(u, vPrime):
            exists[3] = True
        if G.has_edge(v, uPrime):
            exists[4] = True
        if G.has_edge(uPrime, vPrime):
            exists[5] = True

        # Determine sub-graph:
        if exists[3] == False and exists[4] == False and exists[5] == False:
            return "3-path"
        elif exists[5] == True and exists[4] == False and exists[3] == False:
            return "4-cycle"
        elif exists[5] == True and exists[4] == True and exists[3] == True:
            return "4-clique"
        elif exists[5] == True and ((exists[4] == True and exists[3] == False) or (exists[4] == False and exists[3] == True)):
            return "chordal-4-cycle"
        elif exists[5] == False and ((exists[4] == True and exists[3] == False) or (exists[4] == False and exists[3] == True)):
            return "tailed-triangle"
    else:
        return "triangle"




# algorithm 2
def three_path_sampler(G, samples_k, W, dict_of_graph):
    N1 = 0
    # A2,i
    global A2
    global A3
    global A4
    global A5
    global A6
    # Array of subgraphs from sample algorithm
    global S
    # Initialize count_i = 0 for i in [2,6]
    global count
    # Cb_i
    global Cb
    # Run sample samples_k times
    for i in range(samples_k):
        sub = sample(G, edge_pE, dict_of_graph)
        node_list = [sub[0][0], sub[0][1], sub[1][1], sub[2][1]]
        S.append(node_list)
    # For l in [1,samples_k]
    for l in range(samples_k):
        if determineSubgraphMotif(S[l], G) == "4-clique":
            count[1] += A2[5]
            count[2] += A3[5]
            count[3] += A4[5]
            count[4] += A5[5]
            count[5] += A6[5]


        if determineSubgraphMotif(S[l], G) == "chordal-4-cycle":
            count[1] += A2[4]
            count[2] += A3[4]
            count[3] += A4[4]
            count[4] += A5[4]
            count[5] += A6[4]


        if determineSubgraphMotif(S[l], G) == "4-cycle":
            count[1] += A2[3]
            count[2] += A3[3]
            count[3] += A4[3]
            count[4] += A5[3]
            count[5] += A6[3]
        if determineSubgraphMotif(S[l], G) == "tailed-triangle":

            count[1] += A2[2]
            count[2] += A3[2]
            count[3] += A4[2]
            count[4] += A5[2]
            count[5] += A6[2]
        if determineSubgraphMotif(S[l], G) == "3-path":

            count[1] += A2[1]
            count[2] += A3[1]
            count[3] += A4[1]
            count[4] += A5[1]
            count[5] += A6[1]
    # For each i in [2,6]
    for i in range(1,6):
        Cb[i] = (count[i]/samples_k) * (W/A2[i])
    vertices = list(G.nodes)
    for v in vertices:
        if G.degree(v) > 3:
            N1 += combination(G.degree(v),3)
    Cb[0] = N1 - Cb[2] - 2*Cb[4] - 4*Cb[5]

    return Cb

# Calculates combination (n chooses r)
def combination(n ,r):
    return math.factorial(n) / (math.factorial(r)* (math.factorial(n-r)))


def L_xy(x, y, G):
    count = 0
    for neighbour in G.neighbors(x):
        if(neighbour != y):
            if (G.degree[neighbour] >= G.degree[y]):
                count += 1
            elif (int(neighbour) > int(y)):
                count += 1
    return count






# Computes lambda and lambdaE
def findLambda(G):
    L = 0.0
    for edge in edge_pE:
        lambdaE = L_xy(edge[0], edge[1], G) * L_xy(edge[1], edge[0], G)
        L += lambdaE
        edge_pE[edge] = lambdaE

    # Calculate pE:
    for edge in edge_pE:
        edge_pE[edge] = edge_pE[edge] / L
    return L

# Algorithm 3:
def sample_centered( dict_of_graph,edge_pE):
    randomEdge = random.choices(list(edge_pE.keys()), weights=list(edge_pE.values()))[0]
    uPrime = random.choices([i for i in list(dict_of_graph[randomEdge[0]]) if
                             i not in [randomEdge[1] and len(dict_of_graph[i]) > len(dict_of_graph[randomEdge[0]])]])[0]

    vPrime = random.choices([i for i in list(dict_of_graph[randomEdge[1]]) if
                             i not in [randomEdge[0] and len(dict_of_graph[i]) > len(dict_of_graph[randomEdge[1]])]])[0]
    return ((randomEdge[0], randomEdge[1]), (randomEdge[0], uPrime), (randomEdge[1], vPrime))


def isGreater(x, y, G):
    if(int(x) > int(y) or G.degree[x] > G.degree[y]):
        return True
    return False

def is_centered_3_Path(u, v, t, w, G):
    if (isGreater(t, v, G) and isGreater(w, u, G)):
        if(G.has_edge(t, w)):
            return True
    return False


# Algorithm 4:
def centered_sample (G, samples_k, L, dict_of_graph):

    # Array of subgraphs from sample centered algorithm
    global S2
    # Initialize count_i = 0 for i in [4,6]

    count[3] = 0
    count[4] = 0
    count[5] = 0
    # Cb_i
    global Cb
    B = [0, 0, 0, 1, 1, 3]
    # Run sample samples_k times
    for i in range(samples_k):
        sub = sample_centered(dict_of_graph,edge_pE)
        S2.append(sub)

    # For l in [1,samples_k]
    for l in range(samples_k):
        if(is_centered_3_Path(S2[l][0][0], S2[l][0][1], S2[l][1][1], S2[l][2][1], G)):
            if determineSubgraphMotif(S[l], G) == "4-clique":
                count[1] += A2[5]
                count[2] += A3[5]
                count[3] += A4[5]
                count[4] += A5[5]
                count[5] += A6[5]
            if determineSubgraphMotif(S[l], G) == "chordal-4-cycle":
                count[1] += A2[4]
                count[2] += A3[4]
                count[3] += A4[4]
                count[4] += A5[4]
                count[5] += A6[4]
            if determineSubgraphMotif(S[l], G) == "4-cycle":
                count[1] += A2[3]
                count[2] += A3[3]
                count[3] += A4[3]
                count[4] += A5[3]
                count[5] += A6[3]
            if determineSubgraphMotif(S[l], G) == "tailed-triangle":
                count[1] += A2[2]
                count[2] += A3[2]
                count[3] += A4[2]
                count[4] += A5[2]
                count[5] += A6[2]
            if determineSubgraphMotif(S[l], G) == "3-path":

                count[1] += A2[1]
                count[2] += A3[1]
                count[3] += A4[1]
                count[4] += A5[1]
                count[5] += A6[1]
    for i in range(3,6):
        Cb[i] = (count[i]/samples_k) * (L/B[i])
    return Cb
print("deneme",flush=True)
f = open("amazon0312.txt","rb")
G = nx.read_edgelist(f)
findN1(G)
f.close()
scores = nx.adamic_adar_index(G, G.edges)
k = 0.3


dict_of_graph = nx.to_dict_of_dicts(G, nodelist=None, edge_data=None)

enum_time = time.time()
W = findW(scores, k, G)
three_path_time =time.time()
three_path_sampler(G, 10000, W, dict_of_graph)

print("3-path time: ")
print(time.time()-three_path_time)
print(Cb)
L = findLambda(G)
centered_sample(G, 10000, L, dict_of_graph)
print(Cb)
print("enum. time ")
print((time.time() - enum_time))
