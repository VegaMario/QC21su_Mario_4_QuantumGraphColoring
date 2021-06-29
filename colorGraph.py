import networkx as nx
from dwave.system import EmbeddingComposite, DWaveSampler
import neal
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from dimod import BinaryQuadraticModel
import math
import sys
import copy
from dimod.generators.constraints import combinations
from hybrid.reference import KerberosSampler
from pyqubo import Spin
from tabu import TabuSampler
from mpl_toolkits.mplot3d import Axes3D
import greedy
from dwave_qbsolv import QBSolv


# function for generating variables in a particular format "name_color"
def generate_var(name, col):
    return f"{name}_{col}"


# function used for reading file containing the region names
def read_nodes(filename):
    with open(filename, "r") as f:
        content = f.readlines()

    nodes = []  # we will store the nodes in this list

    for line in content:
        new_line = line.rstrip()  # remove the empty space at the end of the name

        if new_line:  # if we have a non-empty line
            nodes.append(new_line)

    return nodes


# function used for reading the file containing the bordering pair of regions
def read_edges(filename):
    with open(filename, "r") as f:
        content = f.readlines()

    edges = []  # we will store the edges in this list

    for line in content:
        new_line = line.rstrip()  # remove the empty space at the end of the name

        if new_line:
            new_line = tuple(map(str, new_line.split(' ')))  # make a tuple from the pair of nodes sharing edges
            edges.append(new_line)  # append the tuple to the edges list

    return edges


# function used for generating the QUBO. Based on my own equation.
def gen_QUBO(nodes, edges, colors, gamma):
    Q = defaultdict(int)  # create the QUBO matrix
    offset = 1 * len(nodes) * gamma  # calculate the offset

    # constraint
    for i in nodes:  # for each node
        for j in range(colors):  # for each color
            var_j = generate_var(i, j)  # generate the variable
            Q[(var_j, var_j)] += -1 * gamma  # fill in the linear diagonal terms of the QUBo matrix
            for k in range(j + 1, colors):  # loop for the off-diagonal terms of the QUBO matrix
                var_k = generate_var(i, k)  # generate variable
                Q[(var_j, var_k)] += 2 * gamma  # fill in the quadratic off-diagonal terms of the QUBo matrix

    # objective
    for i in range(colors):  # for each color
        for j, k in edges:  # for each edge
            offset += 1  # increment the offset
            var_j = generate_var(j, i)  # generate the variables
            var_k = generate_var(k, i)
            Q[(var_j, var_j)] += -1  # fill in the linear diagonal terms
            Q[(var_k, var_k)] += -1
            Q[(var_j, var_k)] += 2  # fill in the quadratic off-diagonal terms

    return Q, offset  # return both the QUBO matrix and the offset


# convert the QUBO matrix into a BQM
def qubo_to_bqm(Q, o):
    bqm = BinaryQuadraticModel.from_qubo(Q, offset=o)
    return bqm


# function for checking the solution of the best answer found
def check_soln(sample, edges):
    is_correct = 1  # 1 for correct, 0 for incorrect
    constraint_error = 0  # 0 if no error, 1 if error

    for node1, node2 in edges:  # for each pair of nodes forming an edge

        col1 = -1  # color 1 and color 2 to keep track of the colors in the pair of nodes forming the edge
        col2 = -1

        print("Comparing " + node1 + " and " + node2 + ":")

        for var, val in sample.items():  # for every item in the sample
            node, col = var.split('_')  # get the node and the color information from the variable name

            if val == 1:  # if the variable had a value of 1
                if node1 == node:  # if the variable node is the same as node1 of the edge
                    print("Prov: {} \nColor: {}\n".format(node, col))
                    if col1 > -1:  # if color 1 has been previously assigned a color
                        print('Incorrect Solution: {} has more than one color assigned'.format(node1))
                        is_correct = 0  # the solution is incorrect
                        constraint_error = 1  # this type of error is a constraint error
                        break

                    col1 = int(col)  # set the node color equal to the variable's color

                if node2 == node:  # if the variable node is the same as node2 of the edge
                    print("Prov: {} \nColor: {}\n".format(node, col))
                    if col2 > -1:  # if the color 1 has been previously assigned a color
                        print('Incorrect Solution: {} has more than one color assigned'.format(node2))
                        is_correct = 0  # the solution is incorrect
                        constraint_error = 1  # this type of error is a constraint error
                        break

                    col2 = int(col)  # set the node color equal to the variable's color

            if col1 > -1 and col2 > -1:  # If both color 1 and color 2 have been assigned a color
                break  # exit from the loop

        if constraint_error:  # if there has been a constraint error, then exit the loop
            break

        if col1 == col2:  # if the colors are the same, then this is an objective function error
            print('Incorrect Solution: There is an occurrence where {} and {} have the same color'.format(node1, node2))
            is_correct = 0
            break

    if is_correct:  # if the solution was correct
        print('The solution is correct')


# Function used for plotting the results
def plot_graph(sample, nodes, edges, colors):
    G = nx.Graph()  # make a graph based on the edges and nodes used
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)

    # credit for this part of the code goes to the original DWave Map Coloring Example
    color_map = {}

    for node in nodes:
        for i in range(colors):
            if sample[node + "_" + str(i)]:
                color_map[node] = i

    node_colors = [color_map.get(node) for node in G.nodes]

    nx.draw_circular(G, with_labels=True, node_color=node_colors, node_size=1000, cmap=plt.cm.rainbow)
    plt.show()


# main function
def main():
    nodes = read_nodes("nodes.txt")  # get the nodes
    edges = read_edges("edges.txt")  # get the edges

    QUBO, offset = gen_QUBO(nodes, edges, 5, 10)  # create the QUBO: gen_QUBO(nodes, edges, colors, gamma)

    sampler = neal.SimulatedAnnealingSampler()  # we are using the neal simulated annealer
    bqm = qubo_to_bqm(QUBO, offset)  # convert our QUBO to a BQM
    sampleset = sampler.sample(bqm, num_reads=1000)  # sample the BQM
    sample = sampleset.first.sample  # store the first/best sample from the sampleset
    print(sampleset)  # print the sampleset
    print(sampleset.first.energy)
    print(sampleset.first.sample)
    check_soln(sample, edges)  # check the solution for errors
    plot_graph(sample, nodes, edges, 5)  # plot it


if __name__ == '__main__':
    main()
