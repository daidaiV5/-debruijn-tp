#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics

__author__ = "DAI_EID"
__copyright__ = "Universite Paris"
__credits__ = ["DAI_EID"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "DAI_EID"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
	nucle=['A', 'T', 'C','G']
	with open (fastq_file, 'rt') as handle:
		for i in handle:
			for j in range(0,1):
				if i[j] in nucle:
					seqn= i.strip("\n")
					yield seqn
				
				
def cut_kmer(read, kmer_size):
	for start in range(0,len(read)-(kmer_size-1),1):
		km= read[start:start+kmer_size]    
		yield km
	
def build_kmer_dict(fastq_file, kmer_size):
	kmer_dic={}
	for i in read_fastq(fastq_file) :
		for kmer in cut_kmer(i, kmer_size):
			if kmer in kmer_dic.keys():
				kmer_dic[kmer]+=1
			else :
				kmer_dic[kmer]=1
	return kmer_dic	

def build_graph(kmer_dict):
	G = nx.DiGraph()
	for kmer in kmer_dict:
    		G.add_edge(kmer[:-1],kmer[1:], weight=kmer_dict[kmer]) 
	return G

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    for path in path_list:
        graph.remove_nodes_from(path[(not delete_entry_node):(None if delete_sink_node else -1)])
    return graph


def std(data):
    std=round(statistics.stdev(data),1)
    return std

def path_average_weight(graph, path):
    list_weight=[]
    for i in graph.subgraph(path).edges(data=True):
	    d=i[2]
	    list_weight.append(d['weight'])
    wm =statistics.mean(list_weight)
    return wm

def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
	std_weight=std(weight_avg_list)
	if std_weight == 0:
		std_value=std(path_length)
		if std_value == 0:
			path_long_index=random.randrange(0,len(path_length),1)
		else:
			path_long_index=path_length.index(max(path_length))
	else:
		path_long_index=weight_avg_list.index(max(weight_avg_list))
	del path_list[path_long_index]
	graph1=remove_paths(graph, path_list, delete_entry_node, delete_sink_node)
	return graph1


def solve_bubble(graph, ancestor_node, descendant_node):
    path_list=nx.all_simple_paths(graph,ancestor_node,descendant_node)
    path_list2=[]
    path_length=[]
    weight_avg_list=[]
    for path in path_list:
        path_length.append(len(path))
        weight_avg_list.append(path_average_weight(graph, path))
        path_list2.append(path)
    graph1=select_best_path(graph, path_list2, path_length, weight_avg_list,delete_entry_node=False, delete_sink_node=False)
    return graph1

def simplify_bubbles(graph):
    start= get_starting_nodes(graph)
    end= get_sink_nodes(graph)
    for ancestor in start:
        for descendant in end:
            successors = list(graph.successors(ancestor))
            predecessors = list(graph.predecessors(descendant))
            while (len(successors) < 2 and successors):
                ancestor = successors[0]
                successors = list(graph.successors(ancestor))
            while len(predecessors) < 2 and predecessors:
                descendant = predecessors[0]
                predecessors = list(graph.predecessors(descendant))
            if list(nx.all_simple_paths(graph, ancestor, descendant)):
                graph = solve_bubble(graph,ancestor , descendant)
    return graph

def solve_entry_tips(graph, starting_nodes):
    graph = simplify_bubbles(graph)
    tips = []
    for start in starting_nodes:
        path = [start]
        successors = list(graph.successors(start))
        predecessors = list(graph.predecessors(start))
        while len(successors) < 2 and len(predecessors) < 2 and successors:
            start = successors[0]
            path.append(start)
            successors = list(graph.successors(start))
            predecessors = list(graph.predecessors(start))
        tips.append(path)
    path_lengths=[]
    avg_path_weights=[]
    for path in tips:
        path_lengths.append(len(path))
        avg_path_weights.append(path_average_weight(graph, path))
    graph = select_best_path(graph, tips, path_lengths, avg_path_weights,
                             delete_entry_node=True)
    return graph

def solve_out_tips(graph, ending_nodes):
    graph = simplify_bubbles(graph)
    tips = []
    for end in ending_nodes:
        path = [end]
        successors = list(graph.successors(end))
        predecessors = list(graph.predecessors(end))
        while len(successors) < 2 and len(predecessors) < 2 and predecessors:
            end = predecessors[0]
            path.append(end)
            successors = list(graph.successors(end))
            predecessors = list(graph.predecessors(end))
        tips.append(path[::-1])
    path_lengths=[]
    avg_path_weights=[]
    for path in tips:
        path_lengths.append(len(path))
        avg_path_weights.append(path_average_weight(graph, path))
    graph = select_best_path(graph, tips, path_lengths, avg_path_weights,
                             delete_sink_node=True)
    return graph

def get_starting_nodes(graph):
	list_start=[]
	for i in graph.nodes:
		nodepre=list(graph.predecessors(i))
		if (nodepre ==[]):
			list_start.append(i)
	return list_start
	
def get_sink_nodes(graph):
	list_sortie=[]
	for i in graph.nodes:
		nodesuc=list(graph.successors(i))
		if (nodesuc ==[]):
			list_sortie.append(i)
	return list_sortie

def get_contigs(graph, starting_nodes, ending_nodes):
	tub=[]
	for start in starting_nodes:
		for end in ending_nodes:
			if nx.has_path(graph,start,end)== True:
				for path in nx.all_simple_paths(graph,start,end):
					tubu=""
					for i,word in enumerate(path):
						if (i==0):
							tubu=word
						else :
							tubu+=word[-1]
					tub.append([tubu, len(tubu)])
	return tub			
			
def save_contigs(contigs_list, output_file):
    with open (output_file, 'w') as fileout:
    	for i in range(0,len(contigs_list)):
    		fileout.write(">contig_"+str(i)+ " len=" + str(contigs_list[i][1]) +'\n'+ fill(contigs_list[i][0]) +'\n')


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def draw_graph(graph, graphimg_file):
    """Draw the graph
    """                                    
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt") as save:
            pickle.dump(graph, save)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    #print(args.kmer_size)
    #sequence=read_fastq(sys.argv[2])
    #print (sequence)
    #cutk=cut_kmer(sequence,args.kmer_size)
    #print (cutk)
    dico= build_kmer_dict(sys.argv[2], args.kmer_size)
    #print(dico)
    graph=build_graph(dico)
    #print(graph)
    start=get_starting_nodes(graph)
    end=get_sink_nodes(graph)
    path= get_contigs(graph,start, end)
    save_contigs(path, './res/resultat1.fasta')


    print('graph tip resolution')
    print('before')
    print(graph.number_of_nodes())
    print(graph.number_of_nodes())
    
    graph = simplify_bubbles(graph)
    graph = solve_entry_tips(graph, get_starting_nodes(graph))
    graph = solve_out_tips(graph, get_sink_nodes(graph))
    print('after')
    print(graph.number_of_nodes())
    print(graph.number_of_nodes())
    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)
    # Save the graph in file
    # if args.graph_file:
    #     save_graph(graph, args.graph_file)


if __name__ == '__main__':
    main()
