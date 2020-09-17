#!/usr/bin/env python

import os
import sys
import ConfigParser
import re
import random
import logging
import numpy
import itertools
import networkx as nx
from itertools import islice
from eden.converter.fasta import fasta_to_sequence
from eden.converter.rna.rnafold import rnafold_to_eden
import workSuit.tools as util


logger = logging.getLogger(__name__)


def basepairs(graph = None):
	"""
	Returns a list of tuples (x,y) of all basepairs in the graph.
	"""
	basepair_list = []
	for line in nx.generate_edgelist(graph):
		if ('basepair' in line):
			x = int(line.split(' ')[0])
			y = int(line.split(' ')[1])
			if not (((x,y) in basepair_list) or ((y,x) in basepair_list)):
				basepair_list.append((x,y))
	return basepair_list


def remove_stems(graph = None):
	"""
	Finds all basepairs in the graph and removes the participating nodes.
	"""
	for x,y in basepairs(graph = graph):
		graph.remove_node(x)
		graph.remove_node(y)


def successor_pair(basepair, start, stop):
	"""
	Calculates node numbers for the probable basepair
	which follows the given basepair in the graph.
	Returns the hypothetical calculated basepair. 
	"""
	x , y = basepair
	if (x + 1 > stop) or (y - 1 < start):
		return (-1,-1)
	else:
		return ( x + 1 , y - 1 )


def predecessor_pair(basepair, start, stop):
	"""
	Calculates node numbers for the probable basepair
	which is followed by the given basepair in the graph.
	Returns the hypothetical calculated basepair. 
	"""
	x , y = basepair
	if (x - 1 < start) or (y + 1 > stop):
		return (-1,-1)
	else:
		return ( x - 1 , y + 1 )


def break_points(graph = None):
	"""
	Finds and returns the list of all basepairs that act as "break point"s.
	A break point occurs when there is a change from one structural motif to another
	(i.e from a stem to a loop).
	"""
	start = 0
	stop = nx.number_of_nodes(graph) -1
	break_points = []
	basepair_list = basepairs(graph = graph)
	for basepair in basepair_list:
		if (successor_pair(basepair, start, stop) == (-1,-1)) or \
		(predecessor_pair(basepair, start, stop) == (-1,-1)) or \
		not (successor_pair(basepair, start, stop) in basepair_list) or \
		not (predecessor_pair(basepair, start, stop) in basepair_list):
			break_points.append(basepair)
	return break_points


def tag_break_points(graph = None):
	"""
	Finds and relabels all break points in the graph.
	"""
	break_points_list = break_points(graph = graph)
	for break_point in break_points_list:
		id1, id2 = break_point
		graph.edge[id1][id2]['label']= '$'
		graph.edge[id1][id2]['type']='breakpoint'


def reset_break_points(graph = None, pairs = None):

	for pair in pairs:
		id1, id2 = pair
		graph.edge[id1][id2]['label']= '='
		graph.edge[id1][id2]['type']='basepair'


def add_block(block1 = None, block2 = None):
	if util.is_null_graph(graph = block1):
		graph = block2
	else:
		if len(scan_magnetic_ends(graph = block1)) > 0 and len(scan_magnetic_ends(graph = block2)) > 0:
			graph = nx.disjoint_union(block1 , block2)
			component1 , component2 = list(nx.connected_component_subgraphs(graph))
			mg_end_list1 =  scan_magnetic_ends(graph = component1)
			mg_end_list2 =  scan_magnetic_ends(graph = component2)
			mg_end1 = random.choice(mg_end_list1)
			mg_end2 = random.choice(mg_end_list2)
			x,y = mg_end1
			w,z = mg_end2
			graph.add_edge(x,w,type='backbone',label='-')
			graph.add_edge(y,z,type='backbone',label='-')
			reset_break_points(graph = graph, pairs = [mg_end1,mg_end2])
		else:
			logger.debug('Error in add_block!! Premature magnetic-end outage.')
			return block1, False
	return graph, True

# *** bigram logic
#def add_block(block1 = None, block2 = None, open_end1 = None , open_end2 = None):
	
#	graph = nx.disjoint_union(block1 , block2)
#	x,y = open_end1
#	w,z = open_end1
#	graph.add_edge(x,w,type='backbone',label='-')
#	graph.add_edge(y,z,type='backbone',label='-')
#	reset_break_points(graph = graph, pairs = [open_end1,open_end2])

#	return graph, True



def stems(graph = None):
	"""
	Accepts a tagged graph as input.
	Outputs a list of stems each as a list.
	"""
	basepair_nodes_list = [x for x,y in basepairs(graph = graph)] + [y for x,y in basepairs(graph = graph)]
	graph_nodes_list = graph.nodes()
	non_stem_nodes_list = []

	for node in graph_nodes_list:
		if not(node in basepair_nodes_list):
			non_stem_nodes_list.append(node)
	g = graph.copy()
	for node in non_stem_nodes_list:
		g.remove_node(node)
	stems_list = []
	for subgraph in nx.connected_component_subgraphs(g):
		if subgraph.number_of_nodes() > 2:
			tag_stem_magnetic_ends(graph = subgraph)
			stems_list.append(subgraph)
	return stems_list


def loops(graph = None):
	"""
	Creates and returns a dictionary containing all loop types and their occurrences.
	"""
	unknown_structs = []
	compound_structs = []
	loops_dict = create_components_dict()
	for subgraph in nx.connected_component_subgraphs(graph):
		if subgraph.number_of_nodes() < 3:
			unknown_structs.append(subgraph)
		else:
			if connectivity_threshold(graph = subgraph) > 2 or loop_type(graph= subgraph) == 'NA':
				compound_structs.append(subgraph)
			else:
				loops_dict[loop_type(graph= subgraph)].append(subgraph)
	return loops_dict


def loop_type(graph = None):
	"""
	Accepts a graph of unknown non-stem strusture.
	Returns the loop type or dangling-end.
	"""
	breakpoints = 0
	for edge in nx.generate_edgelist(graph):
		if ('breakpoint' in edge):
			breakpoints += 1

	if breakpoints == 1:
		if is_loop(graph = graph):
			return 'hairpinloop'
		else:
			return 'dangling_end'

	if is_loop(graph = graph):
		if breakpoints == 2:
			return multiloop_type(graph = graph)
		elif breakpoints == 3:
			return 'multiloop3'
		elif breakpoints == 4:
			return 'multiloop4'
		elif breakpoints == 5:
			return 'multiloop5'
	return 'NA'


def multiloop_type(graph = None):
	"""
	Distingushes bulges from internal_loops and outputs the result.
	"""
	g = graph.copy()
	breakpoint_nodes = []
	for edge in g.edges(data=True):
		id1 ,id2 , data = edge
		if g.edge[id1][id2]['label']== '$':
			if not id1 in breakpoint_nodes:
				breakpoint_nodes.append(id1)
			if not id2 in breakpoint_nodes:
				breakpoint_nodes.append(id2)
	for node in breakpoint_nodes:
		g.remove_node(node)
	if nx.number_connected_components(g) ==1:
		return 'bulge'
	else:
		return 'internal_loop'


def calculate_nodes_connectivitiy(graph = None):
	node_branch = {}
	for node in graph.nodes():
		node_branch.update({node:0})
	for line in nx.generate_edgelist(graph):
		node_branch[int(line.split(' ')[0])] += 1 
		node_branch[int(line.split(' ')[1])] += 1
	return node_branch


def connectivity_threshold(graph = None):
	node_branch = calculate_nodes_connectivitiy(graph = graph)
	threshold = 0
	for node in node_branch:
		if node_branch[node] > threshold:
			threshold = node_branch[node]
	return threshold


def is_loop(graph = None):
	node_branch = calculate_nodes_connectivitiy(graph = graph)
	connectivity = 2
	for node in node_branch:
		if not( node_branch[node] == connectivity):
			return False
	return True


def create_components_dict():
	components_dict = {'hairpinloop':[], 'bulge':[], 'internal_loop':[], 'multiloop3':[], 'multiloop4':[], 'multiloop5':[], 'dangling_end':[], 'stem':[]}
	return components_dict


def create_counter_dict():
	counter_dict = {'hairpinloop':[], 'bulge':[], 'internal_loop':[], 'multiloop3':[], 'multiloop4':[], 'multiloop5':[], 'dangling_end':[], 'size':[], 'stem':[]}
	return counter_dict


def graph_decomposition(graph = None):
	"""
	Outputs a dictionary of all motifs found in the graph.
	Each motif entry keeps all of the kind as individual lists of subgraphs.
	The ID of the original sequence is included.
	The empty dictionary looks as follows:

	"""
	g = graph.copy()
	tag_break_points(graph = g)
	stems_list = stems(graph = g)
	remove_stems(graph = g)

	component_dict = loops(graph = g)
	component_dict['stem'] = stems_list

	hairpinloops = len(component_dict['hairpinloop'])
	bulges = len(component_dict['bulge'])
	internal_loops = len(component_dict['internal_loop'])
	multiloop3s = len(component_dict['multiloop3'])
	multiloop4s = len(component_dict['multiloop4'])
	multiloop5s = len(component_dict['multiloop5'])
	dangling_ends = len(component_dict['dangling_end'])
	stem = len(component_dict['stem'])
	seq_length = nx.number_of_nodes(graph)

	component_counter = {'hairpinloop':hairpinloops , 'bulge':bulges , 'internal_loop':internal_loops , \
						'multiloop3':multiloop3s , 'multiloop4':multiloop4s , 'multiloop5':multiloop5s, 'dangling_end':dangling_ends , \
						'stem':stem, 'size':seq_length}
	return component_dict, component_counter


def global_decomposition(iterable = None):
	"""
	Main function for decomposing all members in iterable.
	Generator of RNA-motif dictionaries over a given rfam dataset.
	"""
	for graph in iterable:
		yield graph_decomposition(graph = graph)


def scan_magnetic_ends(graph = None):

	magnetic_ends = []
	for line in nx.generate_edgelist(graph):
		if ('$' in line):
			x = int(line.split(' ')[0])
			y = int(line.split(' ')[1])
			magnetic_ends.append((x,y))
	return magnetic_ends


def count_magnetic_ends(graph = None):
	return len(scan_magnetic_ends(graph = graph))


def find_stem_magnetic_ends(graph = None):
	node_branch = {}
	for node in graph.nodes():
		node_branch.update({node:0})
	for line in nx.generate_edgelist(graph):
		if ('backbone' in line):
			node_branch[int(line.split(' ')[0])] += 1 
			node_branch[int(line.split(' ')[1])] += 1 
	ends = []
	for node in node_branch:
		if node_branch[node] == 1:
			if not(node in ends):
				ends.append(node)
	magnetic_ends = []
	for x,y in basepairs(graph = graph):
		if x in ends and y in ends:
			magnetic_ends.append((x,y))

	return magnetic_ends


def tag_stem_magnetic_ends(graph = None):
	"""
	Finds and relabels all magnetic ends in the stem.
	"""
	magnetic_ends = find_stem_magnetic_ends(graph = graph)
	for magnetic_end in magnetic_ends:
		id1, id2 = magnetic_end
		graph.edge[id1][id2]['label']= '$'
		graph.edge[id1][id2]['type']='breakpoint'
	return


def reset_unimportant_nodes(iterable = None, threshold = None, label = None):
	"""
	Accepts annotated graphs as input.
	Resets labels of the nodes with importance below 'threshold' to label.
	"""
	for graph in iterable:
		for node, data in graph.nodes_iter(data=True):
			if float(data['importance']) < threshold:
				data['label'] = label
		yield graph


if __name__ == "__main__":
	
	logging.basicConfig(level=logging.INFO)
	logger = logging.getLogger(__name__)
	logger.info('Call to RNAmotifs module.')

	#iterable = itertools.islice(fasta_to_sequence(rfam_uri('RF00005')),20)
	iterable = fasta_to_sequence(util.rfam_uri('RF00005'))
	graphs = rnafold_to_eden(iterable)
	iterable = RNAmotifs(rfam_id  = 'RF00005')
	for item in iterable:
		print item

