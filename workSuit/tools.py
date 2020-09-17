#!/usr/bin/env python
import logging
from eden.util.display import draw_graph, serialize_graph
from eden.modifier.graph import vertex_attributes
from itertools import islice
import networkx as nx
import pickle


# Downloads the fasta file of family_if from rfam.
def rfam_uri(family_id):
	return 'http://rfam.xfam.org/family/%s/alignment?acc=%s&format=fastau&download=0'%(family_id,family_id)


def plot(graphs, num=1):
	vis_opts={'size':10, 'node_border':False, 'node_size':200, 'font_size':9, 'vertex_alpha':0.6 , 'vertex_color':'weight', 'colormap':'autumn'}   
	graphs = islice(graphs,num)
	for graph in graphs: draw_graph(graph, **vis_opts)
	

def pre_process(graphs):
    graphs = vertex_attributes.colorize(graphs, output_attribute = 'level', labels = ['A','U','C','G'])
    return graphs


def is_null_graph(graph = None):
	if nx.number_of_nodes(graph) == 0:
		return True
	else:
		return False

def exp_growth( a , b , c , x ):
	return a * ( b ** (c * x))


def save_pickle( struct = None , path = None):
	try:
		afile = open(path, 'wb')
		pickle.dump(struct , afile)
		afile.close()
	except pickle.PicklingError:
		print 'Error when serializing data.'


def load_pickle(path = None):
	try:
		file = open(path, 'rb')
		struct = pickle.load(file)
		file.close()
		return struct
	except pickle.PicklingError:
		print 'Error when opening pickle file %s' %path


logging.basicConfig(level=logging.INFO)


if __name__ == "__main__":
	
	logger = logging.getLogger(__name__)
	logger.info('Call to toolSuit module.')
	



