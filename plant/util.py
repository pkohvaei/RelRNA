#!/usr/bin/env python
import logging
import numpy as np
from scipy.stats import norm
from eden.util.display import draw_graph, serialize_graph
from eden.modifier.graph import vertex_attributes
from itertools import islice
import networkx as nx
import pickle
import math


# Downloads the fasta file of family_if from rfam.
def rfam_uri(family_id):
	return 'http://rfam.xfam.org/family/%s/alignment?acc=%s&format=fastau&download=0'%(family_id,family_id)


# Fits a normal distribution on the given iterable.
def Gaussian(iterable = None):
	(mu,sigma) = norm.fit(iterable)
	return (mu,sigma)


def phi(x):
    return math.exp(-x * x / 2.0) / math.sqrt(2.0 * math.pi)


# Return the value of the cumulative Gaussian distribution function
# with mean 0.0 and standard deviation 1.0 at the given z value.
def Phi(z):
    if z < -8.0:
        return 0.0
    if z > 8.0:
        return 1.0
    sum = 0.0
    term = z
    i = 3
    while sum != sum + term:
        sum += term
        term *= z * z / float(i)
        i += 2
    return 0.5 + phi(z) * sum


# Return the value of the cumulative Gaussian distribution function
# with mean mu and standard deviation sigma at the given z value.
def cdf(z, mu=0.0, sigma=1.0):
    return Phi((z - mu) / sigma)


def plot(graphs, num=1):
	vis_opts={'size':14, 'node_border':False, 'node_size':200, 'font_size':9, 'vertex_alpha':0.6 , 'vertex_color':'weight', 'colormap':'autumn'}   
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


def priors( counter = None):
	sigma = 0
	for key in counter:
		sigma = sigma + counter[key]
	priors_dict = counter.copy()
	for key in priors_dict:
		priors_dict[key] = float(priors_dict[key]) / float(sigma)
	return priors_dict


def cdfs( prior_distributions = None , prior_counter = None , zero_substitute = 0.001 ):
	cdfs_dict = {}
	for key in prior_distributions:
		mu , sigma = prior_distributions[key]
		turn_out = prior_counter[key]
		if turn_out == 0 :
			turn_out = zero_substitute
		cdfs_dict.update({key:cdf(turn_out , mu = mu , sigma = sigma)})
	return cdfs_dict


def get_multinomial_dist( prior_distributions = None , prior_counter = None , zero_substitute = 0.001 ):
	probs = {}
	multinomial_dist = {}
	
	priors_dict = priors( counter = prior_counter )
	cdfs_dict = cdfs( prior_distributions = prior_distributions , prior_counter = prior_counter , zero_substitute = zero_substitute )
	
	for key in priors_dict:
		probs.update({key:priors_dict[key]*cdfs_dict[key]})
	actions = [key for key in probs]
	values = [probs[key] for key in probs]
	norm = [float(value)/sum(values) for value in values]
	
	norm_action_values = zip(actions,norm)
	
	for action , norm in norm_action_values:
		multinomial_dist.update({action:norm})
	
	return multinomial_dist


def sample_from_multinomial( multinomial_probs = None , sample= 1 , size = 1):
	turn_out = np.random.multinomial( sample , multinomial_probs , size = size )
	return turn_out.tolist()[0]


def soft_select(prior_distributions = None , prior_counter = None , zero_substitute = 0.001):
	multinomial_dist = get_multinomial_dist( prior_distributions = prior_distributions , prior_counter = prior_counter , zero_substitute = zero_substitute )
	actions = []
	probs = []
	for key in multinomial_dist:
		actions.append(key)
		probs.append(multinomial_dist[key])
	turn_out = sample_from_multinomial( multinomial_probs = probs)
	print actions
	print probs
	print turn_out
	soft_select_dict = {}
	selected_action = ''
	n = zip(actions , turn_out)

	for action , selected in n:
		if selected:
			selected_action = action
	return selected_action
	

logging.basicConfig(level=logging.INFO)


if __name__ == "__main__":
	
	logger = logging.getLogger(__name__)
	logger.info('Call to util module.')
	



