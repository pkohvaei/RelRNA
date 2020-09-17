#!/usr/bin/env python
import logging
import numpy as np
from scipy.stats import norm
import random
from itertools import islice
import numpy as np
import math


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

	

	
	
cp_counter = {'stem':4 , 'multiloop':6 , 'hairpin':7 , 'bulge':3}

cd_dist = {'stem':(3,2) , 'multiloop':(2,1) , 'hairpin':(3,1) , 'bulge':(3,1.5)}


logging.basicConfig(level=logging.INFO)


if __name__ == "__main__":
	
	logger = logging.getLogger(__name__)
	logger.info('Call to softSelect module.')
	
	cp = [key for key in cp_counter]
	
	#print cp
	#print priors( counter = cp_counter )
	#print cdfs(distributions = cd_dist , counter = cp_counter)
		
	print soft_select(prior_distributions = cd_dist , prior_counter = cp_counter , zero_substitute = 0.001)
	
	


