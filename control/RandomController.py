#!/usr/bin/env python

from plant.Grammar import Grammar
import plant.RNAdecomposition
from plant.BlockInverseFolding import BlockInverseFolding
import logging




logging.basicConfig(level=logging.INFO)


opts = {'rfam_id':'RF00005', 'episodes':100, 'log_file':'./BlockInverseFolding.log'}


class RandomController():
	"""
	For later use with learning with Q-table:
	qtable_init_name = "qtable_init.txt"
	epsilon  = 0.0
	alpha    = 1.0
	gamma    = 1.0
	training = true
	minimize = true
	this->qtable_save_freq = 100
	"""
	
	#Initializes the State.
	def __init__(self, **opts):
		self._number_of_episodes = opts['episodes']

"""
def do_episode()
def init_modules()
def get_initial_state()
def create_modules()
def read_options()
"""

"""
Modules:
1- Grammar -> rfam_id
2- BlockInverseFolding -> Grammar
3- 
"""

def observation_to_string( observation = None):
	str_observation = ''
	for ob in observation:
		str_observation = str_observation + str(ob) + ' '
	return str_observation


def save_observation_base(**opts):
	with open(opts['log_file'], 'a') as f:
		for observation in opts['observation_base']:
			f.write(observation)
			f.write('\n')


"""
All functions would be bound into one object in kernel
"""

# do episode
# return observations + action history
def do_episode(**opts):
	# Instantiate BlockInverseFolding
	BIF = BlockInverseFolding(Grammar = opts['grammar'])
	while True:
		flag = BIF.get_next_state()
		if not(flag) :
			break
	observation = observation_to_string( observation = BIF.get_observation() )
	biGramTable = BIF.get_biGramTable()
	print BIF.get_observation_vector()
	return observation , biGramTable


# 1 - create objects
# 2- initialize objects
# 3- initialize loop counter
# 4- repeat do_episode
# 5- save observation base to disk
def main_loop(**opts):
	# main loop for doing all episodes
	# Create an instance of Grammar object
	iGrammar = Grammar(opts['rfam_id'])
	opts.update({'grammar':iGrammar})
	opts.update({'observation_base':[]})
	opts.update({'biGramTable':[]})
	for episode in range(opts['episodes']):
		ob , bgt = do_episode(**opts)
		opts['observation_base'].append(ob)
		opts['biGramTable'].append(bgt)
	save_observation_base(**opts)


if __name__ == "__main__":
	
	logger = logging.getLogger(__name__)
	logger.info('Call to Random Controller module.')
	main_loop(**opts)


