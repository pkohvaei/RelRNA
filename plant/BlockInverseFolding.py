import os
import sys
import logging
from math import ceil
import networkx as nx
import random
from plant.Grammar import Grammar
from plant.Grammar import generate_bigram_table
from plant.Grammar import update_bigram_table , biGramTable_to_vector
from plant.RNAdecomposition import add_block
import plant.RNAdecomposition as rnadc
import workSuit.tools as util
import workSuit.probs as prob
from control.Params import Params

"""
State representation:
1- networkx graph
2- component counter: dictionary of elements and their permitted number of occurrence: 
   this is extracted from components's Gaussian distribution contained in the Grammar.
3- biGram component connectivity table.
"""

logger = logging.getLogger(__name__)


def sample_counter_from_gaussian(grammar_statistics_dict = None , component = None):
	mu, sigma = grammar_statistics_dict[component]
	sample_member = random.gauss(mu, sigma)
	if sample_member < 0 :
		return 0
	else:
		return int(ceil(sample_member))


class State():
	#Initializes the State.
	def __init__(self, Grammar):
		self._graph = nx.Graph()
		self.biGramTable = generate_bigram_table( Grammar = Grammar.Words)
		self.previous_component = ''
		# *** bigram logic
		# self.open_ends = []

		n_hairpinloop = sample_counter_from_gaussian(grammar_statistics_dict = Grammar.Statistics , component = 'hairpinloop')
		n_bulge = sample_counter_from_gaussian(grammar_statistics_dict = Grammar.Statistics , component = 'bulge')
		n_internal_loop = sample_counter_from_gaussian(grammar_statistics_dict = Grammar.Statistics , component = 'internal_loop')
		n_multiloop3 = sample_counter_from_gaussian(grammar_statistics_dict = Grammar.Statistics , component = 'multiloop3')
		n_multiloop4 = sample_counter_from_gaussian(grammar_statistics_dict = Grammar.Statistics , component = 'multiloop4')
		n_multiloop5 = sample_counter_from_gaussian(grammar_statistics_dict = Grammar.Statistics , component = 'multiloop5')
		n_dangling_end = sample_counter_from_gaussian(grammar_statistics_dict = Grammar.Statistics , component = 'dangling_end')
		if n_dangling_end > 1:
			n_dangling_end = 1
		n_stem = sample_counter_from_gaussian(grammar_statistics_dict = Grammar.Statistics , component = 'stem')

		self._seq_length_threshold = sample_counter_from_gaussian(grammar_statistics_dict = Grammar.Statistics , component = 'size')
		self._component_counter = {'hairpinloop':n_hairpinloop, 'bulge':n_bulge, 'internal_loop':n_internal_loop, 'multiloop3':n_multiloop3, \
									'multiloop4':n_multiloop4, 'multiloop5':n_multiloop5, 'dangling_end':n_dangling_end , 'stem':n_stem}

	def get_state_counters(self):
		return self._component_counter, self._graph.number_of_nodes()

	def _check_initial_state(self, Grammar):
		#Checks for a valid initial state.
		# *** bigram logic 
		#block2 , block2_open_ends = Grammar.query_chapter( component = 'stem' )
		block2 = random.choice(Grammar.Words['stem'])
		self._graph, self._composability = add_block(block1 = self._graph , block2 = block2)
		#  *** bigram logic
		# self._graph, self._composability , new_open_ends = add_block(block1 = self._graph , block2 = block2)
		# self.open_ends = self.open_ends + new_open_ends
		self.update_component_counter('stem')
		self.previous_component = 'stem'
		#Checks if dangling end is allowed.
		#Adds a randomely-chosen dangling end to the initial graph.
		if self._component_counter['dangling_end'] > 0:
			# *** bigram logic
			#block2 , block2_open_ends = Grammar.query_chapter( component = 'dangling_end' )
			self._graph, self._composability = add_block(block1 = self._graph , block2 = random.choice(Grammar.Words['dangling_end']))
			self.update_component_counter('dangling_end')
			update_bigram_table('dangling_end' , 'stem' , 1 , biGramTable = self.biGramTable)
		

	def update_component_counter(self, component):
		if self._component_counter[component] > 0:
			self._component_counter[component] = self._component_counter[component] - 1
		else:
			pass
			# Raise exception

	def check_state_graph_length(self):
		if self._graph.number_of_nodes() < self._seq_length_threshold:
			return False
		else:
			return True

	def component_counter_sum(self):
		cp_sum = 0
		for key in self._component_counter:
			cp_sum = cp_sum + self._component_counter[key]
		return cp_sum

	def available_components(self):
		available_component_list = []
		for key in self._component_counter:
			if self._component_counter[key] > 0:
				available_component_list.append(key)
		return available_component_list

	def get_state_vector(self):
		state_vector = biGramTable_to_vector( biGramTable = self.biGramTable )
		open_magnetic_ends = (rnadc.count_magnetic_ends(graph = self._graph),)
		return state_vector + open_magnetic_ends


# *** bigram logic
#def merge_graphs( opts , graph1 = None , graph2 = None , open_end1 = None, open_end2 = None ):
	#
	#self._state._graph, self._composability = add_block(block1 = self._state._graph , block2 = component)




class BlockInverseFolding():

	# Initialize Block Inverse Folding plant.
	def __init__(self, Grammar = None):
		# Initializes the Grammer.
		self._Grammar = Grammar
		# Initializes the state.
		self._state = None
		self._composability = None
		self._action_history = None
		logger.info('An instance of BlockInverseFolding environment initialized.')
		
	# Repeats at the beginning of each episode.
	def notify_episode_starts(self , opts):
		# Initializes class variables and objects which are used episode-wise.
		self._state = State(self._Grammar)
		previous_state = self.get_observation_vector()
		opts.update({'previous_state':previous_state})
		self._state._check_initial_state(self._Grammar)
		current_state = self.get_observation_vector()
		opts.update({'current_state':current_state})
		opts.update({'action':'stem'})	
		self._composability = True
		self._action_history = []
		
	# Repeats at the end  of each episode.
	def notify_episode_stops(self , opts):
		if opts['visualize']:
			opts['products'].append(self._state._graph)
			graphs = [self._state._graph]
			graphs = util.pre_process(graphs)
			util.plot(graphs, num=1)
		return

	#Determins if the given plant_state is a terminal state and optionally 
	#returns terminal costs. In the Block Inverse Folding plant, 
	#terminal costs are the covariance score of the folded sequence.
	def is_terminal(self):
		if (self._state.check_state_graph_length()) or (self._state.component_counter_sum() <= 0) or (self._composability == False):
			return True
		else:
			return False

	# Computes the next state of the environment, given a current state and an action (state transition).
	# Returns True if the next state is non-terminal.
	def get_next_state(self , opts):
		
		if self.is_terminal():
			return False
		else:
			component_name = opts['action']
			# *** bigram logic
			#component , component_open_ends = self._Grammar.query_chapter( component = component_name )
			component = self._Grammar.query_chapter( component = component_name )
			# Adds the component to state graph.
			self._state._graph, self._composability = add_block(block1 = self._state._graph , block2 = component)
			# Updates the state counter.
			self._state.update_component_counter(component_name)
			update_bigram_table(self._state.previous_component , component_name , 1 , biGramTable = self._state.biGramTable)
			self._state.previous_component = component_name
			self._action_history.append(component_name)
			previous_state = opts['current_state']
			current_state = self.get_observation_vector()
			opts.update({'previous_state':previous_state})
			opts.update({'current_state':current_state})
		return True

	def get_observation_vector(self):
		state_vector = self._state.get_state_vector()
		return state_vector

	def available_actions(self):
		return self._state.available_components()

	"""
	def random_sample_component(self):
		
		#Plant random action selector, not used with other controllers.
		
		available_components_list = self._state.available_components()
		component = random.choice(available_components_list)
		return component, random.choice(self._Grammar.Words[component])
	"""

	def get_observation(self , reward = False):
		"""
		Mainly for debugging purposes and reward calculation
		"""
		# magnetic_ends ****
		open_magnetic_ends = rnadc.count_magnetic_ends(graph = self._state._graph)
		# sequence_length difference factor
		instance_counter, seq_length = self._state.get_state_counters()
		seq_length_difference = abs(seq_length - self._state._seq_length_threshold)
		# composability factor
		if reward:
			return open_magnetic_ends , seq_length_difference
		else:
			return [open_magnetic_ends] + [seq_length_difference] + self.get_action_history()

	def get_action_history(self):
		"""
		Action trajectory recorder for Random Controller
		"""
		return self._action_history

	def get_biGramTable(self):
		"""
		Mainly for debugging purposes
		"""
		return self._state.biGramTable


if __name__ == "__main__":
	
	logging.basicConfig(level=logging.INFO)
	logger = logging.getLogger(__name__)
	logger.info('Call to Block Inverse Folding plant module.')
	
	grammar = Grammar(rfam_id ='RF00005')

	bif = BlockInverseFolding(grammar)

	

