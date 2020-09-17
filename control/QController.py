#!/usr/bin/env python

from plant.Grammar import Grammar
from plant.BlockInverseFolding import BlockInverseFolding
import logging
import random
import workSuit.tools as util
import workSuit.probs as prob
from control.Params import Params
from math import exp


"""
Implementation comments:
	Q-table as a dictionary with (state,action) tuples as keys
	A method for encoding states and assigning IDs
	A method for assigning IDs to actions
"""

logger = logging.getLogger(__name__)


class QTable():
	def __init__(self , opts):
		self.QT = {}
		if opts['load']:
			self.load(opts['qt_save_path'])
			logger.info('An instance of Q-Table initialized with from file: %s' %opts['qt_save_path'])
		else :
			logger.info('An instance of Q-Table initialized.')

	def update( self , state , action , value ):
		state_action = (state,action)
		self.QT.update({state_action:value})

	def get_value( self , state , action , opts ):
		if not (state , action) in self.QT:
			q_value = opts['initial_q_value']
			self.update( state , action , q_value)
			return opts['initial_q_value']
		else:
			return self.QT[(state , action)]

	def load( self , path):
		self.QT = util.load_pickle( path = path )

	def save( self , opts ):
		path = opts['qt_save_path']
		util.save_pickle( struct = self.QT , path = path)


class QTableController():
	def __init__(self, opts):

		# Initilaizes learning paramters.
		epsilon = opts['epsilon']
		alpha = opts['alpha']
		gamma = opts['gamma']

		# Initilaizes the QTable.
		self._QT = opts['q_table']
		logger.info('An instance of Q-Table controller initialized.')

		
	def notify_episode_starts(self):
		return

	def get_action(self , **opts):
		plant = opts['plant']
		state = opts['current_state']
		epsilon = opts['epsilon']
		
		action = ''
		available_actions = plant.available_actions()
		
		if(available_actions):
			if (epsilon == 0 or epsilon < random.uniform(0,1)):
				action = self.get_opt_action(state , available_actions)
			else:
				action = random.choice(available_actions)
		return action

	def read_options():
		"""
		For the time being not used.
		"""
		return

	def notify_transition( self , reward , opts ):
		plant = opts['plant']
		training = opts['training']
		if training :
			current_state = opts['current_state']
			previous_state = opts['previous_state']
			action = opts['action']
			is_terminal = plant.is_terminal()
			self.train( reward , is_terminal , opts )
		return

	def notify_episode_stops(self):
		return

	def train( self , reward , is_terminal , opts ):
		
		gamma = opts['gamma']
		alpha = opts['alpha']
		this_state = opts['previous_state']
		next_state = opts['previous_state']
		action = opts['action']
		# Fetch the current Q value from table (action, state).
		current_q_value = self._QT.get_value( this_state , action , opts )
		
		if is_terminal:
			new_q_value = reward
		else:
			max_q_next = self.get_opt_q_value( next_state , opts )
			new_q_value = gamma * max_q_next + reward
		new_q_value = (1-alpha) * current_q_value + alpha * new_q_value
		# Update Q-table with new_q_value.
		self._QT.update( this_state , action , new_q_value )
		q_diff = new_q_value - current_q_value
		return q_diff

	def get_opt_action(self , state , actions):
		action_value = []
		for action in actions:
			if (state , action) in self._QT.QT:
				action_value.append((action , self._QT.QT[( state , action)]))

		if not(action_value):
			# IMPORTANT : replace with probabilistic version
			return random.choice(actions)
		else:
			action , value = sorted(action_value , key=lambda x: x[1], reverse=True)[0]
			return action

	def get_opt_q_value( self , state , opts ):
		state_actions = [key for key in self._QT.QT if key[0] == state]
		q_values = [self._QT.QT[state_action] for state_action in state_actions]
		if not(q_values):
			return opts['initial_q_value']
		else:
			return sorted(q_values , reverse=True)[0]


def consensus_variance( flip = False , **opts ):
	
	consensusBigram = opts['consensusBigram']
	stateBigram = opts['current_state']
	maskBigram = []
	
	if not(flip) :
		for x,y in zip(list(consensusBigram),list(stateBigram)):
			if x != 0:
				maskBigram.append(abs(x-y))
			else:
				maskBigram.append(0)
	else:
		for x,y in zip(list(consensusBigram),list(stateBigram)):
			if x == 0:
				maskBigram.append(abs(x-y))
			else:
				maskBigram.append(0)
	return maskBigram


def get_reward( **opts ):
	consensus_var = consensus_variance(flip = False , **opts)
	plant = opts['plant']
	default_reward = opts['default_reward']
	
	if plant.is_terminal():
		open_ends , seq_length_difference = plant.get_observation(reward = True)
		reward = (2 / exp(sum(consensus_var))) 
		# - ( util.exp_growth( 0.5 , 2 , 0.5 , open_ends )) # + util.exp_growth( 0.5 , 2 , 0.5 , seq_length_difference ))
		return reward
	else:
		return default_reward


if __name__ == "__main__":
	
	logger = logging.getLogger(__name__)
	logger.info('Call to QTableController module.')

	iGrammar = Grammar('RF00005')
	BIF = BlockInverseFolding( Grammar = iGrammar)
	params = Params()
	opts = params.opts
	opts.update({'plant':BIF})
	qtc = QTableController(**opts)

