#!/usr/bin/env python

import logging


logger = logging.getLogger(__name__)


actions = ['multiloop4' , 'bulge' , 'internal_loop' , 'stem' , 'multiloop5' , 'multiloop3', 'hairpinloop']

consensusBigram = (1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 4 , 0 , 0 , 0 , 0 , 3 , 0 , 0 , 0 , 0 , 0 , 0)


""" Train
options = { 
		# *** bigram logic
		# *******frontier should be added to params
		# Options which are initialized once and remain fixed:
		'rfam_id':'RF00005', 'episodes':300000 , 'log_file':'./BlockInverseFolding.log', 'epsilon':0.5 , 
		'alpha':1.0 , 'gamma':1.0  , 'default_reward': 0.0 , 'training':True , 'minmize':True , 
		'qt_save_frequency':100 , 'initial_q_value':0.0 , 'save':True  , 'qt_save_path':'./QTable.txt' , 
		'action_set': actions , 'load':True , 'visualize':False , 'consensusBigram':consensusBigram ,
		# Options which are initialized once and change persist through the whole experiment:
		'products':[] ,
		# Options which are initialized once from inside the main loop:
		'plant':None , 'controller':None , 'q_table':None , 'grammar':None , 
		# Options which are initialized at the beginning of each episode:
		'observation_base':None , 'action':None , 'current_state':None , 'previous_state':None 
		}
"""

		
class Params():
	def __init__(self):
		self.opts = options
		logger.info('An instance of parameter container initialized.')
		return


if __name__ == "__main__":
	
	logger = logging.getLogger(__name__)
	logger.info('Call to Params module.')
	params = Params()
	print params.opts
	print params.actions

#""" Test

options = { 
		#  *** bigram logic
		# *******frontier should be added to params
		# Options which are initialized once and remain fixed:
		'rfam_id':'RF00005', 'episodes':25 , 'log_file':'./BlockInverseFolding.log', 'epsilon':0.0 , 
		'alpha':1.0 , 'gamma':1.0  , 'default_reward': 0.0 , 'training':False , 'minmize':True , 
		'qt_save_frequency':100 , 'initial_q_value':0.0 , 'save':False  , 'qt_save_path':'./QTable.txt' , 
		'action_set': actions , 'load':True , 'visualize':True , 'consensusBigram':consensusBigram ,
		# Options which are initialized once and change persist through the whole experiment:
		'products':[] , 
		# Options which are initialized once from inside the main loop:
		'plant':None , 'controller':None , 'q_table':None , 'grammar':None , 
		# Options which are initialized at the beginning of each episode:
		'observation_base':None , 'action':None , 'current_state':None , 'previous_state':None 
		}
#"""
#/home/kohvaeip/RLS/kernel/QT/Version 1
