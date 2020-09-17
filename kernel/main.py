#!/usr/bin/env python

from plant.Grammar import Grammar
from plant.BlockInverseFolding import BlockInverseFolding
import logging
import random
#import workSuit.tools as util
#import workSuit.probs as prob
from control.Params import Params
from control.QController import QTable
from control.QController import QTableController
from control.QController import get_reward
import time


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

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


# do episode
# return observations + action history
def do_episode( opts ):
	
	# episode_wise variables
	plant = opts['plant']
	controller = opts['controller']
	
	# plant notify_episode_starts
	plant.notify_episode_starts( opts )
	# controller notify_episode_starts
	controller.notify_episode_starts()
	# reward notify_episode_starts
	
	while True:
		action = ''
		# Calculate the next action.
		opts.update({'action': controller.get_action(**opts)})
		# Calculate the next state.
		flag = plant.get_next_state( opts )
		if not(flag) :
			break
		else:
			# Calculate the reward.
			reward = get_reward( **opts )
			# Record transition and train.
			controller.notify_transition( reward , opts )
			
	# plant notify_episode_stops
	plant.notify_episode_stops(opts)
	# controller notify_episode_stops
	controller.notify_episode_stops()
	# reward notify_episode_stops



def main_loop():
	"""
	Creates and initializes grammar and q-table objects.
	Creates and initializes plant and controller and params modules.
	Runs the main loop and repeats episodes.
	Saves the Q-Table and policy map.
	"""
	start_stamp = time.time()
	
	logger.info('Starting Block Inverse Folding learning experiment...')

	# Create an instance of parameter container object.
	params = Params()
	opts = params.opts

	# Fetch system specific parameters from the container.
	episodes = opts['episodes']
	rfam_id = opts['rfam_id']
	
	# Create an instance of Q-Table object.
	q_table = QTable( opts )
	opts.update({'q_table':q_table})
	
	# Create an instance of grammar object.
	grammar = Grammar(rfam_id)
	opts.update({'grammar':grammar})
	
	# Create an instance of plant.
	plant = BlockInverseFolding( Grammar = grammar )
	opts.update({'plant':plant})
	
	# Create an instance of controller object.
	controller = QTableController( opts )
	opts.update({'controller':controller})

	for episode in range(episodes):
		logger.info('Episode %d started.' %episode)
		do_episode( opts )
	
	# Final save operations.
	if opts['save']:
		q_table.save( opts )
		logger.info('Q-Table saved in %s' %opts['qt_save_path'])
		
	elapsed_time = time.time() - start_stamp
	
	logger.info('End of experiment.')
	logger.info('%d episodes took %f seconds to run.' %(episodes , elapsed_time))
		
	return opts


if __name__ == "__main__":
	
	logger = logging.getLogger(__name__)
	logger.info('Call to main module.')

	opts = main_loop()
