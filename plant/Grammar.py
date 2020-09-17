#!/usr/bin/env python

import logging
import numpy
import itertools
import networkx as nx
import plant.RNAdecomposition as rd
from eden.converter.fasta import fasta_to_sequence
from eden.converter.rna.rnafold import rnafold_to_eden
import workSuit.tools as util
import workSuit.probs as prob
import random


logger = logging.getLogger(__name__)


def produce_Grammar(rfam_id = None):
	"""
	Generates Grammar for a given rfam family.
	"""
	iterable = fasta_to_sequence('/home/kohvaeip/RLS/kernel/RF00005.fa') #util.rfam_uri(rfam_id))
	iterable = rnafold_to_eden(iterable)

	iterable = rd.global_decomposition(iterable = iterable)

	Grammar_dict = rd.create_components_dict()
	Raw_Statistics_dict = rd.create_counter_dict()
	
	for components, counters in iterable:
		for component in components:
			Grammar_dict[component] = Grammar_dict[component] + components[component]
		for counter in counters:
			Raw_Statistics_dict[counter].append(counters[counter])
	Statistics_dict = {}
	for key in Raw_Statistics_dict:
		Statistics_dict.update({key: prob.Gaussian(iterable = Raw_Statistics_dict[key])})

	return Grammar_dict, Statistics_dict


def generate_bigram_table( Grammar = None):
	biGram_table = {}
	for key1 in Grammar:
		for key2 in Grammar:
			if not ((key2,key1) in biGram_table):
				biGram_table.update({(key1,key2):0})
	return biGram_table


def update_bigram_table(  key1 , key2 , val , biGramTable = None ):
	if (key1,key2) in biGramTable:
		val = biGramTable[(key1,key2)] + val
		biGramTable.update({(key1,key2):val})
	elif (key2,key1) in biGramTable:
		val = biGramTable[(key2,key1)] + val
		biGramTable.update({(key2,key1):val})
	else:
		sys.stdout.write('Error in updating the biGram table!')


def biGramTable_to_vector( biGramTable = None ):
	biGram_list = []
	biGrams = [('dangling_end', 'stem') , ('multiloop5', 'dangling_end') , ('internal_loop', 'internal_loop') , \
	 ('bulge', 'hairpinloop') , ('dangling_end', 'dangling_end') , ('multiloop3', 'stem') , \
	 ('hairpinloop', 'multiloop4') , ('multiloop5', 'multiloop5') , ('internal_loop', 'hairpinloop') , \
	 ('bulge', 'multiloop5') , ('bulge', 'multiloop3') , ('internal_loop', 'multiloop4') , \
	 ('multiloop3', 'internal_loop') , ('multiloop3', 'dangling_end') , ('multiloop3', 'multiloop4') , \
	 ('bulge', 'stem') , ('multiloop4', 'multiloop4') , ('dangling_end', 'multiloop4') , \
	 ('internal_loop', 'stem') , ('internal_loop', 'dangling_end') , ('bulge', 'bulge') , \
	 ('multiloop5', 'multiloop4') , ('bulge', 'internal_loop') , ('hairpinloop', 'dangling_end') , \
	 ('stem', 'multiloop4') , ('bulge', 'multiloop4') , ('bulge', 'dangling_end') , ('multiloop5', 'internal_loop') ,\
	 ('hairpinloop', 'hairpinloop') , ('hairpinloop', 'stem') , ('multiloop3', 'multiloop5') , \
	 ('multiloop3', 'multiloop3') , ('multiloop3', 'hairpinloop') , ('multiloop5', 'stem') , ('stem', 'stem') , \
	 ('multiloop5', 'hairpinloop')]

	for biGram in biGrams:
		biGram_list.append(biGramTable[biGram])
	return tuple(biGram_list)


class Grammar():
	"""
	Grammar is a class containing two main objects:
	1- a dictionary which holds RNA 2D-components as keys and a list of all generated graph components as values
	2- a dictionary which holds RNA 2D-components and RNAsequence size as keys and the Gaussian model of the 
	occurrence probability of each in the form of mean-covarinace pairs as values
	Initializes an instance of the Grammar.
	"""
	def __init__(self, rfam_id):
		self.Words, self.Statistics = produce_Grammar(rfam_id = rfam_id)
		logger.info('Initialized an instance of Grammar for %s family.' %rfam_id)

	def query_chapter(self, component = None):
		"""
		method for queriying words from chapters in Grammar.
		"""
		graph = nx.Graph()
		graph = random.choice(self.Words[component])
		return graph
		
	def select_scan_component(self, component = None):
		"""
		Method for queriying words from chapters in Grammar.
		Returns a graph of type component and the list of its magnetic ends.
		"""
		graph = nx.Graph()
		graph = random.choice(self.Words[component])
		open_ends = rd.scan_magnetic_ends(graph = graph)
		return graph , open_ends


if __name__ == "__main__":

	logger = logging.getLogger(__name__)
	logger.info('Call to Grammar module.')
	rfam_id = 'RF00005'
	bgt = {('dangling_end', 'stem'): 0, ('multiloop5', 'dangling_end'): 2, ('internal_loop', 'internal_loop'): 1, ('bulge', 'hairpinloop'): 0, ('dangling_end', 'dangling_end'): 0, ('multiloop3', 'stem'): 3, ('hairpinloop', 'multiloop4'): 0, ('multiloop5', 'multiloop5'): 1, ('internal_loop', 'hairpinloop'): 1, ('bulge', 'multiloop5'): 1, ('bulge', 'multiloop3'): 2, ('internal_loop', 'multiloop4'): 1, ('multiloop3', 'internal_loop'): 2, ('multiloop3', 'dangling_end'): 1, ('multiloop3', 'multiloop4'): 3, ('bulge', 'stem'): 3, ('multiloop4', 'multiloop4'): 1, ('dangling_end', 'multiloop4'): 6, ('internal_loop', 'stem'): 0, ('internal_loop', 'dangling_end'): 0, ('bulge', 'bulge'): 0, ('multiloop5', 'multiloop4'): 3, ('bulge', 'internal_loop'): 0, ('hairpinloop', 'dangling_end'): 1, ('stem', 'multiloop4'): 2, ('bulge', 'multiloop4'): 4, ('bulge', 'dangling_end'): 0, ('multiloop5', 'internal_loop'): 2, ('hairpinloop', 'hairpinloop'): 0, ('hairpinloop', 'stem'): 3, ('multiloop3', 'multiloop5'): 3, ('multiloop3', 'multiloop3'): 0, ('multiloop3', 'hairpinloop'): 2, ('multiloop5', 'stem'): 0, ('stem', 'stem'): 0, ('multiloop5', 'hairpinloop'): 1}
	print biGramTable_to_vector( biGramTable = bgt )
	iGrammar = Grammar(rfam_id)
	print type(iGrammar.query_chapter(component = 'stem'))

