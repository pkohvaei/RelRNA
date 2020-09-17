#!/usr/bin/env python

import logging

from antaParams import antaParams
from eden.converter.rna.rnashapes import rnashapes_to_eden
from eden import graph
import ConstraintFinder as cf
import antaRNA_v109


def DesignRNA(param_file=None , iterable=None, vectorizer=None, estimator=None, \
			  Cstr_threshold=0, Cstr_adjaceny=1, DotBr_threshold=0, DotBr_adjacency=1, unpaired_adjacency=1, sequence_counter=1):
	"""
	Generator for synthesizing RNA sequences.
	Takes as input an iterator over networkx graphs and outputs an iterator 
	over fasta-seq,original-fasta-id pairs.
	Uses antaRNA for sequence synthesis and EDeN for annotating networkx graphs.
	"""
	op = antaParams(param_file)
	iterable = vectorizer.annotate(iterable, estimator=estimator)
	iterable = cf.generate_antaRNA_constraints(iterable, Cstr_threshold, Cstr_adjaceny, DotBr_threshold, DotBr_adjacency, unpaired_adjacency)
	
	for line in iterable:

		result = ', '.join(antaRNA_v109.findSequence(line[0],line[1],line[2],**op.to_dict()))
		folding = result.split("\n")[1]
		designed_seq = result.split("\n")[2]
		fasta_id = line[3]

		yield designed_seq, fasta_id


if __name__ == "__main__":
	
	logging.basicConfig(level=logging.INFO)
	logger = logging.getLogger(__name__)
	
	logger.info('Call to RNADesignTools package.')
	
