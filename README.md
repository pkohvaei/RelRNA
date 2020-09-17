# RelRNA

RelRNA is a framework for solving the problem of "RNA inverse folding", as used in my [Master's thesis](http://www.bioinf.uni-freiburg.de/Lehre/Theses/MA_Parastou_Kohvaei.pdf). It consists of two subsystems: one responsible for decomposing secondary structures of sample RNAs into structural features and building a structural features corpus.  It also extracts neighborhood connectivity models of structural features in the form of N-grams. The other subsystem is a reinforcement learning framework  which  uses  the  corpus  and  connectivity  rules  to  produce  models for generating structures which are similar to the samples.

**This code is merely for demonstration purposes.**
Programming language: Python 2.7

![alt text](./data/sample_structure.png)

### Summary ###

A non-coding RNA molecule functionality depends on its structure, which in turn, is determined by the specific arrangement of its nucleotides.  The inverse folding of an RNA refers to the problem of designing an RNA sequence which will fold into a desired structure.  This work introduces a basic system that given a set of sample RNA secondary structures, produces models which generate structures similar to the sample set.  The objectives and constraints are automatically extracted from samples. For doing this, RelRNA system is designed which generates models by performing learning on families of RNA sequences.
