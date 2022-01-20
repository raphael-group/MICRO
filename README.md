## Overview

MPCR computes the chromosome number c of a subset of the novel adjacencies introduced by a complex rearrangement and outputs c if c is lower than the number of chromosomes affected by that complex rearrangement in the reference genome. MPCR takes as input a [JaBbA graph](https://www.sciencedirect.com/science/article/pii/S0092867420309971?via%3Dihub).

## Installing MPCR

```
git clone https://github.com/raphael-group/MPCR.git
```

MPCR is written in python3 and requires the following dependencies: networkx, json, requests.

## Running MPCR

Running
```
python3 mpcr.py 
```
applies MPCR to the JaBbA graphs that are present in the directory 'jabba_data/graphs'. Upon instalation 'jabba_data/graphs' holds 22 JaBbA graphs from the [gGnome.js portal](http://mskilab.com/gGraph/) for which MPCR identifies a subset b of the novel adjacencies introduced by a chromothripsis event with the chromosome number of b being lower than the number of chromosomes that b affects in the reference genome. 

The rest of the JaBbA graphs from [gGnome.js portal](http://mskilab.com/gGraph/) can be downloaded to the directory 'jabba_data/graphs' by running  

```
python3 request.py 

```





