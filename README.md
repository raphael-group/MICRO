## Overview

Minimum Chromosome Number (MICRO) computes the chromosome number C of a subset of the novel adjacencies introduced by a complex rearrangement, and outputs C if it is lower than the number of chromosomes affected by that complex rearrangement in the reference genome. MPCR takes as input a [JaBbA graph](https://www.sciencedirect.com/science/article/pii/S0092867420309971?via%3Dihub).

## Installing MICRO

```
git clone https://github.com/raphael-group/MICRO.git
```
MICRO is written in python3 and requires the following dependencies: networkx, json, requests.

## Running MICRO

Running
```
python3 micro.py 
```
applies MICRO to the JaBbA graphs that are present in the directory 'jabba_data/graphs'. Upon instalation 'jabba_data/graphs' holds 22 JaBbA graphs from the [gGnome.js portal](http://mskilab.com/gGraph/) discussed in the paper.

## Data 

The rest of the JaBbA graphs can be downloaded from [gGnome.js portal](http://mskilab.com/gGraph/). These graphs can be downloaded to the directory 'jabba_data/graphs' by running  

python3 request.py 






