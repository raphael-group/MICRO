## Overview

MPCR computes the chromosome number of a set of adjacencies introduced by a complex rearrangement. MPCR takes as input a JaBbA graph and outputs the chromosome numbers of the subsets of the novel adjacencies annotated as introduced by complex rearrangements.  

## Installing MPCR

git clone https://github.com/raphael-group/MPCR.git

MPCR is written in python3 and requires the following dependencies: networkx, json, requests.

## Running MPCR

MPCR can be run using the following command:
```
python3 mpcr.py 
```
## Input

Directory 'jabba_data/graphs' holds 22 JaBbA graphs downloaded from http://mskilab.com/gGraph/, the rest of the JaBbA graphs from http://mskilab.com/gGraph/ can be downloaded to the directory 'jabba_data/graphs' by using the following command:  

```
python3 request.py 

```




