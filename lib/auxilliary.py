#!/usr/bin/env python3
# Pijus Simonaitis                                               Spring 2021
#

"""
Useful auxilliary functions 
"""

import itertools


def sign(i):
  if i < 0:
    return 't'
  return 'h'

def mod(i):
  if i < 0:
    return -i
  return i

#        ________________________________
#_______/       Ordered Partitions       \_____________________________________________



def partition(collection):
  """
  Iterates throug all the partitions of a list.

  I have copied the code from internet and would need to rewrite it myself.  
  
  """

  if len(collection) == 1:
    yield [ collection ]
    return

  first = collection[0]
  for smaller in partition(collection[1:]):
    # insert `first` in each of the subpartition's subsets
    for n, subset in enumerate(smaller):
        yield smaller[:n] + [[ first ] + subset]  + smaller[n+1:]
    # put `first` in its own subset 
    yield [ [ first ] ] + smaller

def order(collection):
  """
  Iterates through all permutations of a list.

  itertools.permutations(unordered_partition) yields (a,) for input ['a'],
  so I just avoid this comma here.

  """
  
  for unordered_partition in partition(collection):
    if len(unordered_partition) != 1:
      for ordered_partition in itertools.permutations(unordered_partition):
        yield ordered_partition
        

def ordered_partition_of_components(ordered_partition,components):
  """
  Given an ordered partition of indexes outputs an ordered partition of elements in a list
  
  """

  scenario = []
  for subset in ordered_partition:
    operation = []
    for element in subset:
      operation.append(components[element])
    scenario.append(operation)
  return scenario


def order_without_partition(components):
  """
  Outputs ordered paritions where every subset is a singleton.

  """

  if len(components) == 1:
    yield [components]
  else:
    for ordering in itertools.permutations(components):
      partition = []
      for element in ordering:
        partition.append([element])
      yield partition
        


#        ___________________________
#_______/      Set Operations       \_____________________________________________




#input: a list of sets
def disjoint_unions(sets):
  
  #output: parition of the elements in the input sets 
  disjoint = []
  #for every subset in sets we check if it
  #intersects any of the sets in disjoint
  for subset in sets:
    remove = []
    #we continue until we fin a non empty intersection
    for disjoint_set in disjoint:
    #if intersection is not empty, then we replace disjoint with 
    #a union of two sets
      if not subset.isdisjoint(disjoint_set):
        subset = subset.union(disjoint_set)
        remove.append(disjoint_set)

    for i in remove:
      disjoint.remove(i)
    disjoint.append(subset)
    
  return disjoint 