#!/usr/bin/env python3
# Pijus Simonaitis                                               Winter 2018
#

import networkx as nx
import os
from . import auxilliary

OUT_DIR = "data/breakpointgraphs"


def length_of_deletions(indexes, components):
    length = 0
    for i in indexes:
        for edge in components[i].edges.data():
            if (edge[2]['color'] == 'red'):
                length+=edge[2]['length']
    return length

def length_of_duplications(indexes, components):
    length = 0
    for i in indexes:
        for edge in components[i].edges.data():
            if (edge[2]['color'] == 'brown'):
                length+=edge[2]['length']
    return length


#        ___________________________
#_______/       Construction        \_____________________________________________


# Vetex format is ((block_id, strand), chromosome)
def breakpointgraphTOgenomegraphJABBA(bg):
  """
  Created a genome graph, which is a 2-edge-colored multigraph
  with black edges representing the adjacencies in g1 (thus equal to the black edges of the breakpoint graph),
  while the pink edges join the 0 and 1 extremities of every syntenic block.  
  """
  genome_graph = nx.MultiGraph()

  #copy the black edges
  for edge in bg.edges.data():
    if (edge[2]['color'] == 'black'):
      genome_graph.add_edge(edge[0],edge[1], color = 'black')

  #obtain the set of the blocks 
  nodes = bg.nodes()
  blocks = set()
  for node in nodes:
      blocks.add((node[0][0],node[1]))


  #add block edges
  for block in blocks:
    genome_graph.add_edge(((block[0],0),block[1]),((block[0],1),block[1]), color = 'pink')

  return genome_graph


def breakpointgraphTOcomponents(bg):
  """
  
  Transforms a graph into a list of connected components.
  These connecteted components are just the "views" of the graph,
  they are nots supposed to be touched, 
  nor used after the corresponding edges in the graph were touched. 

  """
  nxcomponents = list()
  for i in nx.connected_components(bg):
    H = nx.MultiGraph(bg.subgraph(i))
    nxcomponents.append(H)
  return nxcomponents



#        ___________________________________________________________
#_______/       Multi-break scearios without breakpoint reuse        \_______________________

def transformbreakpointgraph(bg,components):
  """
  Takes a list of edge-disjoint Eulerian subgraphs of bg, 
  and sort these subgraphs in bg. 


  1) as implemented this only works graphs without parallel edges,
  and the componented are supposed to be non-trivial.
  This is because remove_edge(u,v,color = 'c') does not seem to work 
  and I would need to implement it myself, or use 'key' in the edges. 
  
  2) Yellow edges keep the components together, but they are not adjacencies

  """
  for component in components:
    for edge in component.edges.data():
      if edge[2]['color'] == 'black':
        bg.remove_edge(edge[0], edge[1])
      elif edge[2]['color'] != 'yellow':
        bg.add_edge(edge[0], edge[1], color='black')


def replace_extremity(duplicated_iid, extremity):
  end = extremity[0][1]
  iid = extremity[0][0]
  chromosome = extremity[1]
  if iid in duplicated_iid:
    if end == 0:
      return ((iid-1,1),chromosome)
    else:
      return ((iid+1,0),chromosome)
  return extremity


#in this case we will pretend that the duplicated segment is absent 
#and the chromosomes are joined as during a translocation.
#this will require renaming the extremities that belong to a duplicated segment 
#which is done by the replace_extremity

def transformbreakpointgraph_duplications(bg,components):
  for component in components:
    duplicated_iid = set()
    for edge in component.edges.data():
      if edge[2]['color'] == 'black': 
        bg.remove_edge(edge[0], edge[1])
      elif edge[2]['color'] == 'brown':
        duplicated_iid.add(min(edge[0][0][0], edge[1][0][0])+1)
    for edge in component.edges.data():
      if edge[2]['color'] == 'gray':
        extremity1 = replace_extremity(duplicated_iid, edge[0])
        extremity2 = replace_extremity(duplicated_iid, edge[1])
        bg.add_edge(extremity1, extremity2, color='black')
    


def number_of_affected_chromosomes_by_nodes(genome_graph, nodes):
  affected = 0
  nodes_of_gg = set()

  for node in nodes:
      if node not in nodes_of_gg:
        nodes_of_gg = nodes_of_gg.union(nx.node_connected_component(genome_graph, node))
        affected+=1
  return affected


def run_a_partition_jabba(ordered_partition,nxcomponents, bg, TEMPLATED_INSERTIONS):
  scenario = auxilliary.ordered_partition_of_components(ordered_partition,nxcomponents)
  rearranged = []
  genome_graph = breakpointgraphTOgenomegraphJABBA(bg)

  #I do not need to generate the last breakpoint graph 
  for components in scenario:
    nodes = set()
    for component in components:
      nodes = nodes.union(component.nodes())

    rearranged.append(number_of_affected_chromosomes_by_nodes(genome_graph, nodes))
    if TEMPLATED_INSERTIONS:
      transformbreakpointgraph(genome_graph,components)
    else:
      transformbreakpointgraph_duplications(genome_graph,components)
  return rearranged, genome_graph

  

#        _________________________
#_______/       Auxilliary        \_____________________________________________


def degreeonevertices(G):
  """
  Outputs a list of degree 1 vertices of G.
  """
  ends = []
  for u in G.nodes(): 
    if G.degree(u) == 1:
      ends.append(u)
  return ends 

def plotgraphs(graphs, name):
  """
  Plots graphs from a given list in OUT_DIR

  """

  if not os.path.exists(OUT_DIR):
    os.makedirs(OUT_DIR)  

  index = 0
  for G in graphs:
    graph = nx.nx_pydot.to_pydot(G)
    graph.write_png('{}/{}-{}.png'.format(OUT_DIR,name,str(index)))
    index +=1

def plotgraphs_with_indexes(graphs, name, indexes):
  """
  Plots graphs from a given list in OUT_DIR

  """

  if not os.path.exists(OUT_DIR):
    os.makedirs(OUT_DIR)  

  for i in range(len(graphs)):
    graph = nx.nx_pydot.to_pydot(graphs[i])
    graph.write_png('{}/{}-{}.png'.format(OUT_DIR,name,str(indexes[i])))

