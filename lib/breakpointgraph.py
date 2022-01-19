#!/usr/bin/env python3
# Pijus Simonaitis                                               Winter 2018
#
"""
Library for breakpoint graph.

Breakpoint graph is a networkx multigrpah with black and gray edges.
It's connected components are either circles, odd paths or even paths.

1) Telomeric nodes are strings starting with letter 't'


The endponts of the odd paths are merged to form a circle.
Resulting telomer is named "t1"

The even paths are either completed with an edge to form a circle,
or two even paths are merged to form a circle.
Resulting telomers are named "t1" and "t2"  

I distinguish between nxcircles that are netwrokx graphs and 
circles, that are the lists of vertices.

In the latter case I take a convention that for a circle [v_{0},v_{1},\ldots, v_{n}], 
the edge (v_{0},v_{1}) is black. 

"""

import sys
import networkx as nx
import random 
import os
import itertools
import copy

from collections import defaultdict

from . import jabba_preporcessing 

from . import auxilliary

OUT_DIR = "data/breakpointgraphs"

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















#        _____________________________________
#_______/       Even paths to nxcircles        \_____________________________________________

def componentsTOevenpaths(components):
  """
  Given connected components of a breakpoint graph 
  outputs AA and BB paths.

  """
  aapaths, bbpaths = list(), list()
  for component in components: 
    if len(component.nodes())%2 == 0:
       ends = degreeonevertices(component)
       if len(ends) == 2: 
         #take one edge incident to an end of a path 
         [u,v] = list(component.edges(ends[0]))[0]
         if component[u][v][0]['color'] == 'black':
           aapaths.append(component)
         elif component[u][v][0]['color'] == 'gray':
           bbpaths.append(component)
         else:
           sys.exit('Color must be black or gray')
       if len(ends) > 2:
         sys.exit('One of the components is neither a path nor a circle')
  return (aapaths, bbpaths)

def pairofevenpathsTOnxcircles(path1, path2):
  """
  Given a pair of an AA path and a BB paths,
  outputs the two circles obtained after merging these paths.

  """
  ends1, ends2 = degreeonevertices(path1), degreeonevertices(path2)
  
  #compose both paths into a single graph
  #the nodes are uniquely named, thus the new graph has two
  #connected components 
  composed_paths = nx.compose(path1, path2)

  #contract the telomers to obtain the two possible circles 
  c1 = copy.deepcopy(composed_paths)
  c1 = nx.contracted_nodes(c1, ends1[0], ends2[1], self_loops=True)
  c1 = nx.contracted_nodes(c1, ends1[1], ends2[0], self_loops=True)
 
  c2 = copy.deepcopy(composed_paths)
  c2 = nx.contracted_nodes(c2, ends1[0], ends2[0], self_loops=True)
  c2 = nx.contracted_nodes(c2, ends1[1], ends2[1], self_loops=True)

  #rename the telomers
  mapping = {ends1[0]: 't1', ends1[1]: 't2'}
  c1 = nx.relabel_nodes(c1, mapping, copy=False)
  c2 = nx.relabel_nodes(c2, mapping, copy=False)

  return [c1, c2]

def evenpathTOnxcircle(path, col):
  """
  Adds an edge of color col between the ends of an even path,
  renames the ends and outputs the circle.

  """
  ends = degreeonevertices(path)
  
  circle = copy.deepcopy(path)
  circle.add_edge(ends[0],ends[1], color = col)
  mapping = {ends[0]: 't1', ends[1]: 't2'}
  circle = nx.relabel_nodes(circle, mapping, copy=False)

  return circle
       
def evenpathsTOnxcircles(aapaths, bbpaths):
  """
  Given AA paths and BB paths returns 
  a dictionary storing the pairs of circles obtained merging AA and BB paths. 

  'empty_a' and 'empty_b' denotes a path of length 1 joining two telomers.  
  """

  AApathTOBBpathTOnxcircles = defaultdict(dict)
 
  for A in aapaths:
    AApathTOBBpathTOnxcircles[A]['empty_b'] = [evenpathTOnxcircle(A, 'gray')]
  for B in bbpaths:
    AApathTOBBpathTOnxcircles['empty_a'][B] = [evenpathTOnxcircle(B, 'black')]
  
  for A in aapaths:
    for B in bbpaths:
      AApathTOBBpathTOnxcircles[A][B] = pairofevenpathsTOnxcircles(A, B) 
  
  return AApathTOBBpathTOnxcircles 
   

def plotevenpathnxcircles(AApathTOBBpathTOnxcircles,name):
  nxcircles = list()
  for A in AApathTOBBpathTOcircles:
    for B in AApathTOBBpathTOcircles[A]:
      nxcircles += AApathTOBBpathTOcircles[A][B]
  plotgraphs(nxcircles, name+'_even_paths')

#        ___________________________________________________
#_______/       odd paths and nxcircles to nxcircles        

def breakpointgraphTOnxcircles(bg):
  """
  Transforms a breakpoint graph into a list of circles that are still networx graphs.
  The ends of an odd path are contracted to form a circle.

  """
  nxcircles = list()
  for i in nx.connected_components(bg):
    H = nx.MultiGraph(bg.subgraph(i))
    #odd paths having more than 2 edges
    if len(H.nodes())%2 == 1 and len(H.nodes())>3:
      #Transform an odd path into a circles
      #by contracting its ends and renaming them by 't'
      ends = list()
      for u in H.nodes(): 
        if H.degree(u) == 1:
          ends.append(u)
      H = nx.contracted_nodes(H, ends[0], ends[1], self_loops=True)
      mapping = {ends[0]: 't1'}
      H = nx.relabel_nodes(H, mapping, copy=False)
      nxcircles.append(H)
    #circles having more than 2 edges
    elif len(degreeonevertices(H)) == 0 and len(H.nodes)>2:
      nxcircles.append(H) 
  return nxcircles

#        __________________________________
#_______/       nxgraphs to circles        \_____________________________________________


def nxcircleTOcircle(nxcircle):
  """
  Returns a list of nodes respecting their order on the circle and the first edge being black.

  """

  #finds a gray edge 
  edge = None
  for i in nxcircle.edges():
    if (nxcircle[i[0]][i[1]][0]['color'] == 'gray' and edge == None):
      edge = i

  #removes a chosen gray edge
  nxcircle.remove_edge(edge[0], edge[1], 0)

  #finds a path joining the ends of the gray edge
  circle = nx.shortest_path(nxcircle, source=edge[0], target=edge[1], weight=None)

  return circle

def genomesTOcircles(genomesfile, show = False):
  #loads a pair of genomes
  (g1,g2) = pickle.load(open(genomesfile,"rb"))
  #constructs a breakpoint graph
  bg = genomesTObreakpoinggraph(g1, g2)
  #extracts the connected components of bg
  nxcircles = breakpointgraphTOnxcircles(bg)
  #plots the circles in OUT_GRAPHS
  if show:
    plotgraphs(nxcircles, g1.name+'-'+g2.name)
    print("Finished plotting the graphs.")

  #transforms networkx graphs of circles into the lists of vertices

  circles = []
  for nxcircle in nxcircles:
    circle = nxcircleTOcircle(nxcircle)
    if circle != False:
      circles.append(circle)
      
  return circles 






#        ___________________________________________________________
#_______/       BLOCKS        \_______________________

def transformbreakpointgraph_blocks(bg,components):
  """
  Takes a list of edge-disjoint Eulerian subgraphs of bg, 
  and sort these subgraphs in bg. 


  1) as implemented this only works graphs without parallel edges,
  and the componented are supposed to be non-trivial.
  This is because remove_edge(u,v,color = 'c') does not seem to work 
  and I would need to implement it myself, or use 'key' in the edges. 
  
  2) Yellow edges keep the components together, but they are not adjacencies

  """
  transformed_bg = copy.deepcopy(bg)

  for component in components:
    for edge in component.edges.data():
      if edge[2]['color'] == 'black':
        transformed_bg.remove_edge(edge[0], edge[1])
      elif edge[2]['color'] != 'yellow':
        transformed_bg.add_edge(edge[0], edge[1], color='black')

  return transformed_bg 



def run_a_scenario_blocks(bg, ordered_components):
  """
  A list of lists of edge-disjoint eulerian subgraphs of a breakpoint graph corresponds to a multi-break 
  scenario without breakpoint reuse.  

  We run this scenario and output its statistics. 
  
  A breakpoint graph that is being transformed is an element in a list of length 1.
  I pop it out, transform it, and append it back. 
  """

  breakpoint_graphs = [bg] 
  components_of_ggs = [breakpointgraphTOcomponents(breakpointgraphTOgenomegraph(bg))]


  #breakpointgraphTOcomponents(genome_graph)


  #tests if no circular chromosomes are formed during the scenario. 
  good = True
  #a list of (the number of affected chromosomes, the number of resulting chromosomes). 
  scenario_stats = []
  
  for components in ordered_components:
    #performs a multibreak
    breakpoint_graphs.append(transformbreakpointgraph(breakpoint_graphs.pop(),components))
    components_of_ggs.append(breakpointgraphTOcomponents(breakpointgraphTOgenomegraph(breakpoint_graphs[0])))

    #collects statistics on the number of resulting chromosomes
    starting = len(components_of_ggs[0])
    resulting = len(components_of_ggs[1])
    if resulting != starting:
      good = False

    #collects statistics on the number of affected chromosomes
    affected = len(number_of_affected_chromosomes(components_of_ggs.pop(0), components))
    scenario_stats.append([affected, resulting])

  return [good, scenario_stats]






def run_a_partition(ordered_partition,nxcomponents, bg, min_affected, last_affected):
  scenario = auxilliary.ordered_partition_of_components(ordered_partition,nxcomponents)
  [good, scenario_stats] = run_a_scenario_blocks(bg, scenario)
  #print(scenario_stats)

  minimum = 0
  for [affected, resulting] in scenario_stats:
    if affected > minimum:
      minimum = affected
  if scenario_stats[-1][0] < last_affected:
    last_affected = scenario_stats[-1][0]
  if minimum <=  min_affected:
    min_affected = minimum
    #print(min_affected, scenario_stats[-1][0], ordered_partition)


  return min_affected, last_affected 
 




def iterate_through_scenarios(bg, components, singletons, max_index):
  
  nxcomponents = breakpointgraphTOcomponents(bg)
  genome_graph = breakpointgraphTOgenomegraph(bg)
  gg_components = breakpointgraphTOcomponents(genome_graph)

  min_affected = len(gg_components)
  last_affected = len(gg_components)

  

  if singletons:
    for ordered_partition in auxilliary.order_without_partition(components): 
      ordered_partition = list(ordered_partition)
      ordered_partition.append([max_index])
      min_affected = run_a_partition(ordered_partition,nxcomponents, bg, min_affected)
  else: 
    for ordered_partition in auxilliary.order(components):
      ordered_partition = list(ordered_partition)
      ordered_partition.append([max_index])
      min_affected, last_affected = run_a_partition(ordered_partition,nxcomponents, bg, min_affected, last_affected)
   

  #print(good_scenarios, good_scenarios/number_of_scenarios)
  #print("min_affected", min_affected)
  return min_affected, last_affected 


def breakpointgraphTOgenomegraph(bg):
  """
  Created a genome graph, which is a 2-edge-colored multigraph
  with black edges representing the adjacencies in g1 (thus equal to the black edges of the breakpoint graph),
  while the red edges join the 0 and 1 extremities of every syntenic block.  

  Relies on the fact that a telomeric node is a string starting with letter 't'
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
    #check if the node is not a telomere
    if node[0] != 't':
      blocks.add(node[0])


  #add block edges
  for block in blocks:
    genome_graph.add_edge((block,0),(block,1), color = 'yellow')
    

  return genome_graph

def genomesTObreakpoinggraph(g1, g2):
  
  """
  Generates a breakpoint graph of genomes g1 and g2 with a separate vertex for every telomere.
  
  This means that every even path, odd path and circle of a breakpoint graph is a separate
  connected component.

  Vertices of a breakpoint graph have format (ext, c1, c2), where ci is the chromosome
  containing ext in gi.
  
  ext is of a format (block ID, extremity (either 0 or 1))

  """
 
  bg=nx.MultiGraph()
  
  #Black edges are the adjacencies of g1
  for c1 in g1.chrTOext:
    #Telomeric adjacencies
    chrom1 = g1.chrTOext[c1]
    ext = chrom1[0]
    c2 = g2.extTOchrpos[ext][0]
    bg.add_edge('t{}1b'.format(c1),(ext, c1, c2), color='black')

    ext = chrom1[-1]
    c2 = g2.extTOchrpos[ext][0]
    bg.add_edge('t{}2b'.format(c1),(ext, c1, c2), color='black')
   
    #Internal adjacencies 
    for i in range(1,len(chrom1)-1,2):
      ext1 = chrom1[i]
      ext2 = chrom1[i+1]
      c21 = g2.extTOchrpos[ext1][0]
      c22 = g2.extTOchrpos[ext2][0]
      bg.add_edge((ext1, c1, c21),(ext2, c1, c22), color='black')    

  #Gray edges are the adjacencies of g2
  for c2 in g2.chrTOext:
    #Telomeric adjacencies
    chrom2 = g2.chrTOext[c2]
    ext = chrom2[0]
    c1 = g1.extTOchrpos[ext][0]
    bg.add_edge('t{}1g'.format(c2),(ext, c1, c2), color='gray')

    ext = chrom2[-1]
    c1 = g1.extTOchrpos[ext][0]
    bg.add_edge('t{}2g'.format(c2),(ext, c1, c2), color='gray')

    #Internal adjacencies 
    for i in range(1,len(chrom2)-1,2):
      ext1 = chrom2[i]
      ext2 = chrom2[i+1]
      c11 = g1.extTOchrpos[ext1][0]
      c12 = g1.extTOchrpos[ext2][0]
      bg.add_edge((ext1, c11, c2),(ext2, c12, c2), color='gray')

  return bg



def NEWgenomesTObreakpoinggraph(g1, g2, nodes):
  """
  a variant of genomesTObreakpoinggraph that orders the nodes in a graph 
  according to their physical location and ommits chromosome names from the names.

  Generates a breakpoint graph of genomes g1 and g2 with a separate vertex for every telomere.
  
  This means that every even path, odd path and circle of a breakpoint graph is a separate
  connected component.

  Vertices of a breakpoint graph have format (block ID, extremity (either 0 or 1)).

  """
 
  bg=nx.MultiGraph()

  for i in nodes:
    bg.add_node(i)

  
  #Black edges are the adjacencies of g1
  for c1 in g1.chrTOext:
    #Telomeric adjacencies
    chrom1 = g1.chrTOext[c1]
    ext = chrom1[0]
    c2 = g2.extTOchrpos[ext][0]
    bg.add_edge('t{}1b'.format(c1),ext, color='black')

    ext = chrom1[-1]
    c2 = g2.extTOchrpos[ext][0]
    bg.add_edge('t{}2b'.format(c1),ext, color='black')
   
    #Internal adjacencies 
    for i in range(1,len(chrom1)-1,2):
      ext1 = chrom1[i]
      ext2 = chrom1[i+1]
      c21 = g2.extTOchrpos[ext1][0]
      c22 = g2.extTOchrpos[ext2][0]
      bg.add_edge(ext1,ext2, color='black')    

  #Gray edges are the adjacencies of g2
  for c2 in g2.chrTOext:
    #Telomeric adjacencies
    chrom2 = g2.chrTOext[c2]
    ext = chrom2[0]
    c1 = g1.extTOchrpos[ext][0]
    bg.add_edge('t{}1g'.format(c2),ext, color='gray')

    ext = chrom2[-1]
    c1 = g1.extTOchrpos[ext][0]
    bg.add_edge('t{}2g'.format(c2),ext, color='gray')

    #Internal adjacencies 
    for i in range(1,len(chrom2)-1,2):
      ext1 = chrom2[i]
      ext2 = chrom2[i+1]
      c11 = g1.extTOchrpos[ext1][0]
      c12 = g1.extTOchrpos[ext2][0]
      bg.add_edge(ext1,ext2, color='gray')

  return bg
