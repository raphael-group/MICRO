#!/usr/bin/env python3
# Pijus Simonaitis                                               Winter 2018
#

import networkx as nx
import os
from collections import defaultdict


from . import auxilliary


OUT_DIR = "data/breakpointgraphs"


#        ___________________________
#_______/       Construction        \_____________________________________________


# Vetex format is ((block_id, strand), chromosome)
def breakpointgraphTOgenomegraph(bg):
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


def graphTOcomponents(bg):
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


 
def graphTOcircles(bg):
  """
  
  Transforms a graph into a list of connected components.
  These connecteted components are just the "views" of the graph,
  they are nots supposed to be touched, 
  nor used after the corresponding edges in the graph were touched. 

  """
  circles = list()
  for i in nx.connected_components(bg):
    H = nx.MultiGraph(bg.subgraph(i))
    if circle(H):
        circles.append(H)
  return circles 



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
  

  """
  for component in components:
    for edge in component.edges.data():
      if edge[2]['color'] == 'black':
        bg.remove_edge(edge[0], edge[1])
      if edge[2]['color'] == 'gray':
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


def run_a_partition(ordered_partition,nxcomponents, bg, templated_insertions):
  scenario = auxilliary.ordered_partition_of_components(ordered_partition,nxcomponents)
  rearranged = []
  genome_graph = breakpointgraphTOgenomegraph(bg)

  #I do not need to generate the last breakpoint graph 
  for components in scenario:
    nodes = set()
    for component in components:
      nodes = nodes.union(component.nodes())

    rearranged.append(number_of_affected_chromosomes_by_nodes(genome_graph, nodes))
    if templated_insertions:
      transformbreakpointgraph(genome_graph,components)
    else:
      transformbreakpointgraph_duplications(genome_graph,components)
  return rearranged, genome_graph

  

#        _________________________________________________________________
#_______/      Connected components of the multi-genome graph            \_____________________________________________

# Two chromosomes are 'related' if they there is a connected component of the breakpoint graph that includes vertices from both of these chromsoomes. 
# We construct a graph M whose vertices are chromosomes and they are joined by an edge iff they are 'related'
# disjoint_unions(component_chr) outputs the vertices of the connected components of this graph M.
# and related_chromosomes(component_chr, subset) outputs the vertices of a connected component of M
# that includes the chromosomes from subset.

def related_chromosomes(component_chr, subset):
    connected_components_of_subset = set()
    number_of_components = 0
    for disjoint in auxilliary.disjoint_unions(component_chr):
        if subset & disjoint:
            number_of_components += 1
            connected_components_of_subset = connected_components_of_subset.union(disjoint)
    if number_of_components > 1:
        print("Some event spans multiple components of the multi-genome graph")
    return connected_components_of_subset

#selects the components of the graph M that include 
#also outputs the maximum number of affected chromosomes by an event in that component. 

def subsets_of_related_chromosomes(anotationsTOchr, component_chr, min_number_of_chromosomes):
    annotations_of_maximums = []
    maximums = []
    for i in anotationsTOchr:
        if len(anotationsTOchr[i]) >= min_number_of_chromosomes:
            maximum = related_chromosomes(component_chr, anotationsTOchr[i])
            add = True
            for j in range(len(maximums)):
                set = maximums[j]
                if set & maximum:
                    add = False
                    annotations_of_maximums[j].append(i)
            if add:
                maximums.append(maximum)
                annotations_of_maximums.append([i])
    return maximums, annotations_of_maximums 


def annotation_of_bg_components(nxcomponents):
    
    component_chr = []
    annotationTOnodes = defaultdict(dict)
    annotationTOcomponents = defaultdict(set)


    for index in range(len(nxcomponents)):
        graph = nxcomponents[index]
        chromosomes = set()
        for i in graph.edges(keys = True):
            if i[2] != '' and i[2] != 0:
                for mutation in str(i[2]).split("|"):
                    if mutation not in annotationTOnodes:
                        annotationTOnodes[mutation] = set()
                    annotationTOnodes[mutation].add(i[0])
                    annotationTOnodes[mutation].add(i[1])
                    annotationTOcomponents[mutation].add(index)

        for node in graph.nodes():
            chromosomes.add(node[1])

        component_chr.append(chromosomes)

    return component_chr, annotationTOnodes, annotationTOcomponents




def find_components_with_annotation(component_annotation, annotation):
    components_with_annotation = []
    for i in range(len(component_annotation)):
        if annotation in component_annotation[i]:
            components_with_annotation.append(i)
    return components_with_annotation





#        ______________________________
#_______/       Study components       \_____________________________________________



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



def number_of_gray_edges(component):
    gray = 0
    for node in component.nodes():
        for edge in component.edges(node, 'color'):
            if edge[2] == 'gray':
                gray +=1
    return gray/2


    
def circle(component):
    for node in component.nodes():
        if component.degree(node) != 2:
            return False
    return True    


def extremities_in_components(nxcomponents, indexes):
    nodes = set()
    for index in indexes:
        for node in nxcomponents[index].nodes():
            nodes.add(node[0])
    return nodes


#        __________________________________________________________________________
#_______/       Analyze the connected components of the breakpoint graph            \_____________________________________________



#checks if the component "resembles" a deletion or a duplication of a single segment. 
#for now this means that a component has 4 nodes that correspond to 4 consecutive block extremities,
#i would need to also check the edges and copy numbers to be sure that this corresponds to a deletion/duplication

def deletion_or_duplication(component):
    if len(component.nodes()) == 4:
        nodes = []
        chromosomes = set()
        for node in component.nodes():
            nodes.append(node[0])
            chromosomes.add(node[1])
        if len(chromosomes) == 1:
            nodes.sort()
            if nodes[0][0] + 1 == nodes[1][0] == nodes[2][0]  == nodes[3][0] - 1 and nodes[0][1] == 1 and nodes[3][1] == 0:
                return True

    return False

# Compute the separation of a component from others
def other_extremity(node):
    if node[0][1] == 0:
        return((node[0][0],1),node[1])
    else:
        return((node[0][0],0),node[1])

def surrounding_blocks(component, iidTOlength):
    nodes = set(component.nodes())
    min_length = float('inf')
    for node in nodes:
        if other_extremity(node) not in nodes:
            min_length = min(iidTOlength[node[0][0]], min_length)
    return min_length

def surrounding_blocks_components(nxcomponents, indexes, iidTOlength):
    nodes = set()
    for index in indexes:
        nodes = nodes.union(set(nxcomponents[index].nodes()))
    min_length = float('inf')
    for node in nodes:        
        if other_extremity(node) not in nodes:
            min_length = min(iidTOlength[node[0][0]], min_length)
    return min_length




#        _________________________
#_______/       Auxilliary        \_____________________________________________



#        _____________________________________________
#_______/      Modify the breakpoint graph            \_____________________________________________



#joins the open ends of the components together
#to obtain an Eulerian graph
def join_paths_to_circles(good_ones):
    #loose_ends store -1 if the component is already a circle
    #the number vertices of degree one that are "loose" in JabbA otherwise.
    loose_ends = []
    for component in good_ones:
        if len(component.nodes()) > 2:
            degree_one = []
            loose = 0 
            for node in component.nodes('loose'):
                degree = 0
                for edge in component.edges(node, 'color'):
                    if edge[2] != 'yellow':
                        degree+=1
                if degree%2 == 1:
                    degree_one.append(node)
                    if node[1] == True:
                        loose+=1

            if len(degree_one) == 0:
                loose_ends.append(-1)
            else: 
                loose_ends.append(loose)
            j = 0
            while j < len(degree_one)-1:
                component.add_edge(degree_one[j][0],degree_one[j+1][0], color = 'blue')
                j = j+2
            
    return loose_ends

# Colors the edges annotated with an annotation in blue
def recolor(graph, annotation):
    for edge in graph.edges(keys=True):
        if edge[2] != '' and edge[2] != 0:
            for mutation in str(edge[2]).split("|"):
                if mutation == annotation:
                    node1 = edge[0]
                    node2 = edge[1]
                    key = edge[2]
                    nx.set_edge_attributes(graph, {(node1, node2, key): {'color': 'blue'}})


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

