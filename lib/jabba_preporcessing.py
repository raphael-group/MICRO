#WE ASSUME THAT NEIGHBORING SEGMENTS HAVE iid's THAT DIFFER BY ONE!


import json
import copy
import networkx as nx
import os
from collections import defaultdict

from networkx.algorithms.operators.binary import disjoint_union

from lib import breakpointgraph
from lib import auxilliary
import itertools


import numpy as np

import sys



#        _______________________________
#_______/       CONSTANTS               \_____________________________________________


MIN_NUMBER_OF_CHROMOSOMES = 2


#THESE HAVE BEEN REPLACED:
#extremities of shorter blocks will be joined with an edge if  MIN_BLOCK
MIN_BLOCK_LENGTH = 100
MIN_BLOCK = True
YELLOW = False

#join extremities of a block if their gray degrees are equal to 1 and black degrees are equal to 1
#MIN_BLOCK must be true for FILL to work
FILL = True





#        ___________________________
#_______/       INPUT               \_____________________________________________

def read_json(file_name):
    # Opening JSON file
    f = open(file_name)
    # returns JSON object as a dictionary
    data = json.load(f)
    f.close()
    return data

#read all the graphs and keep those
#that have an annotated mutation
#spanning more than some constant number of chromosomes

def read_type_and_dataset(jabba_graph):
    info = jabba_graph['settings']['description']
    type_start = info.find("type=")  
    dataset_start = info.find("dataset=")

    if type_start == -1:
        type = "NA"  
    else:
        type = info[type_start+5:].split(',')[0]

    if dataset_start == -1:
        dataset = "NA"
    else:
        dataset = info[dataset_start+8:].split(',')[0]

    if dataset == "BACA" and type == "NA":
        type = "PRAD"
    if dataset == "KLUAD" and type =="NA":
        type = "LUAD"

    return type, dataset


def read_annotations(dir_name):

    typeTOtotalnumber  = defaultdict(int)
    typeTOnumberwithchromothripsis = defaultdict(int)


    files_of_interest = [] 


    for json_file in os.listdir(dir_name):
        
        data = read_json(dir_name+json_file)
        type, dataset = read_type_and_dataset(data)


        typeTOtotalnumber[type]+=1

        anotationsTOchr = anotations(data)

        add = False
        for mutation in anotationsTOchr:
            if len(anotationsTOchr[mutation]) >= MIN_NUMBER_OF_CHROMOSOMES:
                add = True
        if add:
            typeTOnumberwithchromothripsis[type] +=1
            files_of_interest.append(json_file)
    
    #plt.bar(*np.unique(values, return_counts=True))
    #plt.xlabel("Number of affected chromosomes")
    #plt.ylabel("Number of chromoplexy")
    #plt.show()
    #plt.savefig("chromoplexy_bar_total")

    return files_of_interest, typeTOtotalnumber, typeTOnumberwithchromothripsis

#Scans the input graph for the mutation annotations.
#For every a mutation, for example a "chromothripsis=1",
#it finds all the chromosomes including block DO225398.jsonextremities 
#incident to edges in this mutation.

#The code outputs a mutation (from the given list) that affects the larger 
#number of chromosomes

def anotations(data):
    iidTOchr = {}
    anotationsTOchr = {}
    
    for interval in data['intervals']:
        iidTOchr[interval['iid']] = interval['chromosome']

    for connection in data['connections']:
        if connection['type'] != 'LOOSE':

            source_ext = end_to_ext(connection['source'])
            source_chr = iidTOchr[source_ext[0]]
            sink_ext = end_to_ext(connection['sink'])
            sink_chr = iidTOchr[sink_ext[0]]

            annotation = connection['annotation']
            #Check if annotation is not empty
            if annotation != '' and annotation != 0:
                #An edge might have multiple annotations separated by dashes
                for mutation in annotation.split("|"):
                    #Select particular mutations
                    if mutation.split("=")[0] == 'chromothripsis':
                        if mutation not in anotationsTOchr:        
                            anotationsTOchr[mutation] = set({source_chr,sink_chr})
                        else: 
                            anotationsTOchr[mutation].add(source_chr)
                            anotationsTOchr[mutation].add(sink_chr)

    return anotationsTOchr



#        ________________________________________________
#_______/       Constructing Breakpoint Graph            \_____________________________________________

#Reads in the JABBA graph and outputs a breakpoint graph
#with possibly some additional edges for short syntenic blocks
#and those that do not have gray edges incident to any of their extremities. 

def jabba_to_breakpoint_graph(data, DUPLICATION_LENGTH, file_name, rearrange, deletion_length):
    
    jabba_type, jabba_dataset = read_type_and_dataset(data)
    bg=nx.MultiGraph(type= jabba_type, dataset = jabba_dataset, filename = file_name)
    iidTOchr = defaultdict(dict)
    iidTOlength = defaultdict(int)

    #reads in the segments
    for interval in data['intervals']:

        chr = interval['chromosome']
        iid = interval['iid']
        iidTOchr[iid] = chr
        cn =  interval['y']
        iidTOlength[iid] = interval['endPoint']-interval['startPoint'] + 1

        #adds tail and head extremity of every syntenic block 
        bg.add_node(((iid,0),chr), loose = False, y = cn)
        bg.add_node(((iid,1),chr), loose = False, y = cn)

    black = 0
    gray = 0 
    loose = 0

    #read in the adjacencies
    for connection in data['connections']:

        source_ext = end_to_ext(connection['source'])
        source_chr = iidTOchr[source_ext[0]]

        if connection['type'] != 'LOOSE':
            sink_ext = end_to_ext(connection['sink'])
            sink_chr = iidTOchr[sink_ext[0]]

            #checks for self loops
            if (source_ext, source_chr) == (sink_ext, sink_chr):
                print('Extremity {} has a self loop'.format((source_ext, source_chr)))

            #adds refferece edge, annotation stores information about an SV type (deletion, chromothripsis...)
            if connection['type'] == 'REF':
                black+=1
                bg.add_edge((source_ext, source_chr),(sink_ext, sink_chr), color ='black', key = connection['annotation'])
            #adds alternative edge, annotation stores information about an SV type (deletion, chromothripanotationsTOchrsis...)
            elif connection['type']  == 'ALT':
                gray +=1
                bg.add_edge((source_ext, source_chr),(sink_ext, sink_chr), color ='gray', key = connection['annotation'])
        else:
            nx.set_node_attributes(bg, {(source_ext, source_chr):{'loose':True}})
            loose+=1
            
    #joins head and tail extremities of every block shorter than a given length. 



    added_edges = 0
    lengths, lengths_strange = [], [] 

    fill(data, bg, DUPLICATION_LENGTH, rearrange, deletion_length)

    #if not rearrange:   
    #    fill_proximity(data, bg, DUPLICATION_LENGTH)
    #else: 
    #    added_edges += fill_duplications(data, bg, DUPLICATION_LENGTH)

    #    added_edges_del, lengths, lengths_strange =  fill_deletions(data, bg, DUPLICATION_LENGTH)
    #    added_edges += added_edges_del

    return bg, iidTOlength, black, gray, added_edges, lengths, lengths_strange


#tail and head extremities are encoded as -iid and +iid
#we transform them into (iid, 0) and (iid, 1)  
def end_to_ext(end):
    if end < 0:
        return (end*(-1),0)
    else:
        return (end, 1)


#black and gray degrees of a vertex in a breakpoint graph
def vertex_degrees(bg, extremity):
    black, gray = 0, 0 
    #scans edges incident to a given vertex and outputs its black and gray degree
    for edge in bg.edges(extremity, 'color'):
        if edge[2] == 'black':
            black+=1
        if edge[2] == 'gray':
            gray+=1
    
    return gray, black

def color_degree(bg, extremity, color):
    degree = 0 
    #scans edges incident to a given vertex and outputs its black and gray degree
    for edge in bg.edges(extremity, 'color'):
        if edge[2] == color:
            degree+=1
    return degree

def fill_duplications(data, bg, constant):
    added_edges = 0
    for i in range(len(data['intervals'])):
        interval = data['intervals'][i]
        iid = interval['iid'] 
        chromosome = interval['chromosome']
        length_interval = get_interval_length(interval)

        tail, tail_black, tail_gray = get_extremity(0, iid, chromosome, bg)
        head, head_black, head_gray = get_extremity(1, iid, chromosome, bg)

        if  tail_gray == 1 and head_gray == 1 and tail_black + head_black == 2:
            
            cn = interval['y']
            
            interval_left = data['intervals'][i-1]
            iid_left = interval_left['iid']
            cn_left = interval_left['y']
            length_left = get_interval_length(interval_left)
            left_head, left_black, left_gray = get_extremity(1, iid_left, chromosome, bg)

            interval_right = data['intervals'][i+1]
            iid_right = interval_right['iid']
            cn_right = interval_right['y']
            length_right = get_interval_length(interval_right)
            right_tail, right_black, right_gray = get_extremity(0, iid_right, chromosome, bg)

            if right_gray + left_gray == 0 and length_interval < constant:# and cn > cn_left and cn > cn_right:
                bg.add_edge(left_head, right_tail, color ='brown', length = length_interval)
                added_edges+=1
 
    return added_edges


def get_extremity(end, iid, chromosome, bg):
    extremity = ((iid,end), chromosome)
    gray_degree, black_degree = vertex_degrees(bg, extremity)
    return extremity, black_degree, gray_degree

def get_interval_length(interval):
    return interval['endPoint']-interval['startPoint'] + 1

def fill_proximity(data, bg, constant):
    for i in range(len(data['intervals'])):
        interval = data['intervals'][i]
        iid = interval['iid'] 
        chromosome = interval['chromosome']
        length_interval = get_interval_length(interval)

        tail, tail_black, tail_gray = get_extremity(0, iid, chromosome, bg)
        head, head_black, head_gray = get_extremity(1, iid, chromosome, bg)

        if tail_black == 1 and head_black == 1:
            interval_left = data['intervals'][i-1]
            iid_left = interval_left['iid']
            left_head, left_black, left_gray = get_extremity(1, iid_left, chromosome, bg)

            interval_right = data['intervals'][i+1]
            iid_right = interval_right['iid']
            right_tail, right_black, right_gray = get_extremity(0, iid_right, chromosome, bg)

            if right_gray + head_gray >=1 and left_gray + tail_gray >=1 and length_interval < constant:
                bg.add_edge(tail, head, color ='yellow', length = length_interval)


def fill(data, bg, constant, rearrangement, deletion_length):
    for i in range(len(data['intervals'])):
        interval = data['intervals'][i]
        iid = interval['iid'] 
        chromosome = interval['chromosome']
        length_interval = get_interval_length(interval)

        tail, tail_black, tail_gray = get_extremity(0, iid, chromosome, bg)
        head, head_black, head_gray = get_extremity(1, iid, chromosome, bg)

        if tail_black == 1 and head_black == 1:
            interval_left = data['intervals'][i-1]
            iid_left = interval_left['iid']
            left_head, left_black, left_gray = get_extremity(1, iid_left, chromosome, bg)

            interval_right = data['intervals'][i+1]
            iid_right = interval_right['iid']
            right_tail, right_black, right_gray = get_extremity(0, iid_right, chromosome, bg)

            if right_gray + left_gray == 0 and tail_gray == 1 and head_gray == 1 and length_interval < constant:
                bg.add_edge(left_head, right_tail, color ='brown', length = length_interval)
            elif right_gray == 1 and left_gray == 1 and tail_gray + head_gray == 0:
                if rearrangement and length_interval < deletion_length:
                    if color_degree(bg, tail, 'brown') == 0 and color_degree(bg, head, 'brown') == 0 and color_degree(bg, left_head, 'yellow') == 0 and color_degree(bg, right_tail, 'yellow') == 0:
                        bg.add_edge(tail, head, color ='red', length = length_interval)
                elif length_interval < constant:
                    bg.add_edge(tail, head, color ='red', length = length_interval)
            elif right_gray + head_gray >=1 and left_gray + tail_gray >=1 and length_interval < constant:
                bg.add_edge(tail, head, color ='yellow', length = length_interval)


def fill_deletions(data, bg, constant):
    added_edges = 0
    lengths = [] 

    lengths_strange = []
    for i in range(len(data['intervals'])):
        interval = data['intervals'][i]
        iid = interval['iid'] 
        chromosome = interval['chromosome']
        length_interval = get_interval_length(interval)

        tail, tail_black, tail_gray = get_extremity(0, iid, chromosome, bg)
        head, head_black, head_gray = get_extremity(1, iid, chromosome, bg)

        if  tail_gray + head_gray == 0 and tail_black + head_black == 2:
            cn = interval['y']
        
            interval_left = data['intervals'][i-1]
            iid_left = interval_left['iid']
            cn_left = interval_left['y']
            length_left = get_interval_length(interval_left)
            left_head, left_black, left_gray = get_extremity(1, iid_left, chromosome, bg)


            interval_right = data['intervals'][i+1]
            iid_right = interval_right['iid']
            cn_right = interval_right['y']
            length_right = get_interval_length(interval_right)
            right_tail, right_black, right_gray = get_extremity(0, iid_right, chromosome, bg)

            if right_gray == 1 and left_gray == 1:
                #if length_interval < constant and (color_degree(bg, tail, 'brown') > 0 or color_degree(bg, head, 'brown') > 0):
                #    bg.add_edge(tail, head, color ='red', length = length_interval)
                #    print('WTF', bg.graph['filename'], iid, chromosome, length_interval)
                if color_degree(bg, tail, 'brown') == 0 and color_degree(bg, head, 'brown') == 0:
                    bg.add_edge(tail, head, color ='red', length = length_interval)
    return added_edges, lengths, lengths_strange


#not all of the iids are used by JaBbA
#but in what follows I suppose that consecutive blocks have consecutive iids 
def fill_block_edges(data, bg):
    
    #Working on detecting small blocks duplicated during chromoplexy
    #Their both extremities have gray degree equal to one
    #neighboring extremities have gray degree equal to zero
    #and their copy number is larger than that of the neighbors
    
    added_edges = 0
    for i in range(len(data['intervals'])):
        interval = data['intervals'][i]
        iid = interval['iid'] 
        chromosome = interval['chromosome']

        tail, tail_black, tail_gray = get_extremity(0, iid, chromosome, bg)
        head, head_black, head_gray = get_extremity(1, iid, chromosome, bg)

        if MIN_BLOCK:
            if interval['endPoint']-interval['startPoint'] < MIN_BLOCK_LENGTH:
                if  tail_gray + head_gray == 0 and tail_black + head_black == 2:
                    bg.add_edge(tail, head, color ='green', length = interval['endPoint']-interval['startPoint'] + 1)
                    added_edges+=1
                elif YELLOW: 
                    bg.add_edge(tail, head, color ='yellow', length = interval['endPoint']-interval['startPoint'] + 1)
            elif FILL:
                #tests if the interval is not at the end of the chromosome
                if  tail_gray + head_gray == 0 and tail_black + head_black == 2:
                    iid = interval['iid']
                    #checks the copy numbers
                    if  1 < iid < len(data['intervals']):
                        cn = interval['y']
                        #the intervals to the left and to the right are data['intervals'][iid-2] and data['intervals'][iid]['iid']
                        cn_left = data['intervals'][i-1]['y']
                        cn_right = data['intervals'][i+1]['y']
                        if cn < cn_left and cn < cn_right:
                            bg.add_edge(tail, head, color ='red', length = interval['endPoint']-interval['startPoint'] + 1)
                            added_edges+=1
    return added_edges



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


#checks if a component has more than one gray edge

def number_of_gray_edges(component):
    gray = 0
    for node in component.nodes():
        for edge in component.edges(node, 'color'):
            if edge[2] == 'gray':
                gray +=1
    return gray/2
    
 

#should be addapted if to be used with yellow edges
def circle(component):
    for node in component.nodes():
        if component.degree(node) != 2:
            return False
    return True        

    
#   for every connected components outputs
#   1) the referece chromosomes to which the vertices belongs.
#   2) annotations of the edges 

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

def extremities_in_components(nxcomponents, indexes):
    nodes = set()
    for index in indexes:
        for node in nxcomponents[index].nodes():
            nodes.add(node[0])
    return nodes


def find_components_with_annotation(component_annotation, annotation):
    components_with_annotation = []
    for i in range(len(component_annotation)):
        if annotation in component_annotation[i]:
            components_with_annotation.append(i)
    return components_with_annotation






#        _________________________________________________________________
#_______/      Connected components of the multi-genome graph            \_____________________________________________

#selects the connected components of the multi-genome graph to which the chromosomes in "subset" belong. 
#component_chr contains a list of chromosomes for every connected component of the breakpoint graph.
#disjoint_unions outputs the chromosomes of the connected components of the multi-genome graph 

def component_of_the_genome_graph(component_chr, subset):
    connected_components_of_subset = set()
    number_of_components = 0
    for disjoint in auxilliary.disjoint_unions(component_chr):
        if subset & disjoint:
            number_of_components += 1
            connected_components_of_subset = connected_components_of_subset.union(disjoint)
    if number_of_components > 1:
        print("Some event spans multiple components of the multi-genome graph")
    return connected_components_of_subset

#selects the components of the multi-genome graph that include annotated events affecting more than MIN_NUMBER_OF_CHROMOSOMES of chromosomes.
#also outputs the maximum number of affected chromosomes by an event in that component. 

def components_with_annotations(anotationsTOchr, component_chr):
    annotations_of_maximums = []
    maximums = []
    for i in anotationsTOchr:
        if len(anotationsTOchr[i]) >= MIN_NUMBER_OF_CHROMOSOMES:
            maximum = component_of_the_genome_graph(component_chr, anotationsTOchr[i])
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


def component_with_high_gray_degree(component):
    nodes = component.nodes()
    for node in nodes:
        gray, black = vertex_degrees(component, node)
        if gray > 1:
            return True
    return False






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










#        __________________________________
#_______/       MISCELIANOUS               \_____________________________________________




#constructs a subgraph of the breakpoint graph only containing 
#those componentes that belong to the connected component "maximum"
#of the multi-genome graph
def subgraph_of_breakpoint_graph(nxcomponents, maximum, component_chr):

    graph=nx.MultiGraph()
    for i in range(len(nxcomponents)):
        if component_chr[i] & maximum:
            graph = nx.union(graph, copy.deepcopy(nxcomponents[i]))
    return graph

def bfb(component, iidTOchr):
    for edge in component.edges(component.nodes(),'color'):
        if edge[2] == 'gray':
            copy_numbers = nx.get_node_attributes(component,'y')
            cn1 = copy_numbers[edge[0]]
            cn2 = copy_numbers[edge[1]]
            if edge[0][0][1]==edge[1][0][1]==1 and edge[1][0][0]-edge[0][0][0] == 1 and cn2 == cn1 - 1 and cn2  == 1:
                next = ((edge[1][0][0]+1,0), edge[1][1])
                if copy_numbers[next] == 0 and iidTOchr[edge[1][0][0]+2]!= edge[1][1]:
                    print(edge)
                    return True
            if edge[0][0][1]==edge[1][0][1]==1 and edge[0][0][0]-edge[1][0][0] == 1 and cn1 == cn2 - 1 and cn1  == 1:
                next = ((edge[0][0][0]+1,0),edge[1][1])
                if copy_numbers[next] == 0 and iidTOchr[edge[0][0][0]+2]!= edge[1][1]:
                    print(edge)
                    return True

def search_for_bfb_in_sample(data,filename):
    bg, black, gray, added_edges, loose, iidTOchr, iidTOlength = jabba_to_breakpoint_graph(data, FILL, MIN_BLOCK)
    nxcomponents = breakpointgraph.breakpointgraphTOcomponents(bg)
    component_chr, annotationTOnodes =  annotation_of_bg_components(nxcomponents)

    components_of_multi_genome_graph = auxilliary.disjoint_unions(component_chr)
    
    for i in range(len(nxcomponents)):
        for subset in components_of_multi_genome_graph:
            if component_chr[i] & subset and len(subset)<23:
                if bfb(nxcomponents[i], iidTOchr):
                    print(filename, subset)


def search_for_bfb(dir_name):
    for json_file in os.listdir(dir_name):
        data = read_json(dir_name+json_file)
        search_for_bfb_in_sample(data,json_file)



class scenario:

    def __init__(self, surrounding, max_rearranged, last_rearranged, rearranged = [], ordered_partition = None, gray = None, proportion_chromo = None, proportion_components = None, deletions = None, duplications = None):
        self.surrounding = surrounding
        self.proportion_chromo = proportion_chromo
        self.proportion_components = proportion_components
        self.max_rearranged = max_rearranged
        self.last_rearranged = last_rearranged
        self.rearranged = rearranged
        self.ordered_partition = ordered_partition
        self.gray = gray
        self.deletions = deletions
        self.duplications = duplications

        self.filename = None
        self.type = None
        self.dataset = None

        self.title = None
        self.annotation = None



class sample:


    def __init__(self, filename, type, dataset):
        self.filename = filename
        self.type = type
        self.dataset = dataset
        self.annotations = []

      
class annotation:
    
    def __init__(self, title):
        self.title = title
        self.scenario_max_sur = None
        self.scenario_removal_max_sur = None


def prox(nodes1, nodes2, constant, iidTOlength):
    proximity = 0
    if nodes1 & nodes2:
        sys.exit("components must be disjoint")
    for node in nodes1:
        if other_extremity(node) in nodes2 and iidTOlength[node[0][0]] < constant:
            proximity+=1
    return proximity


