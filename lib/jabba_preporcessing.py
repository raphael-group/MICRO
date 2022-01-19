#WE ASSUME THAT NEIGHBORING SEGMENTS HAVE iid's THAT DIFFER BY ONE!

import json
import copy
import networkx as nx
import os
from collections import defaultdict

from lib import auxilliary

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


def read_annotations(dir_name, min_number_of_chromosomes):
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
            if len(anotationsTOchr[mutation]) >= min_number_of_chromosomes:
                add = True
        if add:
            typeTOnumberwithchromothripsis[type] +=1
            files_of_interest.append(json_file)

    return files_of_interest, typeTOtotalnumber, typeTOnumberwithchromothripsis

#Scans the input graph for the mutation annotations.
#For every a mutation, for example a "chromothripsis=1",
#it finds all the chromosomes including block extremities 
#incident to edges in this mutation.

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

def jabba_to_breakpoint_graph(data, duplication_length, file_name, deletion_length):
    
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
            
    #identifies blocks deleted and duplicated by multi-breaks and augment the graph with the corresponding edges 
    fill(data, bg, duplication_length, deletion_length)

    return bg, iidTOlength


def fill(data, bg, duplication_length, deletion_length):
    for i in range(len(data['intervals'])):
        interval = data['intervals'][i]
        iid = interval['iid'] 
        chromosome = interval['chromosome']
        length_interval = get_interval_length(interval)

        tail, tail_black, tail_gray = get_extremity(0, iid, chromosome, bg)
        head, head_black, head_gray = get_extremity(1, iid, chromosome, bg)

        #makes sure that not an interval at an end of a chromosome
        if tail_black == 1 and head_black == 1:
            interval_left = data['intervals'][i-1]
            iid_left = interval_left['iid']
            left_head, left_black, left_gray = get_extremity(1, iid_left, chromosome, bg)

            interval_right = data['intervals'][i+1]
            iid_right = interval_right['iid']
            right_tail, right_black, right_gray = get_extremity(0, iid_right, chromosome, bg)

            if right_gray + left_gray == 0 and tail_gray == 1 and head_gray == 1 and length_interval < duplication_length:
                bg.add_edge(left_head, right_tail, color ='brown', length = length_interval)
            elif right_gray == 1 and left_gray == 1 and tail_gray + head_gray == 0:
                if length_interval < deletion_length:
                    if color_degree(bg, tail, 'brown') == 0 and color_degree(bg, head, 'brown') == 0:
                        bg.add_edge(tail, head, color ='red', length = length_interval)


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

def get_extremity(end, iid, chromosome, bg):
    extremity = ((iid,end), chromosome)
    gray_degree, black_degree = vertex_degrees(bg, extremity)
    return extremity, black_degree, gray_degree

def get_interval_length(interval):
    return interval['endPoint']-interval['startPoint'] + 1



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

def components_with_annotations(anotationsTOchr, component_chr,min_number_of_chromosomes):
    annotations_of_maximums = []
    maximums = []
    for i in anotationsTOchr:
        if len(anotationsTOchr[i]) >= min_number_of_chromosomes:
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





#        __________________________________________________________________________
#_______/       Analyze the connected components of the breakpoint graph            \_____________________________________________


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


def find_components_with_annotation(component_annotation, annotation):
    components_with_annotation = []
    for i in range(len(component_annotation)):
        if annotation in component_annotation[i]:
            components_with_annotation.append(i)
    return components_with_annotation




def component_with_high_gray_degree(component):
    nodes = component.nodes()
    for node in nodes:
        gray, black = vertex_degrees(component, node)
        if gray > 1:
            return True
    return False

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


class scenario:

    def __init__(self, surrounding, max_rearranged, last_rearranged, rearranged = [], ordered_partition = None, gray = None, proportion_chromo = None, deletions = None, duplications = None):
        self.surrounding = surrounding
        self.proportion_chromo = proportion_chromo
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


def subgraph_of_breakpoint_graph(nxcomponents, maximum, component_chr):

    graph=nx.MultiGraph()
    for i in range(len(nxcomponents)):
        if component_chr[i] & maximum:
            graph = nx.union(graph, copy.deepcopy(nxcomponents[i]))
    return graph

