#! /bin/python3

#WE ASSUME THAT NEIGHBORING SEGMENTS in a JaBbaA graph HAVE iid's THAT DIFFER BY ONE!

import json
import networkx as nx
import os
from collections import defaultdict



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


def read_annotations(dir_name, min_number_of_chromosomes, list_of_annotations):

    files_of_interest = [] 
    total_files = 0
    total_annotations = 0

    for json_file in os.listdir(dir_name):
        
        total_files +=1 

        data = read_json(dir_name+json_file)

        anotationsTOchr = anotations(data, list_of_annotations)

        add = False
        for mutation in anotationsTOchr:
            if len(anotationsTOchr[mutation]) >= min_number_of_chromosomes:
                add = True
                total_annotations += 1

        if add:
            files_of_interest.append(json_file)

    print("number of JaBbA graphs in {}:  {}".format(dir_name,total_files))
    print("number of JaBbA graphs in {} annotated with a complex rearrangement from a list {} affecting at least {} chromsomes: {}".format(dir_name,list_of_annotations, min_number_of_chromosomes, len(files_of_interest)))
    print("total number of complex rearrangements from a list {} affecting at least {} chromsomes: {}".format(list_of_annotations, min_number_of_chromosomes, total_annotations))
   
    return files_of_interest

#Scans the input graph for the mutation annotations.
#For every a mutation, for example a "chromothripsis=1",
#it finds all the chromosomes including block extremities 
#incident to edges in this mutation.

def anotations(data, list_of_annotations):
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
                    if mutation.split("=")[0] in list_of_annotations:
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

            cn = interval['y']
            cn_left = interval_left['y']
            cn_right = interval_right['y']

            if right_gray + left_gray == 0 and tail_gray == 1 and head_gray == 1 and length_interval < duplication_length and cn_left < cn and cn_right < cn:
                bg.add_edge(left_head, right_tail, color ='brown', length = length_interval)

            elif right_gray == 1 and left_gray == 1 and tail_gray + head_gray == 0  and cn_left > cn and cn_right > cn:
                if length_interval < deletion_length:
                    if color_degree(bg, tail, 'brown') == 0 and color_degree(bg, head, 'brown') == 0:
                        bg.add_edge(tail, head, color ='red', length = length_interval)



def get_interval_length(interval):
    return interval['endPoint']-interval['startPoint'] + 1


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




