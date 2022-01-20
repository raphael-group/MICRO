#! /bin/python3

import itertools
from lib import breakpointgraph
from lib import jabba_preporcessing
from lib import auxilliary

#DUPLICATION_LENGTH = 1000, DELETION_LENGTH = 100000

LIST_OF_ANNOTATIONS = ['chromothripsis']
GRAPHS = 'jabba_data/graphs/'


TEMPLATED_INSERTIONS = False
DUPLICATION_LENGTH = 1000
DELETION_LENGTH = 100000
NUMBER_OF_COMPONENTS = 2


MIN_NUMBER_OF_CHROMOSOMES = 2



#        ______________________________________
#_______/       STUDY COMPONENTS               \_____________________________________________




def annotated_nodes_minus_chosen_nodes(annotated_nodes, indexes, nxcomponents):

    intersecting_nodes = set()
    pruned_annotated_nodes = set()
    circle_nodes = set()
    for index in indexes:
        circle_nodes = circle_nodes.union(nxcomponents[index].nodes())
                            
    for node in annotated_nodes:
        if node not in circle_nodes:
            pruned_annotated_nodes.add(node)
        else:
            intersecting_nodes.add(node)

    return pruned_annotated_nodes, intersecting_nodes

def check_ordered_partition(ordered_partition, nxcomponents, bg, annotated_nodes, number_of_affected_chromosomes):

    ordered_partition = list(ordered_partition)
    rearranged, genome_graph = breakpointgraph.run_a_partition_jabba(ordered_partition, nxcomponents, bg, TEMPLATED_INSERTIONS)
    last_affected = breakpointgraph.number_of_affected_chromosomes_by_nodes(genome_graph, annotated_nodes)
    
    rearranged.append(last_affected)

    if last_affected < number_of_affected_chromosomes:
        return rearranged
    return False
        
#keep only circles
def filter_components(nxcomponents, component_chr, maximum):
    circles_indexes = [] 
    for i in range(len(nxcomponents)):
        number_of_gray = jabba_preporcessing.number_of_gray_edges(nxcomponents[i])
        if component_chr[i] & maximum and number_of_gray > 0:
            if jabba_preporcessing.circle(nxcomponents[i]) and number_of_gray>=2:
                circles_indexes.append(i)
    return circles_indexes


def stats_ordered_partition(nxcomponents, ordered_partition, iidTOlength):
    surrounding = []
    deletions = []
    gray = []
    duplications = []
    for indexes in ordered_partition:
        surrounding.append(jabba_preporcessing.surrounding_blocks_components(nxcomponents, indexes, iidTOlength))
        deletions.append(breakpointgraph.length_of_deletions(indexes, nxcomponents))
        duplications.append(breakpointgraph.length_of_duplications(indexes, nxcomponents))
        number_of_gray = 0
        for index in indexes:
            number_of_gray += jabba_preporcessing.number_of_gray_edges(nxcomponents[index])
        gray.append(number_of_gray) 
    return surrounding, deletions, duplications, gray


def nodes_to_affected_chromosomes(nodes):
    chromosomes = set()
    for node in nodes:
        chromosomes.add(node[1])
    return chromosomes

def ordered_partition_to_indexes(ordered_partition):
    indexes = []
    for subset in ordered_partition:
        for element in subset:
            indexes.append(element)
    return(indexes)


#        ______________________________________
#_______/       STUDY SCNEARIOS               \_____________________________________________


def check_if_scenario_more_optimal(nxcomponents, iidTOlength, ordered_partition, numbers_of_rearranged_chromosomes, scenario_max_surr, proportion_chromo):
    surrounding, deletions, duplications, nodes = stats_ordered_partition(nxcomponents, ordered_partition, iidTOlength)


    last_rearranged = numbers_of_rearranged_chromosomes[-1]
    max_rearranged = max(numbers_of_rearranged_chromosomes)
    if (last_rearranged, max_rearranged) < (scenario_max_surr.last_rearranged, scenario_max_surr.max_rearranged):
        scenario_max_surr = jabba_preporcessing.scenario(surrounding, max_rearranged, last_rearranged, numbers_of_rearranged_chromosomes, ordered_partition, nodes, proportion_chromo, deletions, duplications)
    return scenario_max_surr

def check_if_enough_to_remove(affected, ordered_partition, nxcomponents, iidTOlength, scenario_no_intersection_max_surr, proportion_chromo):
    surrounding, deletions, duplications, gray = stats_ordered_partition(nxcomponents, ordered_partition, iidTOlength)
    if len(affected) < scenario_no_intersection_max_surr.max_rearranged:
        return jabba_preporcessing.scenario(surrounding, len(affected), affected, [len(affected)], ordered_partition, gray, proportion_chromo, deletions, duplications)
    return scenario_no_intersection_max_surr

#        _________________________________
#_______/       STUDY FILES               \_____________________________________________



def study_files_with_chromothripsis(files_of_interest):

    rearrange = []
    remove = []

    for file in files_of_interest:
        data = jabba_preporcessing.read_json(GRAPHS+file)
        sample = study_sample(data, file)
        
        for annotation in sample.annotations:
            if annotation.scenario_max_sur != None:
                rearrange.append(annotation)

            if annotation.scenario_removal_max_sur != None:
                remove.append(annotation)

    print("")
    print("There are {} complex rearrangements for which there exists an ISA multi-break transforming the reference genome into a genome in which that complex rearrangemt affects less chromosomes".format(len(rearrange)))
  
    index = 0
    for annotation in rearrange:
        index+=1
        scenario = annotation.scenario_max_sur
        print("")
        print("{}. filename: {}, type: {}, dataset: {}".format(index, annotation.filename, annotation.type, annotation.dataset))
        print("annotation: {}, affected chromosomes: {}, chromosome number: {}".format(annotation.title, annotation.chromosomes, scenario.rearranged[-1]))
        print("extremities of the ISA multi-break:", annotation.nodes_of_multibreak)
        print("number of the deleted nucleotides: {}, the number of duplicated nucleotides: {}".format(sum(scenario.deletions), sum(scenario.duplications)))

    print("")
    print("There are {} complex rearrangements for which some of the novel adjacencies are missannotated as having been introduced by that complex rearrangenent.".format(len(remove)))

    index = 0
    for annotation in remove:
        index += 1
        scenario = annotation.scenario_removal_max_sur
        print("")
        print("{}. filename: {}, type: {}, dataset: {}".format(index, annotation.filename, annotation.type, annotation.dataset))
        print("annotation: {}, affected chromosomes: {}".format(annotation.title, annotation.chromosomes, scenario.max_rearranged))
        print("proposed affected chromosomes: {}, chromosome number: {}".format(scenario.last_rearranged, scenario.max_rearranged))
        print("extremities of the ISA multi-break:", annotation.nodes_of_multibreak)
        print("number of the deleted nucleotides: {}, the number of duplicated nucleotides: {}".format(sum(scenario.deletions), sum(scenario.duplications)))

  


#        _________________________________
#_______/       STUDY SAMPLE               \_____________________________________________



def study_sample(data, filename):

    #construct breakpoint graph
    bg, iidTOlength = jabba_preporcessing.jabba_to_breakpoint_graph(data, DUPLICATION_LENGTH, filename, DELETION_LENGTH)
    #find connected components of the breakpoint graph
    nxcomponents = breakpointgraph.breakpointgraphTOcomponents(bg)
    #find what chromosomes are affeced by the connected components and what nodes do those annotations affect
    component_chr, annotationTOnodes, annotationTOcomponents =  jabba_preporcessing.annotation_of_bg_components(nxcomponents)
    #for every annotation of an event, find what chromosomes does it affect
    anotationsTOchr = jabba_preporcessing.anotations(data, LIST_OF_ANNOTATIONS)
    #find the connected components of the multi-genome graph that include annotated events
    #maximums contains sets of chromosomes that are affected by these connected components
    #while annotation_size contains the maximum number of chromosomes in these components affected by some annotated event 
    maximums, annotations_of_maximums = jabba_preporcessing.components_with_annotations(anotationsTOchr,component_chr,MIN_NUMBER_OF_CHROMOSOMES )        
    #iterate through the connected components of the multi-genome graph
    sample = jabba_preporcessing.sample(filename, bg.graph['type'],bg.graph['dataset'])
    
    for j in range(len(maximums)):
        indexes = filter_components(nxcomponents, component_chr, maximums[j])
            
        for annotation in annotations_of_maximums[j]:

            chromothripsis = jabba_preporcessing.annotation(annotation, filename, bg.graph['type'],bg.graph['dataset'])

            interesting_annotation = False
            removed_interesting = False
            scenario_max_surr = jabba_preporcessing.scenario([0], float('inf'), float('inf'))
            scenario_no_intersection_max_surr = jabba_preporcessing.scenario([0], float('inf'), float('inf'), 0,0,0,0,0,0)
                
            for number_of_components in range(NUMBER_OF_COMPONENTS):
                for components_of_choice in itertools.combinations(indexes, number_of_components+1):
                    annotated_nodes, intersecting_nodes = annotated_nodes_minus_chosen_nodes(annotationTOnodes[annotation], components_of_choice, nxcomponents)                      
                    #proportion of gray edges in the graph that are annotated as chromothripsis and also included in the components_of_choice
                    proportion_chromo = len(intersecting_nodes)/len(annotationTOnodes[annotation])
                    #proportion of the gray edges in components_of_choice that are annotated as chromothripsis
                        
                    ordered_partition = [components_of_choice]

                    #check what happens if components of choice are removed
                    affected = nodes_to_affected_chromosomes(annotated_nodes)
                    if len(affected) < len(anotationsTOchr[annotation]) and proportion_chromo<=0.2:
                        removed_interesting = True
                        scenario_no_intersection_max_surr = check_if_enough_to_remove(affected, ordered_partition, nxcomponents, iidTOlength, scenario_no_intersection_max_surr, proportion_chromo)
                        
                    if len(intersecting_nodes) == 0:
                        rearranged = check_ordered_partition(ordered_partition, nxcomponents, bg, annotated_nodes, len(anotationsTOchr[annotation]))
                        if rearranged:
                            interesting_annotation = True
                            scenario_max_surr = check_if_scenario_more_optimal(nxcomponents, iidTOlength,ordered_partition, rearranged, scenario_max_surr, proportion_chromo)
                                
                            
                            if number_of_components > 0:
                                #iterate through the rest of ordered partitions
                                for ordered_partition in auxilliary.order(list(components_of_choice)):
                                    rearranged = check_ordered_partition(ordered_partition, nxcomponents, bg, annotated_nodes, len(anotationsTOchr[annotation]))
                                    if rearranged:
                                        scenario_max_surr = check_if_scenario_more_optimal(nxcomponents, iidTOlength,ordered_partition, rearranged, scenario_max_surr, proportion_chromo) 
                            

            if interesting_annotation: 
                chromothripsis.nodes_of_multibreak = jabba_preporcessing.extremities_in_components(nxcomponents,ordered_partition_to_indexes(scenario_max_surr.ordered_partition))
                chromothripsis.chromosomes = anotationsTOchr[annotation]
                chromothripsis.scenario_max_sur = scenario_max_surr


            if removed_interesting:
                chromothripsis.nodes_of_multibreak = jabba_preporcessing.extremities_in_components(nxcomponents,ordered_partition_to_indexes(scenario_no_intersection_max_surr.ordered_partition))
                chromothripsis.chromosomes = anotationsTOchr[annotation]
                chromothripsis.scenario_removal_max_sur = scenario_no_intersection_max_surr



            sample.annotations.append(chromothripsis)

    #breakpointgraph.plotgraphs([nxcomponents[62],nxcomponents[265],nxcomponents[331],nxcomponents[394],nxcomponents[1720]], filename)
    #print("1",jabba_preporcessing.extremities_in_components(nxcomponents, [270]))

    return sample
    
files_of_interest, typeTOtotalnumber, typeTOnumberwithchromothripsis = jabba_preporcessing.read_annotations(GRAPHS, MIN_NUMBER_OF_CHROMOSOMES, LIST_OF_ANNOTATIONS)

scenarios = study_files_with_chromothripsis(files_of_interest)