#! /bin/python3

import itertools
from lib import breakpointgraph
from lib import jabba_preporcessing
from lib import auxilliary
from lib import samples



GRAPHS = 'jabba_data/graphs/'

#        _______________________________
#_______/       CONSTANTS               \_____________________________________________


LIST_OF_ANNOTATIONS = ['chromothripsis']

MIN_NUMBER_OF_CHROMOSOMES = 2

PROPORTION_OF_MISSANOTATION = 0.2

TEMPLATED_INSERTIONS = False

DUPLICATION_LENGTH = 1000

DELETION_LENGTH = 100000

NUMBER_OF_COMPONENTS = 2






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
    rearranged, genome_graph = breakpointgraph.run_a_partition(ordered_partition, nxcomponents, bg, TEMPLATED_INSERTIONS)
    affected = breakpointgraph.number_of_affected_chromosomes_by_nodes(genome_graph, annotated_nodes)
    
    
    rearranged.append(affected)

    if affected < number_of_affected_chromosomes:
        return rearranged, affected
    return False, False
        
#keep only circles
def filter_components(nxcomponents, component_chr, maximum):
    circles_indexes = [] 
    for i in range(len(nxcomponents)):
        number_of_gray = breakpointgraph.number_of_gray_edges(nxcomponents[i])
        if component_chr[i] & maximum and number_of_gray > 0:
            if breakpointgraph.circle(nxcomponents[i]) and number_of_gray>=2:
                circles_indexes.append(i)
    return circles_indexes

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


def stats_ordered_partition(nxcomponents, ordered_partition, iidTOlength):
    surrounding = []
    deletions = []
    gray = []
    duplications = []
    for indexes in ordered_partition:
        surrounding.append(breakpointgraph.surrounding_blocks_components(nxcomponents, indexes, iidTOlength))
        deletions.append(breakpointgraph.length_of_deletions(indexes, nxcomponents))
        duplications.append(breakpointgraph.length_of_duplications(indexes, nxcomponents))
        number_of_gray = 0
        for index in indexes:
            number_of_gray += breakpointgraph.number_of_gray_edges(nxcomponents[index])
        gray.append(number_of_gray) 
    return surrounding, deletions, duplications, gray


def check_if_scenario_more_optimal(nxcomponents, iidTOlength, affected, ordered_partition, rearranged, optimal_scenario):
    surrounding, deletions, duplications, numbers_of_gray_edges = stats_ordered_partition(nxcomponents, ordered_partition, iidTOlength)

    if (affected, max(rearranged)) < (optimal_scenario.affected, max(optimal_scenario.rearranged)):
        optimal_scenario = samples.scenario(rearranged, affected, ordered_partition, numbers_of_gray_edges, deletions, duplications, surrounding)
    return optimal_scenario

def check_if_enough_to_remove(affected, ordered_partition, nxcomponents, iidTOlength, optimal_multibreak, proportion_chromo):
    surrounding, deletions, duplications, numbers_of_gray_edges = stats_ordered_partition(nxcomponents, ordered_partition, iidTOlength)
    if len(affected) < len(optimal_multibreak.affected):
        return samples.multibreak(affected, ordered_partition, numbers_of_gray_edges, proportion_chromo, deletions, duplications, surrounding)
    return optimal_multibreak


#        _________________________________
#_______/       STUDY FILES               \_____________________________________________



def study_files_with_complex_rearrangements(files_of_interest):

    rearrange = []
    remove = []

    for file in files_of_interest:
        data = jabba_preporcessing.read_json(GRAPHS+file)
        sample = study_sample(data, file)
        
        for annotation in sample.annotations:
            if annotation.optimal_scenario != None:
                rearrange.append(annotation)

            if annotation.optimal_multibreak != None:
                remove.append(annotation)

    print("")
    print("There are {} complex rearrangements for which there exists a ISA multi-break that transforms the reference genome into a genome in which the complex rearrangemt affects less chromosomes".format(len(rearrange)))
  
    index = 0
    for annotation in rearrange:
        index+=1
        scenario = annotation.optimal_scenario
        print("")
        print("{}. filename: {}, type: {}, dataset: {}".format(index, annotation.filename, annotation.type, annotation.dataset))
        print("annotation: {}, affected chromosomes: {}, chromosome number: {}".format(annotation.title, annotation.chromosomes, scenario.affected))
        print("extremities of the identified ISA {}-break: {}".format(int(sum(scenario.numbers_of_gray_edges)), annotation.nodes_of_multibreak))
        print("number of the deleted nucleotides: {}, the number of duplicated nucleotides: {}".format(sum(scenario.deletions), sum(scenario.duplications)))
        print("chromosomes in a chromosome graph component: {}".format(annotation.chromosomes_in_a_chromosome_graph_component))

    print("")
    print("There are {} complex rearrangements for which there exists a ISA multi-break that introduces adjacencies potentially missannotated as having been introduced by the complex rearrangenent.".format(len(remove)))

    index = 0
    for annotation in remove:
        index += 1
        scenario = annotation.optimal_multibreak
        print("")
        print("{}. filename: {}, type: {}, dataset: {}".format(index, annotation.filename, annotation.type, annotation.dataset))
        print("annotation: {}, affected chromosomes: {}".format(annotation.title, annotation.chromosomes))
        print("proposed affected chromosomes: {}, chromosome number: {}".format(scenario.affected, len(scenario.affected)))
        print("extremities of the identified ISA {}-break: {}".format(int(sum(scenario.numbers_of_gray_edges)), annotation.nodes_of_multibreak))
        print("number of the deleted nucleotides: {}, the number of duplicated nucleotides: {}".format(sum(scenario.deletions), sum(scenario.duplications)))
        print("chromosomes in a chromosome graph component: {}".format(annotation.chromosomes_in_a_chromosome_graph_component))


  


#        _________________________________
#_______/       STUDY SAMPLE               \_____________________________________________



def study_sample(data, filename):

    #find the chromosomes affected by every complex rearrangement. 
    anotationsTOchr = jabba_preporcessing.anotations(data, LIST_OF_ANNOTATIONS)
    #construct breakpoint graph and the lenghts of the synteny blocks
    bg, iidTOlength = jabba_preporcessing.jabba_to_breakpoint_graph(data, DUPLICATION_LENGTH, filename, DELETION_LENGTH)

    sample = samples.sample(filename, bg.graph['type'],bg.graph['dataset'])

    #identify the connected components of the breakpoint graph
    nxcomponents = breakpointgraph.graphTOcomponents(bg)
    #find what chromosomes are affeced by the connected components and what nodes do those annotations affect
    component_chr, annotationTOnodes, annotationTOcomponents =  breakpointgraph.annotation_of_bg_components(nxcomponents)
    #partition the chromosomes into subsets {M_1,...,M_p} such that for every connected component in the breakpoint grahp
    #all of its vertices belong to the chromosomes in a single subset M_i. 
    maximums, annotations_of_maximums = breakpointgraph.subsets_of_related_chromosomes(anotationsTOchr,component_chr,MIN_NUMBER_OF_CHROMOSOMES )        
    
    for j in range(len(maximums)):
        #find indexes of the connected components of the breakpoint graph that are circles (called cycles in the paper) and 
        #whose vertices are among chromosomes in maximums[j]
        indexes = filter_components(nxcomponents, component_chr, maximums[j])
        for annotation in annotations_of_maximums[j]:
        
            
            complex_rearrangement = samples.annotation(annotation, filename, bg.graph['type'],bg.graph['dataset'],maximums[j])
            optimal_scenario = samples.scenario([float('inf')], len(anotationsTOchr[annotation]))
            optimal_multibreak = samples.multibreak(anotationsTOchr[annotation])
            add = False
                    
            #iterate over subsets of the connectected components of the breakpoint graph of size at most NUMBER_OF_COMPONENTS 
            for number_of_components in range(NUMBER_OF_COMPONENTS):
                for components_of_choice in itertools.combinations(indexes, number_of_components+1):

                    annotated_nodes = annotationTOnodes[annotation]
                    annotated_nodes_minus_multi_break_nodes, intersecting_nodes = annotated_nodes_minus_chosen_nodes(annotationTOnodes[annotation], components_of_choice, nxcomponents)                      
                    #proportion of gray edges in the graph that are annotated as introduced by the complex rearrangement and also included in the components_of_choice
                    proportion_chromo = len(intersecting_nodes)/len(annotationTOnodes[annotation])
                        
                    ordered_partition = [components_of_choice]
                    #check what happens if the components of choice are removed
                    affected = nodes_to_affected_chromosomes(annotated_nodes_minus_multi_break_nodes)
                    if len(affected) < len(anotationsTOchr[annotation]) and proportion_chromo<= PROPORTION_OF_MISSANOTATION:
                        optimal_multibreak = check_if_enough_to_remove(affected, ordered_partition, nxcomponents, iidTOlength, optimal_multibreak, proportion_chromo)
                        
                    if len(intersecting_nodes) == 0:
                        rearranged, affected = check_ordered_partition(ordered_partition, nxcomponents, bg, annotated_nodes, len(anotationsTOchr[annotation]))
                        if rearranged:

                            optimal_scenario = check_if_scenario_more_optimal(nxcomponents, iidTOlength, affected, ordered_partition, rearranged, optimal_scenario)
                
                            if number_of_components > 0:
                                #iterate through the rest of ordered partitions
                                for ordered_partition in auxilliary.order(list(components_of_choice)):
                                    rearranged, affected = check_ordered_partition(ordered_partition, nxcomponents, bg, annotated_nodes, len(anotationsTOchr[annotation]))
                                    if rearranged:
                                        optimal_scenario = check_if_scenario_more_optimal(nxcomponents, iidTOlength, affected, ordered_partition, rearranged, optimal_scenario) 
                            

            if optimal_scenario.affected < len(anotationsTOchr[annotation]):
                complex_rearrangement.nodes_of_multibreak = breakpointgraph.extremities_in_components(nxcomponents,ordered_partition_to_indexes(optimal_scenario.ordered_partition))
                complex_rearrangement.chromosomes = anotationsTOchr[annotation]
                complex_rearrangement.optimal_scenario = optimal_scenario
                add = True


            if len(optimal_multibreak.affected) < len(anotationsTOchr[annotation]):
                complex_rearrangement.nodes_of_multibreak = breakpointgraph.extremities_in_components(nxcomponents,ordered_partition_to_indexes(optimal_multibreak.ordered_partition))
                complex_rearrangement.chromosomes = anotationsTOchr[annotation]
                complex_rearrangement.optimal_multibreak = optimal_multibreak
                add = True

            if add:
                sample.annotations.append(complex_rearrangement)

    return sample
    
#select the JabBaA graphs from GRAPHS that contain a complex rearrangement from LIST_OF_ANNOTATIONS that affects more than MIN_NUMBER_OF_CHROMOSOMES
files_of_interest = jabba_preporcessing.read_annotations(GRAPHS, MIN_NUMBER_OF_CHROMOSOMES, LIST_OF_ANNOTATIONS)


scenarios = study_files_with_complex_rearrangements(files_of_interest)
