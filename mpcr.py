#! /bin/python3

import itertools
from collections import defaultdict

from lib import breakpointgraph
from lib import jabba_preporcessing
from lib import auxilliary
import matplotlib.pyplot as plt

import statistics
import os

import numpy as np


#Rearrange = True, Circles = True, DUPLICATION_LENGTH = 1000, DELETION_LENGTH = 100000


REARRANGE = True
CIRCLES = True

TEMPLATED_INSERTIONS = False

DUPLICATION_LENGTH = 1000
DELETION_LENGTH = 100000

NUMBER_OF_COMPONENTS = 1
SURROUNDING_TRESHOLD = 100000
BOUND_ON_TOTAL_NUMBER_OF_COMPONENTS = 10000

def study_edges(files_of_interest):
    total_ratio = []
    total_circles = []
    total_max_break = [] 
    total_sum = []
    empty = 0
    for file in files_of_interest:
        data = jabba_preporcessing.read_json('jabba_data/graphs/'+file)
        bg, iidTOlength, black, gray, added_edges = jabba_preporcessing.jabba_to_breakpoint_graph(data, file)
        if black >0:
            total_ratio.append((gray+added_edges)/black)
        else:
            empty+=1
        nxcomponents = breakpointgraph.breakpointgraphTOcomponents(bg)
        circles = 0
        big = 0

        sum = 0
        for component in nxcomponents:
            if len(component.nodes()) > 4:
                if jabba_preporcessing.circle(component):
                    circles+=1
                big = max(big, jabba_preporcessing.number_of_gray_edges(component))
                if big>2:
                    sum+=big

        total_circles.append(circles)
        total_max_break.append(big)
        total_sum.append(sum)

       
    total_circles.sort()
    total_ratio.sort()
    total_max_break.sort()
    total_sum.sort()

    print(total_sum)
    print(statistics.median(total_sum))

    #print('median black/gray ratio, number of circles, max k-break:', statistics.median(total_ratio), statistics.median(total_circles), statistics.median(total_max_break))
    #print('min black/gray ratio, number of circles, max k-break:', min(total_ratio), min(total_circles), min(total_max_break))
    #print('max black/gray ratio, number of circles, max k-break:', max(total_ratio), max(total_circles), max(total_max_break))

    no_3_break = 0
    for number in total_max_break:
        if number<=2:
            no_3_break+=1
    print("Number of files without breaks", no_3_break)
    print("number of empty files", empty)


#        ______________________________________
#_______/       STUDY COMPONENTS               \_____________________________________________




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

def annotated_nodes_minus_chosen_nodes(annotated_nodes, indexes, nxcomponents, bg):

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

    #find the nodes of the components that have gray degree > 0
    gray_circle_nodes= []
    for node in circle_nodes:
        gray, black = jabba_preporcessing.vertex_degrees(bg, node)
        if gray > 0:
            gray_circle_nodes.append(node)

    return pruned_annotated_nodes, intersecting_nodes, gray_circle_nodes

def check_ordered_partition(ordered_partition, nxcomponents, bg, annotated_nodes, number_of_affected_chromosomes):

    ordered_partition = list(ordered_partition)
    rearranged, genome_graph = breakpointgraph.run_a_partition_jabba(ordered_partition, nxcomponents, bg, TEMPLATED_INSERTIONS)
    last_affected = breakpointgraph.number_of_affected_chromosomes_by_nodes(genome_graph, annotated_nodes)
    
    rearranged.append(last_affected)


    if last_affected < number_of_affected_chromosomes:# and max(rearranged)<= number_of_affected_chromosomes:
        return rearranged
    else: 
        return False
        

def filter_components(nxcomponents, component_chr, maximum):
    circles_indexes = [] 
    multi_break_indexes = []
    component_indexes = []

    total_gray = 0
    for i in range(len(nxcomponents)):
        number_of_gray = jabba_preporcessing.number_of_gray_edges(nxcomponents[i])
        total_gray += number_of_gray
        if component_chr[i] & maximum and number_of_gray > 0:
            component_indexes.append(i)
            if jabba_preporcessing.circle(nxcomponents[i]) and number_of_gray>=2:
                circles_indexes.append(i)
            if number_of_gray == 1: #multi_break_component(nxcomponents[i]) 
                multi_break_indexes.append(i)
    return circles_indexes, component_indexes, multi_break_indexes, total_gray


def multi_break_component(component):
    for edge in component.edges.data():
        if edge[2]['color'] == 'yellow':
            return False
    for node in component.nodes():
        if component.degree(node) > 2:
            return False
    return True


def stats_ordered_partition(nxcomponents, ordered_partition, iidTOlength):
    surrounding = []
    deletions = []
    gray = []
    duplications = []
    for indexes in ordered_partition:
        surrounding.append(jabba_preporcessing.surrounding_blocks_components(nxcomponents, indexes, iidTOlength))
        deletions.append(length_of_deletions(indexes, nxcomponents))
        duplications.append(length_of_duplications(indexes, nxcomponents))
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


#        ______________________________________
#_______/       STUDY SCNEARIOS               \_____________________________________________

def max_min_surrounding_min_last_rearranged_min_max_rearrnanged(min_surrounding, last_rearranged, max_rearranged, scenario):
    if min_surrounding > min(scenario.surrounding):
        return True

    if min_surrounding == min(scenario.surrounding):
        if (last_rearranged, max_rearranged) < (scenario.last_rearranged, scenario.max_rearranged):
            return True
            
    return False


def min_last_rearranged_min_max_rearrnanged_max_min_surrounding(min_surrounding, last_rearranged, max_rearranged, scenario):

    if (last_rearranged, max_rearranged) < (scenario.last_rearranged, scenario.max_rearranged):
        return True
        
    if (last_rearranged, max_rearranged) == (scenario.last_rearranged, scenario.max_rearranged):
        if min_surrounding > min(scenario.surrounding):
            return True
            
    return False


def min_last_rearranged_max_min_surrounding(min_surrounding, last_rearranged, scenario):

    if last_rearranged < scenario.last_rearranged:
        return True
        
    if last_rearranged == scenario.last_rearranged:
        if min_surrounding > min(scenario.surrounding):
            return True
            
    return False


def check_if_scenario_more_optimal(ordered_partition, numbers_of_rearranged_chromosomes, scenario_max_surr, scenario_min_rear, scenario_min_last_rear, surrounding, deletions, duplications, nodes, proportion_chromo, proportion_components):
    last_rearranged = numbers_of_rearranged_chromosomes[-1]
    max_rearranged = max(numbers_of_rearranged_chromosomes)

    if max_min_surrounding_min_last_rearranged_min_max_rearrnanged(min(surrounding), last_rearranged, max_rearranged, scenario_max_surr):
        scenario_max_surr = jabba_preporcessing.scenario(surrounding, max_rearranged, last_rearranged, numbers_of_rearranged_chromosomes, ordered_partition, nodes, proportion_chromo, proportion_components, deletions, duplications)

    if min_last_rearranged_min_max_rearrnanged_max_min_surrounding(min(surrounding), last_rearranged, max_rearranged, scenario_min_rear):
        scenario_min_rear = jabba_preporcessing.scenario(surrounding, max_rearranged, last_rearranged, numbers_of_rearranged_chromosomes, ordered_partition, nodes, proportion_chromo, proportion_components, deletions, duplications)

    if min_last_rearranged_max_min_surrounding(min(surrounding), last_rearranged, scenario_min_last_rear):
        scenario_min_last_rear = jabba_preporcessing.scenario(surrounding, max_rearranged, last_rearranged, numbers_of_rearranged_chromosomes, ordered_partition, nodes, proportion_chromo, proportion_components, deletions, duplications)


    return scenario_max_surr, scenario_min_rear, scenario_min_last_rear


def check_if_enough_to_remove(affected, ordered_partition, nxcomponents, iidTOlength, scenario_no_intersection_max_surr, proportion_chromo, proportion_components):
    surrounding, deletions, duplications, gray = stats_ordered_partition(nxcomponents, ordered_partition, iidTOlength)
    if min(surrounding) > min(scenario_no_intersection_max_surr.surrounding):
        return jabba_preporcessing.scenario(surrounding, len(affected), len(affected), [len(affected)], ordered_partition, gray, proportion_chromo, proportion_components, deletions, duplications)
    return scenario_no_intersection_max_surr

#        ___________________________
#_______/       PLOTS               \_____________________________________________


def plot_chromothripsis(scenarios):
    gray = []
    surrounding = []
    proportion_chromo = []
    proportion_components = []
    large_surrounding = 0
    for scenario in scenarios:
        gray.append(sum(scenario.gray))
        surrounding.append(min(scenario.surrounding))
        if min(scenario.surrounding)>SURROUNDING_TRESHOLD:
            large_surrounding+=1
        proportion_chromo.append(scenario.proportion_chromo)
        proportion_components.append(scenario.proportion_components)

    
    #plt.figure()
    #plt.scatter(surrounding, gray, c=proportion_chromo, alpha=0.3, cmap='viridis')
    #cbar = plt.colorbar();  # show color scale
    #plt.xlabel("separation of the selected componentes")
    #plt.xscale('log')
    #plt.ylabel("# of nodes in the selected components")
    #plt.title('of chromothripsis of interest = {}, # of chromotripsis with scenario surrounding > {} = {}, \n max # of components selected = {}, max # number of components studied = {}'.format(len(scenarios),SURROUNDING_TRESHOLD,large_surrounding, NUMBER_OF_COMPONENTS, BOUND_ON_TOTAL_NUMBER_OF_COMPONENTS))
    #plt.savefig("chromo-proportion-{}-{}-{}.pdf".format(NUMBER_OF_COMPONENTS, BOUND_ON_TOTAL_NUMBER_OF_COMPONENTS, DELETION_LENGTH), bbox_inches='tight')

    plt.figure()
    plt.scatter(surrounding, gray, c=proportion_components, alpha=0.3, cmap='viridis')
    cbar = plt.colorbar();  # show color scale
    cbar.set_label('% of gray edges annotated as chromothripsis', rotation=270, labelpad=20)
    plt.xlabel("separation")
    plt.xscale('log')
    plt.ylabel("# of gray edges")
    plt.title('# of chromothripsis of interest = {}, # of chromotripsis with scenario surrounding > {} = {}, \n # of components allowed = {}, # number of components studied = {}'.format(len(scenarios),SURROUNDING_TRESHOLD,large_surrounding, NUMBER_OF_COMPONENTS, BOUND_ON_TOTAL_NUMBER_OF_COMPONENTS))
    plt.savefig("annotations-{}-{}-{}.pdf".format(NUMBER_OF_COMPONENTS, BOUND_ON_TOTAL_NUMBER_OF_COMPONENTS, DELETION_LENGTH), bbox_inches='tight')
    

    #plt.figure()
    #plt.scatter(surrounding, gray, c=proportion_components, alpha=0.3, cmap='viridis')
    #plt.colorbar();  # show color scale
    #plt.xlabel("separation")
    #plt.xscale('log')
    #plt.xlim(0, 500000)
    #plt.ylabel("# of gray edges")
    #plt.title('# of chromothripsis of interest = {}, # of chromotripsis with scenario surrounding > {} = {}, \n # of components allowed = {}, # number of components studied = {}'.format(len(scenarios),SURROUNDING_TRESHOLD,large_surrounding, NUMBER_OF_COMPONENTS, BOUND_ON_TOTAL_NUMBER_OF_COMPONENTS))
    #plt.savefig("linear_scale-{}-{}.pdf".format(NUMBER_OF_COMPONENTS, BOUND_ON_TOTAL_NUMBER_OF_COMPONENTS), bbox_inches='tight')


def plot_samples(typeTOtotalnumber,typeTOnumberwithchromothripsis,typeTOnumberwithinterestingchromothripsis,typeTOnumberwithinterestingseparatedchromothripsis): 
    chromothripsis = [sum(typeTOnumberwithchromothripsis.values())/sum(typeTOtotalnumber.values())]
    interestingchromothripsis = [sum(typeTOnumberwithinterestingchromothripsis.values())/sum(typeTOtotalnumber.values())]
    interestingseparatedchromothripsis = [sum(typeTOnumberwithinterestingseparatedchromothripsis.values())/sum(typeTOtotalnumber.values())]
    name = ["ALL_SAMPLES ("+str(sum(typeTOtotalnumber.values()))+','+str(sum(typeTOnumberwithchromothripsis.values()))+','+str(sum(typeTOnumberwithinterestingchromothripsis.values()))+','+str(sum(typeTOnumberwithinterestingseparatedchromothripsis.values()))+')']

    for key in typeTOtotalnumber:
        #meso has large proportion of chromothripsis but only 4 samples or so, messes up with the y axis
        if key != "MESO":
            name.append(key+' ('+str(typeTOtotalnumber[key])+','+str(typeTOnumberwithchromothripsis[key])+','+str(typeTOnumberwithinterestingchromothripsis[key])+','+str(typeTOnumberwithinterestingseparatedchromothripsis[key])+')')
            chromothripsis.append(typeTOnumberwithchromothripsis[key]/typeTOtotalnumber[key])
            interestingchromothripsis.append(typeTOnumberwithinterestingchromothripsis[key]/typeTOtotalnumber[key])
            interestingseparatedchromothripsis.append(typeTOnumberwithinterestingseparatedchromothripsis[key]/typeTOtotalnumber[key])

    list1, list2, list3, list4 = zip(*reversed(sorted(zip(interestingseparatedchromothripsis, interestingchromothripsis, chromothripsis, name))))


    plt.figure()
    width = 0.35 
    fig = plt.subplots(figsize =(10, 7))
    p1 = plt.bar(list4, list1, width)
    p2 = plt.bar(list4, list2, width, bottom = list1)
    p3 = plt.bar(list4, list3, width, bottom = list2)

    plt.ylabel('Proportion of the samples')
    plt.xticks(list4, rotation='vertical')
    plt.legend((p3[0], p2[0], p1[0]), ('samples with chromothripsis affecting > 1 chromosome', 'samples of interest', 'separated samples of interest'))
    plt.title('max # of components selected = {}, max # number of components studied = {}'.format(NUMBER_OF_COMPONENTS, BOUND_ON_TOTAL_NUMBER_OF_COMPONENTS))
    plt.savefig("proportions-{}-{}.pdf".format(NUMBER_OF_COMPONENTS, BOUND_ON_TOTAL_NUMBER_OF_COMPONENTS), bbox_inches='tight')
    plt.close()

#        _________________________________
#_______/       STUDY FILES               \_____________________________________________



def test_deletions(files_of_interest):
    
    multi_break = []
    larger_than_2 = 0
    for filename in files_of_interest:
        data = jabba_preporcessing.read_json('jabba_data/graphs/'+filename)
        bg, iidTOlength, black, gray, added_edges, lengths, lengths_strange = jabba_preporcessing.jabba_to_breakpoint_graph(data, DUPLICATION_LENGTH, filename, REARRANGE, DELETION_LENGTH)
        nxcomponents = breakpointgraph.breakpointgraphTOcomponents(bg)

        max_component = 0
        for component in nxcomponents:
            max_component = max(max_component, jabba_preporcessing.number_of_gray_edges(component))
        
        if max_component > 2:
            larger_than_2 += 1

        multi_break.append(int(max_component))
    


    plt.bar(*np.unique(multi_break, return_counts=True))
    plt.xlabel("Maximum number of SVs in a subset")
    plt.ylabel("Number of samples")
    plt.savefig("deletion_max_number_of_svs-{}.pdf".format(DUPLICATION_LENGTH), bbox_inches='tight')
    plt.close()

    print(multi_break)
    print(len(files_of_interest), larger_than_2)




def study_files_with_chromothripsis(files_of_interest):
    typeTOnumberwithinterestingchromothripsis = defaultdict(int)
    typeTOnumberwithinterestingseparatedchromothripsis = defaultdict(int)

    multi_break, separated_multi_break = 0, 0 
    removal, separated_removal = 0, 0
    total_annotations = 0
    total_files = len(files_of_interest)
    scenarios = []

    interesting_files = 0
    separated_files = 0 

    total_components = 0
    total_multi = 0

    circles = []
    components_in_graph = []
    single_grays = []


    for file in files_of_interest:
        data = jabba_preporcessing.read_json('jabba_data/graphs/'+file)
        sample, components, multi, number_of_circles, number_of_components, single_gray = study_sample(data, file)
        
        circles.append(number_of_circles)
        components_in_graph.append(number_of_components)
        single_grays.append(single_gray)
        
        type = sample.type

        total_components+=components
        total_multi+=multi

        interesting = False
        interesting_separated = False

        for annotation in sample.annotations:
            total_annotations +=1
            if annotation.scenario_max_sur != None:
                scenarios.append(annotation.scenario_max_sur)
                multi_break +=1
                if interesting == False: 
                    typeTOnumberwithinterestingchromothripsis[type]+=1
                    interesting = True
                    interesting_files +=1

                if min(annotation.scenario_max_sur.surrounding) > SURROUNDING_TRESHOLD:
                    separated_multi_break +=1
                    if interesting_separated == False:
                        typeTOnumberwithinterestingseparatedchromothripsis[type]+=1
                        interesting_separated = True
                        separated_files +=1 

                elif annotation.scenario_removal_max_sur != None:
                    if min(annotation.scenario_removal_max_sur.surrounding) > SURROUNDING_TRESHOLD:# and annotation.scenario_removal_max_sur.proportion_chromo < 0.1:
                        separated_removal +=1
            elif annotation.scenario_removal_max_sur != None:
                #if not REARRANGE:
                scenarios.append(annotation.scenario_removal_max_sur)
                removal+=1
                if min(annotation.scenario_removal_max_sur.surrounding) > SURROUNDING_TRESHOLD:# and annotation.scenario_removal_max_sur.proportion_chromo < 0.1:
                    separated_removal +=1


    print("total annotations", total_annotations)    
    print(multi_break, separated_multi_break, removal, separated_removal)
    print("total files", total_files)
    print(interesting_files, separated_files)
    print("components, multi", total_components, total_multi)
    print(len(scenarios))


    return scenarios, typeTOnumberwithinterestingchromothripsis, typeTOnumberwithinterestingseparatedchromothripsis


#        _________________________________
#_______/       STUDY SAMPLE               \_____________________________________________



def check_the_number_of_nodes(indexes, nxcomponents, component_chr):
    for index in indexes:
        if len(nxcomponents[index].nodes())<=4 and len(component_chr[index]) == 1:         
            return False
    return True

def investigate_components(nxcomponents):
    circles = 0
    components = 0
    single_gray = 0
    for component in nxcomponents:
        number_of_gray = jabba_preporcessing.number_of_gray_edges(component)
        if jabba_preporcessing.circle(component) and number_of_gray>=2:
            circles += 1
        if number_of_gray >= 2 and not jabba_preporcessing.component_with_high_gray_degree(component):
            components+=1
        if number_of_gray == 1 and not jabba_preporcessing.component_with_high_gray_degree(component):
            single_gray += 1
    return circles, components, single_gray
        





def study_sample(data, filename):

    #construct breakpoint graph
    bg, iidTOlength, black, gray, added_edges, lengths, lengths_strange = jabba_preporcessing.jabba_to_breakpoint_graph(data, DUPLICATION_LENGTH, filename, REARRANGE, DELETION_LENGTH)
    #find connected components of the breakpoint graph
    nxcomponents = breakpointgraph.breakpointgraphTOcomponents(bg)
    #find what chromosomes are affeced by the connected components and what nodes do those annotations affect

    total_circles, total_components, single_gray = investigate_components(nxcomponents)
    component_chr, annotationTOnodes, annotationTOcomponents =  jabba_preporcessing.annotation_of_bg_components(nxcomponents)
    #for every annotation of an event, find what chromosomes does it affect
    anotationsTOchr = jabba_preporcessing.anotations(data)
    #find the connected components of the multi-genome graph that include annotated events
    #maximums contains sets of chromosomes that are affected by these connected components
    #while annotation_size contains the maximum number of chromosomes in these components affected by some annotated event 
    maximums, annotations_of_maximums = jabba_preporcessing.components_with_annotations(anotationsTOchr,component_chr)        
    #iterate through the connected components of the multi-genome graph
    sample = jabba_preporcessing.sample(filename, bg.graph['type'],bg.graph['dataset'])
    
    components = 0
    multi = 0
    
    for j in range(len(maximums)):
        circles_indexes, component_indexes, multi_break_indexes, total_gray = filter_components(nxcomponents, component_chr, maximums[j])

        components += len(component_indexes)
        multi += len(multi_break_indexes)

        if CIRCLES:
            indexes = circles_indexes
        elif REARRANGE:
            indexes = component_indexes
        else:
            indexes = multi_break_indexes


        if len(indexes) < BOUND_ON_TOTAL_NUMBER_OF_COMPONENTS or not REARRANGE:
            for annotation in annotations_of_maximums[j]:

                chromothripsis = jabba_preporcessing.annotation(annotation)

                interesting_annotation = False
                removed_interesting = False

                scenario_max_surr = jabba_preporcessing.scenario([0], float('inf'), float('inf'))
                scenario_min_rear = jabba_preporcessing.scenario([0], float('inf'), float('inf'))
                scenario_min_last_rear = jabba_preporcessing.scenario([0], float('inf'), float('inf'))
                scenario_no_intersection_max_surr = jabba_preporcessing.scenario([0], float('inf'), float('inf'), 0,0,0,0,0,0,0)
                
                for number_of_components in range(NUMBER_OF_COMPONENTS):
                    for components_of_choice in itertools.combinations(indexes, number_of_components+1):
                        annotated_nodes, intersecting_nodes, gray_circle_nodes = annotated_nodes_minus_chosen_nodes(annotationTOnodes[annotation], components_of_choice, nxcomponents, bg)                      
                        #proportion of gray edges in the graph that are annotated as chromothripsis and also included in the components_of_choice
                        proportion_chromo = len(intersecting_nodes)/len(annotationTOnodes[annotation])
                        #proportion of the gray edges in components_of_choice that are annotated as chromothripsis
                        proportion_components = len(intersecting_nodes)/len(gray_circle_nodes)
                        
                        ordered_partition = [components_of_choice]

                        
                        #check what happens if components of choice are removed
                        affected = nodes_to_affected_chromosomes(annotated_nodes)
                        if len(affected) < len(anotationsTOchr[annotation]) and proportion_chromo<=0.2:
                            removed_interesting = True
                            scenario_no_intersection_max_surr = check_if_enough_to_remove(affected, ordered_partition, nxcomponents, iidTOlength, scenario_no_intersection_max_surr, proportion_chromo, proportion_components)
                        
                        
                        if REARRANGE and len(intersecting_nodes) == 0:
                            if check_the_number_of_nodes(components_of_choice, nxcomponents, component_chr): #and jabba_preporcessing.circle(nxcomponents[components_of_choice[0]]):
                                rearranged = check_ordered_partition(ordered_partition, nxcomponents, bg, annotated_nodes, len(anotationsTOchr[annotation]))
                                if rearranged:
                                    interesting_annotation = True
                                    surrounding, deletions, duplications, gray = stats_ordered_partition(nxcomponents, ordered_partition, iidTOlength)
                                    scenario_max_surr, scenario_min_rear, scenario_min_last_rear = check_if_scenario_more_optimal(ordered_partition, rearranged, scenario_max_surr, scenario_min_rear, scenario_min_last_rear, surrounding, deletions, duplications, gray, proportion_chromo, proportion_components)

                                                                      
                                    #if number_of_components > 0:
                                        #iterate through the rest of ordered partitions
                                    #    for ordered_partition in auxilliary.order(list(components_of_choice)):
                                    #        rearranged = check_ordered_partition(ordered_partition, nxcomponents, bg, annotated_nodes, len(anotationsTOchr[annotation]))
                                    #        if rearranged:
                                    #            surrounding, deletions, duplications, gray = stats_ordered_partition(nxcomponents, ordered_partition, iidTOlength)
                                    #            scenario_max_surr, scenario_min_rear, scenario_min_last_rear  = check_if_scenario_more_optimal(ordered_partition, rearranged, scenario_max_surr, scenario_min_rear, scenario_min_last_rear, surrounding, deletions, duplications, gray, proportion_chromo, proportion_components) 

                        

                if (interesting_annotation or removed_interesting):#and scenario_no_intersection_max_surr.proportion_chromo <0.1: 
                   print('\n', NUMBER_OF_COMPONENTS,filename, bg.graph['type'], bg.graph['dataset'], annotation, anotationsTOchr[annotation])
            
                if interesting_annotation: 
                    print(scenario_max_surr.rearranged, scenario_max_surr.ordered_partition, scenario_max_surr.deletions,  scenario_max_surr.duplications)
                    #print("min_rearranged:", scenario_min_rear.surrounding, scenario_min_rear.rearranged, scenario_min_rear.ordered_partition, round(scenario_min_rear.proportion_chromo,2), round(scenario_min_rear.proportion_components,2), scenario_min_rear.gray,  scenario_min_rear.deletions, scenario_min_rear.duplications)
                    #print("min_last_rearranged:", scenario_min_last_rear.surrounding, scenario_min_last_rear.rearranged, scenario_min_last_rear.ordered_partition, round(scenario_min_last_rear.proportion_chromo,2), round(scenario_min_last_rear.proportion_components,2),scenario_min_last_rear.gray,  scenario_min_last_rear.deletions, scenario_min_last_rear.duplications)
                    chromothripsis.scenario_max_sur = scenario_max_surr
                    #print(total_circles, total_components)

                if removed_interesting:
                    print("scenario_no_intersection_max_surr:", scenario_no_intersection_max_surr.surrounding, scenario_no_intersection_max_surr.rearranged, scenario_no_intersection_max_surr.ordered_partition, round(scenario_no_intersection_max_surr.proportion_chromo,2), round(scenario_no_intersection_max_surr.proportion_components,2), scenario_no_intersection_max_surr.gray, scenario_no_intersection_max_surr.deletions, scenario_no_intersection_max_surr.duplications)
                    chromothripsis.scenario_removal_max_sur = scenario_no_intersection_max_surr
                    #print(total_circles, total_components)

                sample.annotations.append(chromothripsis)

    #breakpointgraph.plotgraphs([nxcomponents[62],nxcomponents[265],nxcomponents[331],nxcomponents[394],nxcomponents[1720]], filename)
    #print("1",jabba_preporcessing.extremities_in_components(nxcomponents, [270]))

    return sample, components, multi, total_circles, total_components, single_gray
    
files_of_interest, typeTOtotalnumber, typeTOnumberwithchromothripsis = jabba_preporcessing.read_annotations('jabba_data/graphs/')

scenarios, typeTOnumberwithinterestingchromothripsis, typeTOnumberwithinterestingseparatedchromothripsis = study_files_with_chromothripsis(files_of_interest)

