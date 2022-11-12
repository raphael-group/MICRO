#! /bin/python3

class sample:
    def __init__(self, filename, type, dataset):
        self.filename = filename
        self.type = type
        self.dataset = dataset
        self.annotations = []

class annotation:
    def __init__(self, title, filename, type, dataset, chromosomes_in_a_chromosome_graph_component):
        self.filename = filename
        self.type = type
        self.dataset = dataset
        self.title = title
        self.optimal_scenario = None
        self.optimal_multibreak = None
        self.nodes_of_multibreak = []
        self.chromosomes = []
        self.chromosomes_in_a_chromosome_graph_component = chromosomes_in_a_chromosome_graph_component



class scenario:

    def __init__(self, rearranged = [], affected = None, ordered_partition = None, numbers_of_gray_edges = None, deletions = None, duplications = None, surrounding = []):
        self.affected = affected
        self.rearranged = rearranged
        self.ordered_partition = ordered_partition
        self.numbers_of_gray_edges = numbers_of_gray_edges
        self.deletions = deletions
        self.duplications = duplications
        self.surrounding = surrounding

class multibreak:
    def __init__(self, affected = [], ordered_partition = None, numbers_of_gray_edges = None, proportion_chromo = None,   deletions = None, duplications = None, surrounding = []):
        self.affected = affected
        self.proportion_chromo = proportion_chromo
        self.ordered_partition = ordered_partition
        self.numbers_of_gray_edges = numbers_of_gray_edges
        self.deletions = deletions
        self.duplications = duplications
        self.surrounding = surrounding

    