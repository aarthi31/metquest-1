# -*- coding: utf-8 -*-

from collections import deque, defaultdict


def forward_pass(graph_object, seedmets):
    """
    This function carries out the Guided Breadth First Search on a directed
    bipartite graph starting from the entries in seed metabolite set.

    Parameters
    ----------
    graph_object : NetworkX DiGraph Object
        Bipartite graph of the metabolic network

    seedmets : set
        Set of seed metabolites including the source

    Returns
    -------
    lower_bound_metabolite : defaultdict
        Minimum number of steps required to reach a metabolite
    status_dict : defaultdict
        Dictionary pertaining to the status of every reaction - whether it
        has been visited or not
    scope : set
        Set of metabolites that can be produced from the given set of
        seed metabolites

    Notes
    -----
    Starting with the set of seed metabolites S, the algorithm first finds
    all the reactions from the set R, whose precursor metabolites are in S.
    Such reactions are marked visited and added to the visited reaction set.
    Metabolites produced by these reactions are checked. The reactions where
    these metabolites participate are then checked for the presence of all its
    predecessors and are added to the queue. This traversal continues in a
    breadth-first manner and stops when there are no further reactions to
    be visited.
    """
    pred = graph_object.predecessors
    succ = graph_object.successors
    seed_metabolite_set = seedmets.copy()
    lower_bound_metabolite = defaultdict(list)
    lower_bound_reaction = defaultdict(list)
    # Defaultdict is used simply because to avoid initialisations
    status_dict = defaultdict(str)
    # Using a deque since deques have O(1) speed for appendleft() and popleft()
    # while lists have O(n) performance for inserting and popping.
    queue = deque([])
    # All seed metabolites are always present, hence require 0 steps
    for seedmetabs in seed_metabolite_set:
        lower_bound_metabolite[seedmetabs].append(0)
    stage = 1
    scope = seed_metabolite_set.copy()
    starting_rxn_node = []
    # First stage where starting_rxn_node list contains all the reactions
    # which require only the seed metabolites as input
    for starting_met_nodes in seed_metabolite_set:
        # Essential when analysing mutiple networks with same seed metabolite
        # set, although would be redundant in case of single network
        if starting_met_nodes in graph_object:
            for startingrxns in succ(starting_met_nodes):
                if set(pred(startingrxns)).issubset(seed_metabolite_set):
                    if startingrxns not in starting_rxn_node:
                        starting_rxn_node.append(startingrxns)
                    for metsprod in succ(startingrxns):
                        scope.add(metsprod)
                        if stage not in lower_bound_metabolite[metsprod]:
                            lower_bound_metabolite[metsprod].append(stage)
                    if stage not in lower_bound_reaction[startingrxns]:
                        lower_bound_reaction[startingrxns].append(stage)
    for rxn in starting_rxn_node:
        for metabs in succ(rxn):
            for nextrxn in succ(metabs):
                if set(pred(nextrxn)).issubset(scope):
                    if nextrxn not in queue:
                        queue.append(nextrxn)
        status_dict[rxn] = 'V'
    while queue:
        stage += 1
        for parentrxn in list(queue):
            if status_dict[parentrxn] == '':
                if stage not in lower_bound_reaction[parentrxn]:
                    lower_bound_reaction[parentrxn].append(stage)
                for mets in succ(parentrxn):
                    scope.add(mets)
                    if stage not in lower_bound_metabolite[mets]:
                        lower_bound_metabolite[mets].append(stage)
                    for progeny in succ(mets):
                        if set(pred(progeny)).issubset(scope):

                            if status_dict[progeny] != 'V':
                                if progeny not in queue:
                                    queue.append(progeny)
                status_dict[parentrxn] = 'V'
            elif status_dict[parentrxn] == 'V':
                for mets in succ(parentrxn):
                    if stage not in lower_bound_metabolite[mets]:
                        lower_bound_metabolite[mets].append(stage)
            queue.popleft()
    return lower_bound_metabolite, status_dict, scope
