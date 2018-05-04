# -*- coding: utf-8 -*-

from __future__ import absolute_import

import os
import itertools
import sys
from pickle import dump
import networkx as nx
from . import fetch_reactions


def _create_graph_with_internal_reaction(organismsdata):
    """
    This function creates a NetworkX DiGraph object which consists of
    reactions and metabolites happening inside the organisms in a community.
    This makes use of the reaction information i.e., irreversible and
    reversible, which is obtained from another script fetch_reactions.

    Parameters
    ----------
    organismsdata : dict
        Dictionary containing the reaction information about organisms

    Returns
    -------
    G : NetworkX DiGraph Object
        Bipartite graph consisting of internal reactions in organisms
    """
    G = nx.DiGraph()
    for modelname in organismsdata:
        G.add_nodes_from(organismsdata[modelname]
                         ['irreversible_rxn_no'], bipartite=1)
        G.add_nodes_from(organismsdata[modelname]
                         ['reversible_rxn_no'], bipartite=1)
        G.add_nodes_from(organismsdata[modelname]
                         ['reversible_back_rxn_no'], bipartite=1)
        irrev_lhs_nodes = list(set(
            [item for sublist in organismsdata[modelname]['irreversible_lhs_nodes']
             for item in sublist]))
        irrev_rhs_nodes = list(set(
            [item for sublist in organismsdata[modelname]['irreversible_rhs_nodes']
             for item in sublist]))
        rev_lhs_nodes = list(set(
            [item for sublist in organismsdata[modelname]['reversible_lhs_nodes']
             for item in sublist]))
        rev_rhs_nodes = list(set(
            [item for sublist in organismsdata[modelname]['reversible_rhs_nodes']
             for item in sublist]))
        G.add_nodes_from(irrev_lhs_nodes, bipartite=0)
        G.add_nodes_from(irrev_rhs_nodes, bipartite=0)
        G.add_nodes_from(rev_lhs_nodes, bipartite=0)
        G.add_nodes_from(rev_rhs_nodes, bipartite=0)
        for irrevidx in range(len(organismsdata[modelname]['irreversible_rxn_no'])):
            for lhsmetidx in range(len(organismsdata[modelname]['irreversible_lhs_nodes'][irrevidx])):
                G.add_edges_from([(organismsdata[modelname]['irreversible_lhs_nodes'][irrevidx]
                                   [lhsmetidx],
                                   organismsdata[modelname]['irreversible_rxn_no'][irrevidx])])
            for rhsmetidx in range(len(organismsdata[modelname]['irreversible_rhs_nodes'][irrevidx])):
                G.add_edges_from([(organismsdata[modelname]['irreversible_rxn_no'][irrevidx],
                                   organismsdata[modelname]['irreversible_rhs_nodes'][irrevidx][rhsmetidx])])
        for revidx in range(len(organismsdata[modelname]['reversible_rxn_no'])):
            for lhsmetidxrev in range(len(organismsdata[modelname]['reversible_lhs_nodes'][revidx])):
                G.add_edges_from([(organismsdata[modelname]['reversible_lhs_nodes'][revidx]
                                   [lhsmetidxrev],
                                   organismsdata[modelname]['reversible_rxn_no'][revidx])])
                G.add_edges_from([(organismsdata[modelname]['reversible_back_rxn_no'][revidx],
                                   organismsdata[modelname]['reversible_lhs_nodes'][revidx][lhsmetidxrev])])
            for rhsmetidxrev in range(len(organismsdata[modelname]['reversible_rhs_nodes'][revidx])):
                G.add_edges_from([(organismsdata[modelname]['reversible_rxn_no'][revidx],
                                   organismsdata[modelname]['reversible_rhs_nodes'][revidx][rhsmetidxrev])])
                G.add_edges_from([(organismsdata[modelname]['reversible_rhs_nodes'][revidx]
                                   [rhsmetidxrev],
                                   organismsdata[modelname]['reversible_back_rxn_no'][revidx])])
    return G


def _create_graph_with_exchange_reactions(G, orgs, namemap):
    """
    This function first identifies the common exchange metabolites
    and the non-common exchange metabolites and adds them to the
    DiGraph object generated above.

    Parameters
    ----------
    G : NetworkX DiGraph Object
        Bipartite graph of reaction network from organisms
    orgs : dict
        Dictionary consisting of irreversible, reversible and exchange
        reactions pertaining to the organisms. If more than one organism
        is used, this dictionary consists of information about all the
        organisms.
    namemap : dict
        Dictionary mapping the adhoc reaction names to reaction names in
        the model

    Returns
    -------
    G : NetworkX DiGraph Object
        Bipartite graph consisting of internal and exchange reactions in organisms
    namemap : dict
        Dictionary mapping the adhoc exchange reaction names to reaction names in
        the model
    """
    metabolite_exchanged = []
    for orgnames in orgs:
        exc_met = orgs[orgnames]['exchange_metab_nodes']
        metabolite_exchanged.append(exc_met)
    # Common exchange metabolites in different organisms
    common_exchange_metabolite = list(
        set.intersection(*list(map(set, metabolite_exchanged))))
    common_exchange_metabolite.sort()
    #  Adding the common exchange metabolites to the graph
    for orgnames in orgs:
        renamed_exc_met = [orgnames + ' ' +
                           comexcmet for comexcmet in common_exchange_metabolite]
        number_exc_met = list(range(0, len(common_exchange_metabolite)))
        mod_exc_rxn_number = ['Org_%s ER' %
                              orgnames + str(num + 1) for num in number_exc_met]
        mod_exc_rev_rxn_number = ['Org_%s ERR' %
                                  orgnames + str(num + 1) for num in number_exc_met]
        G.add_nodes_from(mod_exc_rxn_number, bipartite=1)
        G.add_nodes_from(mod_exc_rev_rxn_number, bipartite=1)
        G.add_nodes_from(common_exchange_metabolite, bipartite=0)
        G.add_nodes_from(renamed_exc_met, bipartite=0)
        for k in range(len(renamed_exc_met)):
            namemap[mod_exc_rxn_number[k]] = common_exchange_metabolite[k]
            namemap[mod_exc_rev_rxn_number[k]] = common_exchange_metabolite[k]
            G.add_edges_from([(renamed_exc_met[k], mod_exc_rxn_number[k])])
            G.add_edges_from(
                [(mod_exc_rxn_number[k], common_exchange_metabolite[k])])
            G.add_edges_from(
                [(common_exchange_metabolite[k], mod_exc_rev_rxn_number[k])])
            G.add_edges_from([(mod_exc_rev_rxn_number[k], renamed_exc_met[k])])
    #  Adding the non common exchange metabolites to the graph
    for orgnames in orgs:
        metitems = orgs[orgnames]['exchange_metab_nodes']
        non_common_exc_met = list(
            set(metitems) - set(common_exchange_metabolite))
        non_common_exc_met.sort()
        renamed_non_common_exc_met = [
            orgnames + ' ' + s for s in non_common_exc_met]
        number_non_common_exc_met = list(range(0, len(non_common_exc_met)))
        mod_non_common_exc_rxn_number = [
            'Org_%s NCER' % orgnames + str(num + 1) for num in number_non_common_exc_met]
        mod_non_common_exc_rev_rxn_number = [
            'Org_%s NCERR' % orgnames + str(num + 1) for num in number_non_common_exc_met]
        G.add_nodes_from(mod_non_common_exc_rxn_number, bipartite=1)
        G.add_nodes_from(mod_non_common_exc_rev_rxn_number, bipartite=1)
        G.add_nodes_from(non_common_exc_met, bipartite=0)
        G.add_nodes_from(renamed_non_common_exc_met, bipartite=0)
        for k in range(len(renamed_non_common_exc_met)):
            namemap[mod_non_common_exc_rxn_number[k]] = non_common_exc_met[k]
            namemap[mod_non_common_exc_rev_rxn_number[k]
                   ] = non_common_exc_met[k]
            G.add_edges_from(
                [(renamed_non_common_exc_met[k], mod_non_common_exc_rxn_number[k])])
            G.add_edges_from(
                [(mod_non_common_exc_rxn_number[k], non_common_exc_met[k])])
            G.add_edges_from(
                [(non_common_exc_met[k], mod_non_common_exc_rev_rxn_number[k])])
            G.add_edges_from(
                [(mod_non_common_exc_rev_rxn_number[k], renamed_non_common_exc_met[k])])
    return G, namemap


def create_graph(path_name_with_models, no_of_orgs):
    """
    This function creates bipartite graph of the organisms based on the
    path provided and the number of organsisms. For instance, if a folder
    has 3 model files, and the number of organisms is 2, 3 (3C2) different
    bipartite graphs are created. The graph objects and the dictionary
    are saved as gpickle and pickle files respectively.

    Parameters
    ----------
    path_name_with_models : str
        Absolute path name of the folder containing the models.
    no_of_orgs : int
        Number of organisms to be used for creating the DiGraph.

    Returns
    -------
    H : NetworkX DiGraph Object
        Bipartite graph consisting of internal and exchange reactions in organisms
    full_name_map : dict
        Dictionary mapping the adhoc reaction names to reaction names in
        the model
    """

    organisms_reaction_data, partial_name_map = \
        fetch_reactions.segregate_reactions_from_models(path_name_with_models)
    if organisms_reaction_data:
        organisms_names = list(organisms_reaction_data.keys())
        all_possible_combis = list(itertools.combinations(
            list(range(len(organisms_names))), int(no_of_orgs)))
        if all_possible_combis:
            for ncom in range(len(all_possible_combis)):
                file_name = ''
                current_combination = {}
                for numincom in range(len(all_possible_combis[ncom])):
                    current_combination[organisms_names[all_possible_combis[ncom][numincom]]] = \
                        organisms_reaction_data[organisms_names[all_possible_combis[ncom][numincom]]]
                    file_name = file_name + \
                        organisms_names[all_possible_combis[ncom]
                                        [numincom]] + '_'
                H = _create_graph_with_internal_reaction(current_combination)
                H, full_name_map = _create_graph_with_exchange_reactions(
                    H, current_combination, partial_name_map)              
                print('Number of edges in graph', len(H.edges()))
                print('Number of nodes in graph', len(H.nodes()))
                if os.access(path_name_with_models, os.W_OK):
                    with open(file_name + 'namemap' + '.pickle', 'wb') as filetodump:
                        dump(full_name_map, filetodump)
                    nx.write_gpickle(H, file_name + '.gpickle')
                    print('Graph and namemap saved for file(s) in', path_name_with_models)
            sys.path.append(path_name_with_models)
        else:
            print(
                'Number of organisms for creating a consortium graph is more than the models given')
            print('Program will now exit')
            sys.exit()
    else:
        print("Cannot create graph")
        sys.exit()
    return H, full_name_map
