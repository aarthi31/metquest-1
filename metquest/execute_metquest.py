# -*- coding: utf-8 -*-
from __future__ import absolute_import
import os
import sys
from collections import Counter
from itertools import combinations
from metquest.pathway_assembler import find_pathways
from metquest.construct_graph import create_graph


def write_output_to_file(pathway_table, currenttarmet, cutoff, cyclic_pathways,
                         folder_to_create, namemap, source_metabolites, G):
    """
    This function writes the pathways of sizes less than or equal to the
    cutoff from source to the target and seed metabolites to target.
    This function also writes cyclic pathways of sizes less than or equal to
    cutoff from the source to target.

    Parameters
    ----------
    pathway_table : dict
        Dictionary of dictionary containing the pathways of different sizes
        identified for every metabolite. This will have only the acyclic/
        branched pathways.
    currenttarmet : str
        Current target metabolite
    cutoff : int
        Maximum pathway length cutoff
    cyclic_pathways : dict
        Dictionary of dictionary containing cyclic pathways of different sizes
        identified for every metabolite.
    folder_to_create : str
        Name of the folder where results have to be written
    namemap : dict
        Dictionary mapping the adhoc reaction names to reaction names in
        the model
    source_metabolites : list
        List of source metabolites
    G : NetworkX DiGraph Object
        Bipartite graph of the metabolic network


    Returns
    -------
    None
    """
    pred = G.predecessors
    succ = G.successors
    all_pathways_count = []
    path_count = []
    cyclic_pathway_count = []
    if currenttarmet in pathway_table:
	    for plen in pathway_table[currenttarmet]:
	        all_pathways_count.append(
	            len(pathway_table[currenttarmet][plen]))
	        if plen <= int(cutoff):
	            path_count.append(
	                len(pathway_table[currenttarmet][plen]))
	    if currenttarmet in cyclic_pathways:
	        for plen in cyclic_pathways[currenttarmet]:
	            if plen <= int(cutoff):
	                cyclic_pathway_count.append(
	                    len(cyclic_pathways[currenttarmet][plen]))
	        cycfname = folder_to_create + 'cyclic_pathways_' + currenttarmet.replace(' ', '') + \
	            '_' + 'leq_plen_' + str(cutoff) + '.txt'
	        pathnumcount = 0
	        with open(cycfname, 'w') as filetowrite:
	            print('\nWriting cyclic pathways to a file')
	            for idx in cyclic_pathways[currenttarmet]:
	                if idx <= int(cutoff):
	                    filetowrite.write('Path length ' +
	                                      str(idx) + '\n')
	                for items in cyclic_pathways[currenttarmet][idx]:
	                    pathnumcount += 1
	                    filetowrite.write(
	                        str(pathnumcount) + '\n')
	                    for entities in list(items):
	                        filetowrite.write(
	                            namemap[entities] + '\t' + ' + '.join(pred(entities)) +
	                            '->' + ' + '.join(succ(entities)))
	                        filetowrite.write('\n')
	                    filetowrite.write('--------------------\n')

	    if int(cutoff) in pathway_table[currenttarmet]:
	        for sourcemets in source_metabolites:
	            seedfname = folder_to_create + 'branched_pathways_from_seed_' + \
	                    currenttarmet.replace(' ', '') + \
	                    '_' + 'leq_plen_' + str(cutoff) + '.txt'
	            pathnumcount = 0
	            with open(seedfname, 'w') as filetowrite:
	                print('Writing branched pathways (from seed) to a file')
	                for idx in pathway_table[currenttarmet]:
	                    if idx <= int(cutoff):
	                        filetowrite.write('Path length ' +
	                                          str(idx) + '\n')
	                        for items in pathway_table[currenttarmet][idx]:
	                            pathnumcount += 1
	                            filetowrite.write(
	                                str(pathnumcount) + '\n')
	                            for entities in list(items):
	                                filetowrite.write(
	                                    namemap[entities] + '\t' +
	                                    ' + '.join(pred(entities))
	                                    + '->' + ' + '.join(succ(entities)))
	                                filetowrite.write('\n')
	                            filetowrite.write(
	                                '--------------------\n')
	                        filetowrite.write('--------------------\n')
	            only_source_to_target = []

	            for sourcemets in source_metabolites:
	                for idx in pathway_table[currenttarmet]:
	                    if idx <= int(cutoff):
	                        for items in pathway_table[currenttarmet][idx]:
	                            if set(succ(sourcemets)).intersection(items):
	                                only_source_to_target.append(
	                                    list(items))
	                if only_source_to_target:
	                    sourcefname = folder_to_create + 'branched_pathways_from_source_' + \
	                                currenttarmet.replace(' ', '') + \
	                                '_' + 'leq_plen_' + str(cutoff) + '.txt'

	                    with open(sourcefname, 'w') as filetowrite:
	                        print('Writing branched pathways (from source) to a file \n')
	                        for currentidx, listentries in enumerate(only_source_to_target):
	                            filetowrite.write(str(currentidx+1) + '\n')
	                            filetowrite.write('Path length ' +
	                                              str(len(listentries)) + '\n')
	                            for entities in listentries:
	                                filetowrite.write(
	                                    namemap[entities] + '\t' +
	                                    ' + '.join(pred(entities))
	                                    + '->' + ' + '.join(succ(entities)))
	                                filetowrite.write('\n')
	                            filetowrite.write('--------------------\n')


def find_pathways_starting_from_source(source_metabolites, pathway_table, currenttarmet, cutoff, G):
    """
    This function finds all pathways starting from the source metabolites

    Parameters
    ----------
    source_metabolites : list
        List of source metabolites
    pathway_table : dict
        Dictionary of dictionary containing the pathways of different sizes
        identified for every metabolite. This will have only the acyclic/
        branched pathways.
    currenttarmet : str
        Current target metabolite
    cutoff : int
        Maximum pathway length cutoff
    G : NetworkX DiGraph Object
        Bipartite graph of the metabolic network

    Returns
    -------
    most_different_paths : dict
        For the given source metabolite, a combination of two most different pathways
        based on minimum Jaccard value is returned.
    only_source_to_target : list
        list of list containing all pathways starting from source metabolite

    """
    succ = G.successors
    most_different_paths = {}
    only_source_to_target = []
    for sourcemets in source_metabolites:
        for idx in pathway_table[currenttarmet]:
            if idx <= int(cutoff):
                for items in pathway_table[currenttarmet][idx]:
                    if set(succ(sourcemets)).intersection(items):
                        only_source_to_target.append(
                            list(items))
        if len(only_source_to_target) > 1: 
            # Sometimes there can be only one pathway producing target
            # To find most different paths from source
            j_value, rxn_comb = find_jaccard_between_paths(
                only_source_to_target)
            min_j_index = j_value.index(min(j_value))
            most_different_paths[sourcemets] = rxn_comb[min_j_index]
    return most_different_paths, only_source_to_target



def print_summary(scope, currenttarmet, pathway_table, cutoff, cyclic_pathways, namemap,
                  source_metabolites, folder_to_create, seed_metabolites,
                  number_of_xml, G):
    """
    This function prints the results summary obtained from the pathways, i.e.,
    1. Number of metabolites in scope
    2. Target metabolite
    3. Pathway size cutoff
    4. Number of all branched pathways found from seed
    5. Number of all branched pathways from seed whose size <= Pathway size cutoff
    6. Minimum number of steps to produce target metabolite
    7. Number of branched pathways from source whose size <= Pathway size cutoff
    8. Target metabolite can be produced using cyclic pathway
    9. Number of cyclic pathways whose size <= Pathway size cutoff
    10. One of the combination of most different pathways producing target metabolite
    11. Important reactions based on the frequency of occurrences

    Parameters
    ----------
    scope : set
        Set of metabolites that can be produced from the given set of
        seed metabolites
    currenttarmet : str
        Current target metabolite
    pathway_table : dict
        Dictionary of dictionary containing the pathways of different sizes
        identified for every metabolite. This will have only the acyclic/
        branched pathways.
    cutoff : int
        Maximum pathway length cutoff
    cyclic_pathways : dict
        Dictionary of dictionary containing cyclic pathways of different sizes
        identified for every metabolite.
    namemap : dict
        Dictionary mapping the adhoc reaction names to reaction names in
        the model
    source_metabolites : list
        List of source metabolites
    folder_to_create : str
        Name of the folder where results have to be written
    seed_metabolites : set
        Set of seed metabolites including the source
    number_of_xml : int
        Number of xml files in the folder
    G : NetworkX DiGraph Object
        Bipartite graph of the metabolic network

    Returns
    -------
    None

    """
    pred = G.predecessors
    succ = G.successors
    all_pathways_count = []
    path_count = []
    cyclic_pathway_count = []
    print('\n')
    print('---------------')
    print('Summary')
    print('---------------')
    print('Number of metabolites in scope : ', len(scope))
    print('Target metabolite : ', currenttarmet)
    print('Pathway size cutoff : ', str(cutoff))
    if currenttarmet in pathway_table:
        for plen in pathway_table[currenttarmet]:
            all_pathways_count.append(
                len(pathway_table[currenttarmet][plen]))
            if plen <= int(cutoff):
                path_count.append(
                    len(pathway_table[currenttarmet][plen]))
        print('Number of all branched pathways found from seed', ':', sum(all_pathways_count))
        print('Number of all branched pathways from seed whose size <=', cutoff, ':', sum(path_count))
        minsteps = min(pathway_table[currenttarmet])
        print('Minimum number of steps to produce ',
              currenttarmet, ' : ', int(minsteps))
        # Finding essential reactions and exchange metabolites
        all_reactions_involved = []
        for plen in pathway_table[currenttarmet]:
            for pathways in pathway_table[currenttarmet][plen]:
                for reactions in pathways:
                    all_reactions_involved.append(reactions)
        exchange_candidates_inverted_dict = find_pathways_involving_exchange_mets(number_of_xml, pathway_table, currenttarmet,
                                              seed_metabolites, namemap, G)
        most_different_paths, only_source_to_target = find_pathways_starting_from_source(source_metabolites, pathway_table,
                                                                                         currenttarmet, cutoff, G)
        # Two most different paths
        if most_different_paths:
            print('Number of branched pathways from source whose size <=',
                  cutoff, ':', len(only_source_to_target))

        if currenttarmet in cyclic_pathways:
            print(currenttarmet, 'can be produced using cyclic pathway')
            for plen in cyclic_pathways[currenttarmet]:
                if plen <= int(cutoff):
                    cyclic_pathway_count.append(
                        len(cyclic_pathways[currenttarmet][plen]))
            print('Number of cyclic pathways whose size <=',
                  cutoff, ':', sum(cyclic_pathway_count))
        else:
            print(currenttarmet, 'cannot be produced using cyclic pathway')
        if most_different_paths:
            print('\n')
            print('One of the combination of most different pathways producing target metabolite')
            print('Note - There can be other combinations that can be found')
            print('For finding all the combinations, please use the function find_jaccard_between_paths')
            for sourcemets in most_different_paths:
                counterpaths = 0
                for pathway in most_different_paths[sourcemets]:
                    counterpaths += 1
                    print('Pathway', str(counterpaths))
                    pathway.sort()
                    for rxns in pathway:
                        print(namemap[rxns], ' + '.join(list(pred(rxns))), '-->',
                              ' + '.join(list(succ(rxns))))
                    print('\n')
        else:
            print('No/only one pathway starting from source')
            print('Two most different paths from source : None')
        find_important_reactions(all_reactions_involved, currenttarmet, seed_metabolites, namemap, G)
    else:
        print(currenttarmet, ': Target could not be found.')
        print('Consider changing the cut-off or the seed metabolite set')


def find_important_reactions(all_reactions_involved, currenttarmet, seed_metabolites, namemap, G):
    """
    This function determines the important reactions based on the pathways
    generated for the target metabolite.

    Parameters
    ----------
    all_reactions_involved : list
        list of all reactions found in all the pathways from source to target
    currenttarmet : str
        Current target metabolite
    seed_metabolites : set
        Set of seed metabolites including the source
    namemap : dict
        Dictionary mapping the adhoc reaction names to reaction names in
        the model
    G : NetworkX DiGraph Object
        Bipartite graph of the metabolic network


    Returns
    -------
    None

    Notes
    -----
    We define important reactions as those reactions which occur in almost all
    the pathways producing the target metabolite (apart from the reactions that
    are involved in the production of target metabolite and the uptake of seed
    metabolite)

    """
    pred = G.predecessors
    important_candidates = Counter(all_reactions_involved)
    important_reactions_inverted_dict = {}
    for keys, values in important_candidates.items():
        important_reactions_inverted_dict[values] = \
            important_reactions_inverted_dict.get(values, [])
        important_reactions_inverted_dict[values].append(keys)
    number_of_occurrences = list(important_reactions_inverted_dict.keys())
    number_of_occurrences.sort()
    top_candidate = number_of_occurrences[:-6:-1]

    important_reactions = []
    for numrepeats in top_candidate:
        for rxns in important_reactions_inverted_dict[numrepeats]:
            if rxns not in pred(currenttarmet):
                if not set(pred(rxns)).issubset(seed_metabolites):
                    important_reactions.append(rxns)
    important_reactions_model_names = []
    for rxns in important_reactions:
        important_reactions_model_names.append(namemap[rxns])
    important_reactions_model_names.sort()
    if important_reactions:
        print('Important reactions based on the frequency of occurrence are')
        print('\n'.join(important_reactions_model_names))
    else:
        print('All reactions pertain to uptake of seed metabolite/ production of target metabolite')


def find_pathways_involving_exchange_mets(number_of_xml, pathway_table, currenttarmet,
                                          seed_metabolites, namemap, G):
    """
    This function identifies the pathways producing the target metabolites,
    which involve exchange metabolites. This function prints output only when
    a community of organisms is considered, i.e., when more than one metabolic
    network is used.

    Parameters
    ----------
    number_of_xml : int
        Number of xml files in the folder
    pathway_table : dict
        Dictionary of dictionary containing the pathways of different sizes
        identified for every metabolite. This will have only the acyclic/
        branched pathways.
    currenttarmet : str
        Current target metabolite
    seed_metabolites : set
        Set of seed metabolites including the source
    namemap : dict
        Dictionary mapping the adhoc reaction names to reaction names in
        the model
    G : NetworkX DiGraph Object
        Bipartite graph of the metabolic network

    Returns
    -------
    None

    """
    pred = G.predecessors
    succ = G.successors
    exchange_reactions = []
    if number_of_xml > 1:
        for plen in pathway_table[currenttarmet]:
            for pathways in pathway_table[currenttarmet][plen]:
                for reactions in pathways:
                    if not set(pred(reactions)).issubset(seed_metabolites):
                        #  'ER' is an adhoc reaction name assigned to exchange
                        #  reactions in the models.
                        if 'ER' in reactions:
                            exchange_reactions.append(reactions)
        exchange_candidates = Counter(exchange_reactions)
        exchange_candidates_inverted_dict = {}
        for keys, values in exchange_candidates.items():
            exchange_candidates_inverted_dict[values] = \
                exchange_candidates_inverted_dict.get(values, [])
            exchange_candidates_inverted_dict[values].append(keys)
        number_of_occurrences = list(exchange_candidates_inverted_dict.keys())
        number_of_occurrences.sort()
        # [:-6:-1] # Taking the top reaction
        top_candidate = number_of_occurrences[:-10:-1]
        exchange_reactions = []
        for numrepeats in top_candidate:
            for rxns in exchange_candidates_inverted_dict[numrepeats]:
                if rxns not in pred(currenttarmet):
                    if not set(pred(rxns)).issubset(seed_metabolites):
                        exchange_reactions.append(rxns)
        if exchange_reactions:
            print('Exchange reactions are')
            for rxns in exchange_reactions:
                print(namemap[rxns], list(pred(rxns)), list(succ(rxns)))
        else:
            print('No metabolite exchanged')
	return exchange_candidates_inverted_dict


def find_jaccard_between_paths(only_source_to_target):
    """
    This function determines the jaccard values between the pathways generated
    from the source to the target.

    Parameters
    ----------
    only_source_to_target : list
        list of lists consisting of all pathways producing the target
        metabolite from the source

    Returns
    -------
    jaccard_values : list
        list of all jaccard values (float) for all the pathway combinations
    path_combinations : list
        list of all pathway combinations corresponding to the jaccard values

    Notes
    -----
    Jaccard value J = (set(A).intersection(set(B)))/(set(A).union(set(B)))
    J = 1 indicates two sets are the same
    J = 0 indicates two sets are different

    """

    jaccard_values = []
    path_combinations = []
    for reactionlists in combinations(only_source_to_target, 2):
        j_value = len(set(reactionlists[0]).intersection(
            set(reactionlists[1])))/len(set(reactionlists[0]).union(set(reactionlists[1])))
        jaccard_values.append(j_value)
        path_combinations.append(reactionlists)
    return jaccard_values, path_combinations


def execute_all_codes():
    """
    This function executes all the codes including constructing graphs and executing metquest.
    
    Parameters
    ----------
    None
    
    Returns
    -------
    None
    
    """
    try:
        inputfoldername = sys.argv[1]
    except:
        inputfoldername = input('Enter folder name with all files \n')

    if '~' in inputfoldername:
        inputfoldername = os.path.expanduser(inputfoldername)
    list_of_files = os.listdir(inputfoldername)
    for foldernames in list_of_files:
        # To go through only folders
        if os.path.isdir(os.path.join(inputfoldername, foldernames)):
            current_evaluation_folder = os.path.join(
                inputfoldername, foldernames)
            number_of_files_in_current_folder = os.listdir(
                current_evaluation_folder)
            number_of_xml = len(
                [filenames for filenames in number_of_files_in_current_folder if '.xml' in filenames])
            print('Currently evaluating files in', foldernames)
            print('Number of networks', number_of_xml)
            G, namemap = create_graph(
                current_evaluation_folder, number_of_xml)
            for files in os.listdir('.'):
                if files.endswith('.txt'):
                    if files.startswith('seed'):
                        with open(files, 'r') as seedfile:
                            seedmetslist = seedfile.read().splitlines()
                        seed_metabolites = set(seedmetslist)
                    elif files.startswith('source'):
                        with open(files, 'r') as sourcefile:
                            source_metabolites = sourcefile.read().splitlines()
                    elif files.startswith('target'):
                        with open(files, 'r') as targetfile:
                            # Target metabolites can be multiple values
                            targetmetabolites = targetfile.read().splitlines()
                    elif files.startswith('cutoff'):
                        with open(files, 'r') as cutofffile:
                            cutoff_list = cutofffile.read().splitlines()  # Cutoff can be multiple values
            for mets in source_metabolites:
                seed_metabolites.add(mets)
            metfoundingraph = True
            for metabs in seed_metabolites:
                if metabs not in G:
                    print(metabs, 'not in G. MetQuest will not be executed')
                    print('Please check the metabolite names')
                    metfoundingraph = False
                    break
            for metabs in targetmetabolites:
                if metabs not in G:
                    print(metabs, 'not in G. MetQuest will not be executed')
                    print('Please check the metabolite names')
                    metfoundingraph = False
                    break
            if metfoundingraph:
                folder_to_create = 'Results/'
                if not os.path.exists(folder_to_create):
                    os.makedirs(folder_to_create)
                for currenttarmet in targetmetabolites:  # multiple target mets
                    for cutoff in cutoff_list:  # multiple cutoffs
                        pathway_table, cyclic_pathways, scope = find_pathways(
                            G, seed_metabolites, int(cutoff))
                        print_summary(scope, currenttarmet, pathway_table, cutoff, cyclic_pathways,
                                      namemap, source_metabolites, folder_to_create,
                                      seed_metabolites, number_of_xml, G)
                        write_output_to_file(pathway_table, currenttarmet, cutoff,
                                             cyclic_pathways, folder_to_create,
                                             namemap, source_metabolites, G)
                print('\n')
            os.chdir('../')
