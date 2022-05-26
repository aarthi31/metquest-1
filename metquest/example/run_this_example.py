# -*- coding: utf-8 -*-
from os.path import abspath, dirname, join
from os import listdir
from metquest import construct_graph
from metquest import pathway_assembler
from metquest import execute_metquest


metquest_directory = abspath(join(dirname(abspath(__file__)), ".."))
metquest_location = abspath(join(metquest_directory, ".."))
data_dir = join(metquest_directory, "example", "data", "")


def run_this_example():
    """
    This function runs the example with E. coli iJO1366 model with the
    seed, source and target metabolite input provided.

    Parameters
    ----------
    None

    Returns
    -------
    None

    """
    number_of_xml = 1
    G, namemap = construct_graph.create_graph(data_dir, number_of_xml)
    for files in listdir(data_dir):
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
        folder_to_create = 'foo'
        for currenttarmet in targetmetabolites:  # multiple target mets
            for cutoff in cutoff_list:  # multiple cutoffs
                pathway_table, cyclic_pathways, scope = \
                                pathway_assembler.find_pathways(G, seed_metabolites, int(cutoff))
                assert currenttarmet in pathway_table
                assert len(pathway_table['iJO1366 pyr_c'][15]) == 806
                number_of_pathways = []
                for plen in pathway_table['iJO1366 pyr_c']:
                    number_of_pathways.append(len(pathway_table['iJO1366 pyr_c'][plen]))
                assert sum(number_of_pathways) == 4787
                assert len(scope) == 885
                assert currenttarmet in cyclic_pathways
                execute_metquest.print_summary(scope, currenttarmet, pathway_table, cutoff,
                                               cyclic_pathways, namemap, source_metabolites,
                                               seed_metabolites, number_of_xml, G)


if __name__ == '__main__':
    run_this_example()
