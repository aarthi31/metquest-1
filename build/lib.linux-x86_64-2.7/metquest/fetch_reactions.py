# -*- coding: utf-8 -*-

from __future__ import absolute_import

import os
import glob
import cobra
from .get_reaction_types import find_different_reaction_types


def segregate_reactions_from_models(path_name):
    """
    This function gets the data pertaining to the reactions and the
    metabolites from the models of multiple organisms.
    This requires as input the pathname where the '.xml' files are located.
    From this path, this function reads all the files using the functions
    in the COBRA toolbox and generates the stoichiometric model for these
    SBML models.

    Parameters
    ----------
    path_name : str
        full path name where the model files are

    Returns
    -------
    all_organisms_info : dict
        Dictionary of all model data (reaction information about all the
                                      organisms)
    namemap : dict
        Dictionary mapping the adhoc reaction names to reaction names in
        the model

    """
    all_organisms_info = {}
    namemap = {}
    os.chdir(path_name)
    file_names = glob.glob('*.xml')
    if not file_names:
        print("There are no .xml files. Please check the path")
    print("Filenames", file_names)
    for model_names in file_names:
        model = cobra.io.read_sbml_model(model_names)
        stoi = cobra.util.array.create_stoichiometric_matrix(model)
        if model.id:
            current_model_name = model.id
        else:
            print("Model ID not found; using file name instead")
            current_model_name = model_names.split('.')[0]
        current_organisms_info = {current_model_name: {'exchange_metab_nodes': [],
                                                       'irreversible_lhs_nodes': [],
                                                       'irreversible_rhs_nodes': [],
                                                       'reversible_rhs_nodes': [],
                                                       'reversible_lhs_nodes': [],
                                                       'irreversible_rxn_no': [],
                                                       'reversible_rxn_no': [],
                                                       'total_nodes': [],
                                                       'model_rxns': [],
                                                       'metabolites': [],
                                                       'exch_rxn_name': [],
                                                       'irrev_rxn_name': [],
                                                       'rev_rxn_name': []}}
        rxns_in_model = []
        mets_in_model = []
        for metab in model.metabolites:
            mets_in_model.append(metab.id)
        for reac in model.reactions:
            rxns_in_model.append(reac.id)
        stoi_matrix = stoi.T
        exchange_nodes, irrev_lhs_nodes, irrev_rhs_nodes, rev_lhs_nodes, rev_rhs_nodes, \
            exc_name, irrev_rxn_name, rev_rxn_name = find_different_reaction_types(
                stoi_matrix, model, current_model_name)
        current_organisms_info[current_model_name][
            'exchange_metab_nodes'] = exchange_nodes
        current_organisms_info[current_model_name][
            'irreversible_lhs_nodes'] = irrev_lhs_nodes
        current_organisms_info[current_model_name][
            'irreversible_rhs_nodes'] = irrev_rhs_nodes
        current_organisms_info[current_model_name][
            'reversible_lhs_nodes'] = rev_lhs_nodes
        current_organisms_info[current_model_name][
            'reversible_rhs_nodes'] = rev_rhs_nodes
        current_organisms_info[current_model_name]['exch_rxn_name'] = exc_name
        current_organisms_info[current_model_name][
            'irrev_rxn_name'] = irrev_rxn_name
        current_organisms_info[current_model_name][
            'rev_rxn_name'] = rev_rxn_name

        irrev_rxn_number = []
        for num in range(len(irrev_lhs_nodes)):
            modified_name_irrev = 'Org_%s IR' % current_model_name + str(num + 1)
            irrev_rxn_number.append(modified_name_irrev)
            namemap[modified_name_irrev] = irrev_rxn_name[num]

        rev_rxn_number = []
        for num in range(len(rev_lhs_nodes)):
            modified_name_rev = 'Org_%s RR' % current_model_name + str(num + 1)
            rev_rxn_number.append(modified_name_rev)
            namemap[modified_name_rev] = rev_rxn_name[num]

        rev_back_rxn_number = []
        for num in range(len(rev_lhs_nodes)):
            modified_name_back_rev = 'Org_%s RevBR' % current_model_name + \
                str(num + 1)
            rev_back_rxn_number.append(modified_name_back_rev)
            namemap[modified_name_back_rev] = rev_rxn_name[num]

        current_organisms_info[current_model_name][
            'reversible_rxn_no'] = rev_rxn_number
        current_organisms_info[current_model_name][
            'irreversible_rxn_no'] = irrev_rxn_number
        current_organisms_info[current_model_name]['total_nodes'] = len(
            exchange_nodes) + len(irrev_lhs_nodes) + len(rev_lhs_nodes)
        current_organisms_info[current_model_name]['model_rxns'] = rxns_in_model
        current_organisms_info[current_model_name][
            'reversible_back_rxn_no'] = rev_back_rxn_number
        current_organisms_info[current_model_name]['metabolites'] = mets_in_model
        all_organisms_info.update(current_organisms_info)
    return all_organisms_info, namemap
