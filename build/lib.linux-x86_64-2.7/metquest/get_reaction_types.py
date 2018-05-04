# -*- coding: utf-8 -*-

import numpy as np


def find_different_reaction_types(stoi_matrix, model, current_model_name):
    """
    This function finds the exchange, irreversible and the reversible reactions
    from the model.

    Parameters
    ----------
    stoi_matrix : numpy array
        full path name where the model files are
    model : COBRA model object
        COBRA model object created from SBML models
    current_model_name : str
        Name which is to be prefixed against every
        reaction/metabolite (to differentiate the entries in multiple organisms,
        when a community model is built)
    Returns
    -------
    exchange_met_ids : list
        Metabolite identifiers of exchange metabolites
    irrev_lhs_nodes : list
        Metabolite identifiers of reactants of irreversible reactions
    irrev_rhs_nodes : list
        Metabolite identifiers of products of irreversible reactions
    rev_lhs_nodes : list
        Metabolite identifiers of reactants of reversible reactions
    rev_rhs_nodes : list
        Metabolite identifiers of products of reversible reactions
    exchange_rxn_ids : list
        Reaction identifers of exchange reactions
    irrev_rxn_ids : list
        Reaction identifiers of irreversible reactions
    rev_rxn_ids : list
        Reaction identifiers of reversible reactions

    Notes
    -----
    Consider a model consisting of four reactions:
        a <=>    -1000<v<1000
        <=> b     -1000<v<1000
        a + e <=> f -1000<v<1000
        b + d -> c 0<v<1000,
    where v denotes the flux through every reaction.
    For this, the stoichiometric matrix S is of the form
        S = [[-1,0,0,0,0,0], [0,1,0,0,0,0], [-1,0,0,0,-1,1], [0,-1,1,-1,0,0]]
    Firstly, from S, the exchange reactions are found by considering the
    reactions which have a single entry in every list. In this process, the
    bulk metabolites are also considered.
    This list of exchange reactions is removed from the master reaction list,
    denoted as internal_reactions. From this, based on the lower and the upper
    bounds, irreversible and reversible reactions are separated. The LHS and
    the RHS metabolite of these reactions are returned as separate list.
    """

    xdim = np.shape(stoi_matrix)
    reactants_of_reaction = []
    total_metabolites_in_reaction = []
    products_of_reaction = []
    number_of_reactants_in_reaction = []
    total_number_of_metabs_in_reaction = []
    number_of_products_in_reaction = []
    exchange_reaction_idx = []
    reaction_identifiers = []
    reaction_in_model = []
    metabolite_identifiers = []
    for metab in model.metabolites:
        metabolite_identifiers.append(metab.id)
    for rxns in model.reactions:
        reaction_identifiers.append(rxns.id)
        reaction_in_model.append(rxns.reaction)
    for rxnidx in range(xdim[0]):
        reactants_of_reaction.append(np.where(stoi_matrix[rxnidx] == -1))
        total_metabolites_in_reaction.append(np.where(stoi_matrix[rxnidx] != 0))
        products_of_reaction.append(np.where(stoi_matrix[rxnidx] == 1))
        number_of_reactants_in_reaction.append(len(reactants_of_reaction[rxnidx][0]))
        total_number_of_metabs_in_reaction.append(len(total_metabolites_in_reaction[rxnidx][0]))
        number_of_products_in_reaction.append(len(products_of_reaction[rxnidx][0]))

        # Case 1 - Presence of bulk metabolites in the medium

        if reaction_in_model[rxnidx][-1] == 'b':  # Assuming the bulk metabolites end in 'b'
            if (number_of_reactants_in_reaction[rxnidx] == 1 and
                    number_of_products_in_reaction[rxnidx] == 1):
                exchange_reaction_idx.append(rxnidx)
        # Case 2 - Presence of exchange metabolites
        elif (number_of_reactants_in_reaction[rxnidx] == 1 and
              total_number_of_metabs_in_reaction[rxnidx] == 1):
            exchange_reaction_idx.append(rxnidx)
        elif (number_of_products_in_reaction[rxnidx] == 1 and
              total_number_of_metabs_in_reaction[rxnidx] == 1):
            exchange_reaction_idx.append(rxnidx)
    exchange_met_ids = []
    exchange_met_index = []
    exchange_rxn_ids = []
    for excentry in exchange_reaction_idx:
        exchange_rxn_ids.append(reaction_identifiers[excentry])
        if reaction_in_model[excentry][-1] == 'b':
            exchange_met_ids.append(metabolite_identifiers
                                    [np.nonzero(stoi_matrix[excentry])[0][0]])
        else:
            exchange_met_index.append(
                np.nonzero(stoi_matrix[excentry])[0].tolist()[0])
    if exchange_met_index:
        for metind in exchange_met_index:
            exchange_met_ids.append(metabolite_identifiers[metind])
    all_rxn_idx = list(range(len(reaction_in_model)))
    internal_rxns = list(set(all_rxn_idx) ^ set(exchange_reaction_idx))
    reversible_rxns = []
    irreversible_rxns = []
    rxns_lowerbound = []
    rxns_upperbound = []
    for rxns in model.reactions:
        rxns_lowerbound.append(rxns.lower_bound)
        rxns_upperbound.append(rxns.upper_bound)
    for idxint in internal_rxns:
        if rxns_lowerbound[idxint] < 0 and rxns_upperbound[idxint] >= 0:
            reversible_rxns.append(idxint)
        elif rxns_lowerbound[idxint] >= 0 and rxns_upperbound[idxint] >= 0:
            irreversible_rxns.append(idxint)
    #  Irreversible reaction nodes
    irrev_lhs_temporary = []
    irrev_rhs_temporary = []
    irrev_lhs_nodes = []
    irrev_rhs_nodes = []
    irrev_rxn_ids = []
    for irridx in irreversible_rxns:
        irrev_rxn_ids.append(reaction_identifiers[irridx])
        irrev_lhs_temporary.append(np.where(stoi_matrix[irridx] < 0)[0].tolist())
        irrev_rhs_temporary.append(np.where(stoi_matrix[irridx] > 0)[0].tolist())
    for lhsirridx in range(len(irrev_lhs_temporary)):
        temp_metab_list_lhs = []
        for met_idx_lhs in irrev_lhs_temporary[lhsirridx]:
            met_namech_lhs = "%s %s" % (current_model_name,
                                        metabolite_identifiers[met_idx_lhs])
            temp_metab_list_lhs.append(met_namech_lhs)
        irrev_lhs_nodes.append(temp_metab_list_lhs)
    for rhsirridx in range(len(irrev_rhs_temporary)):
        temp_metab_list_rhs = []
        for met_idx_rhs in irrev_rhs_temporary[rhsirridx]:
            met_namech_rhs = "%s %s" % (current_model_name,
                                        metabolite_identifiers[met_idx_rhs])
            temp_metab_list_rhs.append(met_namech_rhs)
        irrev_rhs_nodes.append(temp_metab_list_rhs)

    #  Reversible reaction nodes
    rev_lhs_temporary = []
    rev_rhs_temporary = []
    rev_lhs_nodes = []
    rev_rhs_nodes = []
    rev_rxn_ids = []
    for rridx in reversible_rxns:
        rev_rxn_ids.append(reaction_identifiers[rridx])
        rev_lhs_temporary.append(np.where(stoi_matrix[rridx] < 0)[0].tolist())
        rev_rhs_temporary.append(np.where(stoi_matrix[rridx] > 0)[0].tolist())
    for lhsrevidx in range(len(rev_lhs_temporary)):
        temp_metab_list_lhs_rev = []
        for met_idx_lhs in rev_lhs_temporary[lhsrevidx]:
            met_namech_lhs = "%s %s" % (current_model_name,
                                        metabolite_identifiers[met_idx_lhs])
            temp_metab_list_lhs_rev.append(met_namech_lhs)
        rev_lhs_nodes.append(temp_metab_list_lhs_rev)
    for rhsrevidx in range(len(rev_rhs_temporary)):
        temp_metab_list_rhs_rev = []
        for met_idx_rhs in rev_rhs_temporary[rhsrevidx]:
            met_namech_rhs = "%s %s" % (current_model_name,
                                        metabolite_identifiers[met_idx_rhs])
            temp_metab_list_rhs_rev.append(met_namech_rhs)
        rev_rhs_nodes.append(temp_metab_list_rhs_rev)
    return exchange_met_ids, irrev_lhs_nodes, \
        irrev_rhs_nodes, rev_lhs_nodes, rev_rhs_nodes, exchange_rxn_ids, \
        irrev_rxn_ids, rev_rxn_ids
