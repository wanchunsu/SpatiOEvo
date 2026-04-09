import os
import os.path as osp
import json
from dictionary_tools import load_json, dump_json
import numpy as np

### label nodes and edges of interface-specific subgraphs

def get_dnds_cat(tp_well_defined_dnds_for_restype_curr_prot, node):
	if tp_well_defined_dnds_for_restype_curr_prot[node] > 1.0:
		return 'pos'
	elif tp_well_defined_dnds_for_restype_curr_prot[node] < 0.1:
		return 'neg'
	else:
		return 'other'

def categorize_residues_as_pos_neg_other(restype, tp_well_defined_dnds, res_distances_dict, outfi):

	# Go through <restype>_res_distances.json (protein: node1-node2: distance) and extract the nodes
	# Then store whether the node is +vely/-vely selected or other (using targ_prot_res_well_defined_dnds_vals.json # restype: protein: respos: dndsval
	
	# Output: {tp1: {tp1node1: pos/neg/other, tp2node2: pos/neg/other, ...}, tp2: {}, ...}
	nodes_dnds_categories = {}
	if 'specific' in restype:
		tp_well_defined_dnds_for_restype = tp_well_defined_dnds[restype.split('_')[0]]
	elif 'mim' in restype:
		tp_well_defined_dnds_for_restype = tp_well_defined_dnds['mimicry']
	else:	
		tp_well_defined_dnds_for_restype = tp_well_defined_dnds[restype]
	for tp in res_distances_dict:
		tp_well_defined_dnds_for_restype_curr_prot = tp_well_defined_dnds_for_restype[tp]
		if tp not in nodes_dnds_categories:
			nodes_dnds_categories[tp] = {}

		for node1_node2 in res_distances_dict[tp]:
			node1 = node1_node2.split('-')[0]
			node2 = node1_node2.split('-')[1]

			if node1 not in nodes_dnds_categories[tp]:
				nodes_dnds_categories[tp][node1] = get_dnds_cat(tp_well_defined_dnds_for_restype_curr_prot, node1)

			if node2 not in nodes_dnds_categories[tp]:
				nodes_dnds_categories[tp][node2] = get_dnds_cat(tp_well_defined_dnds_for_restype_curr_prot, node2)
	# print(nodes_dnds_categories)
	dump_json(nodes_dnds_categories, outfi)
	return nodes_dnds_categories


def categorize_residues_as_pos_neg_other_refined1(restype, tp_well_defined_dnds, res_distances_dict, res_distances_dict2, outfi):
	'''
	Also stores residue contacts that span two interfaces (i.e. res1 in exo-specific, res2 in mimicked), by bringing in all exo and all endo data
	Go through <restype>_res_distances.json (protein: node1-node2: distance) and extract the nodes
	Then store whether the node is +vely/-vely selected or other (using targ_prot_res_well_defined_dnds_vals.json # restype: protein: respos: dndsval
	
	Output: {tp1: {tp1node1: pos/neg/other, tp2node2: pos/neg/other, ...}, tp2: {}, ...}
	'''

	nodes_dnds_categories = {}
	if 'specific' in restype:
		tp_well_defined_dnds_for_restype = tp_well_defined_dnds[restype.split('_')[0]]
	elif 'mim' in restype:
		tp_well_defined_dnds_for_restype = tp_well_defined_dnds['mimicry']
	else:	
		tp_well_defined_dnds_for_restype = tp_well_defined_dnds[restype]
	for tp in res_distances_dict:
		tp_well_defined_dnds_for_restype_curr_prot = tp_well_defined_dnds_for_restype[tp]
		if tp not in nodes_dnds_categories:
			nodes_dnds_categories[tp] = {}

		for node1_node2 in res_distances_dict[tp]:
			node1 = node1_node2.split('-')[0]
			node2 = node1_node2.split('-')[1]

			if node1 not in nodes_dnds_categories[tp]:
				nodes_dnds_categories[tp][node1] = get_dnds_cat(tp_well_defined_dnds_for_restype_curr_prot, node1)

			if node2 not in nodes_dnds_categories[tp]:
				nodes_dnds_categories[tp][node2] = get_dnds_cat(tp_well_defined_dnds_for_restype_curr_prot, node2)

	# Now we're going to check for residues that are only in contact with residues in other interfaces (they will only be found in the all_exo or all_endo res distances data)
	for tp in res_distances_dict2: 
		if tp not in tp_well_defined_dnds_for_restype: continue
		tp_well_defined_dnds_for_restype_curr_prot = tp_well_defined_dnds_for_restype[tp]
		
		for node1_node2 in res_distances_dict2[tp]:
			node1 = node1_node2.split('-')[0]
			node2 = node1_node2.split('-')[1]
			if node1 in tp_well_defined_dnds_for_restype_curr_prot:
				if tp not in nodes_dnds_categories:
					nodes_dnds_categories[tp] = {}
				if node1 not in nodes_dnds_categories[tp]:
					nodes_dnds_categories[tp][node1] = get_dnds_cat(tp_well_defined_dnds_for_restype_curr_prot, node1)
			if node2 in tp_well_defined_dnds_for_restype_curr_prot:
				if tp not in nodes_dnds_categories:
					nodes_dnds_categories[tp] = {}
				if node2 not in nodes_dnds_categories[tp]:
					nodes_dnds_categories[tp][node2] = get_dnds_cat(tp_well_defined_dnds_for_restype_curr_prot, node2)
	# print(nodes_dnds_categories)
	dump_json(nodes_dnds_categories, outfi)
	return nodes_dnds_categories


def categorize_residues_as_pos_neg_other_refined2(restype, tp_well_defined_dnds, res_distances_dict, res_distances_dict2, res_distances_dict3, outfi):
	'''
	Also stores residue contacts that span two interfaces (i.e. res1 in exo-specific, res2 in mimicked)	by bringing in all exo and all endo data
	Go through <restype>_res_distances.json (protein: node1-node2: distance) and extract the nodes
	Then store whether the node is +vely/-vely selected or other (using targ_prot_res_well_defined_dnds_vals.json # restype: protein: respos: dndsval
	
	Output: {tp1: {tp1node1: pos/neg/other, tp2node2: pos/neg/other, ...}, tp2: {}, ...}
	'''
	nodes_dnds_categories = {}
	if 'specific' in restype:
		tp_well_defined_dnds_for_restype = tp_well_defined_dnds[restype.split('_')[0]]
	elif 'mim' in restype:
		tp_well_defined_dnds_for_restype = tp_well_defined_dnds['mimicry']
	else:	
		tp_well_defined_dnds_for_restype = tp_well_defined_dnds[restype]
	for tp in res_distances_dict:
		tp_well_defined_dnds_for_restype_curr_prot = tp_well_defined_dnds_for_restype[tp]
		if tp not in nodes_dnds_categories:
			nodes_dnds_categories[tp] = {}

		for node1_node2 in res_distances_dict[tp]:
			node1 = node1_node2.split('-')[0]
			node2 = node1_node2.split('-')[1]

			if node1 not in nodes_dnds_categories[tp]:
				nodes_dnds_categories[tp][node1] = get_dnds_cat(tp_well_defined_dnds_for_restype_curr_prot, node1)

			if node2 not in nodes_dnds_categories[tp]:
				nodes_dnds_categories[tp][node2] = get_dnds_cat(tp_well_defined_dnds_for_restype_curr_prot, node2)

	# Now we're going to check for residues that are only in contact with residues in other interfaces (they will only be found in the all_exo or all_endo res distances data)
	for tp in res_distances_dict2: 
		if tp not in tp_well_defined_dnds_for_restype: continue
		tp_well_defined_dnds_for_restype_curr_prot = tp_well_defined_dnds_for_restype[tp]
		
		for node1_node2 in res_distances_dict2[tp]:
			node1 = node1_node2.split('-')[0]
			node2 = node1_node2.split('-')[1]
			if node1 in tp_well_defined_dnds_for_restype_curr_prot:
				if tp not in nodes_dnds_categories:
					nodes_dnds_categories[tp] = {}
				if node1 not in nodes_dnds_categories[tp]:
					nodes_dnds_categories[tp][node1] = get_dnds_cat(tp_well_defined_dnds_for_restype_curr_prot, node1)
			if node2 in tp_well_defined_dnds_for_restype_curr_prot:
				if tp not in nodes_dnds_categories:
					nodes_dnds_categories[tp] = {}
				if node2 not in nodes_dnds_categories[tp]:
					nodes_dnds_categories[tp][node2] = get_dnds_cat(tp_well_defined_dnds_for_restype_curr_prot, node2)

	for tp in res_distances_dict3: 
		if tp not in tp_well_defined_dnds_for_restype: continue
		tp_well_defined_dnds_for_restype_curr_prot = tp_well_defined_dnds_for_restype[tp]
		
		for node1_node2 in res_distances_dict3[tp]:
			node1 = node1_node2.split('-')[0]
			node2 = node1_node2.split('-')[1]
			if node1 in tp_well_defined_dnds_for_restype_curr_prot:
				if tp not in nodes_dnds_categories:
					nodes_dnds_categories[tp] = {}
				if node1 not in nodes_dnds_categories[tp]:
					nodes_dnds_categories[tp][node1] = get_dnds_cat(tp_well_defined_dnds_for_restype_curr_prot, node1)
			if node2 in tp_well_defined_dnds_for_restype_curr_prot:
				if tp not in nodes_dnds_categories:
					nodes_dnds_categories[tp] = {}
				if node2 not in nodes_dnds_categories[tp]:
					nodes_dnds_categories[tp][node2] = get_dnds_cat(tp_well_defined_dnds_for_restype_curr_prot, node2)

	# print(nodes_dnds_categories)
	dump_json(nodes_dnds_categories, outfi)
	return nodes_dnds_categories


def categorize_edges_in_contact_graph(res_distances_dict, nodes_dnds_categories, outfi):
	# Go through <restype>_res_distances.json (protein: node1-node2: distance) and extract edges (i.e. distances < 5Angstroms)
	# Check the labels of the two nodes in contact

	# Output: node1-node2: label

	edges_dnds_categories = {}
	for tp in res_distances_dict:
		if tp not in edges_dnds_categories:
			edges_dnds_categories[tp] = {}
		for node1_node2 in res_distances_dict[tp]:
			node1 = node1_node2.split('-')[0]
			node2 = node1_node2.split('-')[1]
			if (node1 + '-' + node2 not in edges_dnds_categories[tp]) and (node2 + '-' + node1 not in edges_dnds_categories[tp]):
				# store edge if the
				if res_distances_dict[tp][node1_node2] < 5.0: # nearest neighbour cutoff: 5 Angstroms
					if nodes_dnds_categories[tp][node1] == 'pos' and nodes_dnds_categories[tp][node2] == 'pos': # two +vely selected residues are nearest neighbours
						edges_dnds_categories[tp][node1 + '-' + node2] = 'pos'
					elif nodes_dnds_categories[tp][node1] == 'neg' and nodes_dnds_categories[tp][node2] == 'neg': # two -vely selected residues are nearest neighbours
						edges_dnds_categories[tp][node1 + '-' + node2] = 'neg'
					elif nodes_dnds_categories[tp][node1] == 'pos' and nodes_dnds_categories[tp][node2] == 'neg': # one +vely one -vely selected residue are nns
						edges_dnds_categories[tp][node1 + '-' + node2] = 'pos_neg'
					elif nodes_dnds_categories[tp][node1] == 'neg' and nodes_dnds_categories[tp][node2] == 'pos': # one +vely one -vely selected residue are nns
						edges_dnds_categories[tp][node1 + '-' + node2] = 'pos_neg'
					else:
						edges_dnds_categories[tp][node1 + '-' + node2] = 'other'
	# print(edges_dnds_categories)
	dump_json(edges_dnds_categories, outfi)
	return edges_dnds_categories

def categorize_edges_in_contact_graph2(res_distances_dict, nodes_dnds_categories, outfi):
	# Go through <restype>_res_distances.json (protein: node1-node2: distance) and extract edges (i.e. distances < 5Angstroms)
	# Check the labels of the two nodes in contact

	# Output: node1-node2: label

	edges_dnds_categories = {}
	for tp in res_distances_dict:
		if tp not in edges_dnds_categories:
			edges_dnds_categories[tp] = {}
		for node1_node2 in res_distances_dict[tp]:
			node1 = node1_node2.split('-')[0]
			node2 = node1_node2.split('-')[1]
			if (node1 + '-' + node2 not in edges_dnds_categories[tp]) and (node2 + '-' + node1 not in edges_dnds_categories[tp]):
				# store edge if the
				if res_distances_dict[tp][node1_node2] < 5.0: # nearest neighbour cutoff: 5 Angstroms
					if nodes_dnds_categories[tp][node1] == 'pos' and nodes_dnds_categories[tp][node2] == 'pos': # two +vely selected residues are nearest neighbours
						edges_dnds_categories[tp][node1 + '-' + node2] = 'pos'
					elif nodes_dnds_categories[tp][node1] == 'neg' and nodes_dnds_categories[tp][node2] == 'neg': # two -vely selected residues are nearest neighbours
						edges_dnds_categories[tp][node1 + '-' + node2] = 'neg'
					elif nodes_dnds_categories[tp][node1] == 'pos' and nodes_dnds_categories[tp][node2] == 'neg': # one +vely one -vely selected residue are nns
						edges_dnds_categories[tp][node1 + '-' + node2] = 'pos_neg'
					elif nodes_dnds_categories[tp][node1] == 'neg' and nodes_dnds_categories[tp][node2] == 'pos': # one +vely one -vely selected residue are nns
						edges_dnds_categories[tp][node1 + '-' + node2] = 'pos_neg'
					elif (nodes_dnds_categories[tp][node1] == 'pos' and nodes_dnds_categories[tp][node2] == 'other') or (nodes_dnds_categories[tp][node1] == 'other' and nodes_dnds_categories[tp][node2] == 'pos'): # one +vely selected residue one other are nns
						edges_dnds_categories[tp][node1 + '-' + node2] = 'pos_other'
					elif (nodes_dnds_categories[tp][node1] == 'neg' and nodes_dnds_categories[tp][node2] == 'other') or (nodes_dnds_categories[tp][node1] == 'other' and nodes_dnds_categories[tp][node2] == 'neg'): # one -vely selected residue one other are nns
						edges_dnds_categories[tp][node1 + '-' + node2] = 'neg_other'
					else: # two other are nns
						edges_dnds_categories[tp][node1 + '-' + node2] = 'other'
	# print(edges_dnds_categories)
	dump_json(edges_dnds_categories, outfi)
	return edges_dnds_categories


def main():
	script_dir = osp.dirname(__file__)

	tp_well_defined_dnds_fi = osp.join(script_dir, '..', 'data', 'dnds', 'targ_prot_res_well_defined_dnds_vals.json')
	tp_well_defined_dnds = load_json(tp_well_defined_dnds_fi)
	
	outdir = osp.join(script_dir, '..', 'data', 'contact_graphs')

	if not osp.exists(outdir):
		os.makedirs(outdir)
	for restype in ['exo_specific', 'mimicked', 'endo_specific', 'all_exo', 'all_endo']:
		res_distances_dict_fi = osp.join(script_dir, '..', 'data', 'dnds', 'distances_btwn_residues_in_diff_dnds_categories', restype + '_res_distances.json')
		res_distances_dict = load_json(res_distances_dict_fi)
	
		print(f'\n##### Restype: {restype} #####')
		out_nodes_anno_fi = osp.join(outdir, restype + '_nodes_anno.json')
		out_edges_anno_fi = osp.join(outdir, restype + '_edges_anno.json')
		out_edges_anno_fi2 = osp.join(outdir, restype + '_edges_anno_split_other.json')

		print('\tMaking node annotations dict . . .')
		
		if restype == 'all_exo' or restype == 'all_endo':
			nodes_dnds_categories = categorize_residues_as_pos_neg_other(restype, tp_well_defined_dnds, res_distances_dict, out_nodes_anno_fi)
		
		elif restype == 'exo_specific':
			res_distances_dict_2 = load_json(osp.join(script_dir, '..', 'data', 'dnds', 'dnds_for_residue_categories_species_tree', 'distances_btwn_residues_in_diff_dnds_categories', 'all_exo_res_distances.json'))
			nodes_dnds_categories = categorize_residues_as_pos_neg_other_refined1(restype, tp_well_defined_dnds, res_distances_dict, res_distances_dict_2, out_nodes_anno_fi)
		
		elif restype == 'endo_specific':
			res_distances_dict_2 = load_json(osp.join(script_dir, '..', 'data', 'dnds', 'dnds_for_residue_categories_species_tree', 'distances_btwn_residues_in_diff_dnds_categories', 'all_endo_res_distances.json'))
			nodes_dnds_categories = categorize_residues_as_pos_neg_other_refined1(restype, tp_well_defined_dnds, res_distances_dict, res_distances_dict_2, out_nodes_anno_fi)
		
		elif restype == 'mimicked':
			res_distances_dict_2 = load_json(osp.join(script_dir, '..', 'data', 'dnds', 'dnds_for_residue_categories_species_tree', 'distances_btwn_residues_in_diff_dnds_categories', 'all_exo_res_distances.json'))
			res_distances_dict_3 = load_json(osp.join(script_dir, '..', 'data', 'dnds', 'dnds_for_residue_categories_species_tree', 'distances_btwn_residues_in_diff_dnds_categories', 'all_endo_res_distances.json'))
			nodes_dnds_categories = categorize_residues_as_pos_neg_other_refined2(restype, tp_well_defined_dnds, res_distances_dict, res_distances_dict_2, res_distances_dict_3, out_nodes_anno_fi)
			
		print('\tMaking edge annotations dict . . .')
		edges_dnds_categories = categorize_edges_in_contact_graph(res_distances_dict, nodes_dnds_categories, out_edges_anno_fi)
		

		# Part2: where we also categorize 'other' into 'pos_other', 'neg_other', and 'other'
		print('\tMaking edge annotations dict . . .')
		edges_dnds_categories = categorize_edges_in_contact_graph2(res_distances_dict, nodes_dnds_categories, out_edges_anno_fi2)
		

	print('DONE! :)')


if __name__ == '__main__':
	main()
	
