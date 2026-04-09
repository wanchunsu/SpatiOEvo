import os
import os.path as osp
import json
from dictionary_tools import load_json, dump_json
import numpy as np
import random
from collections import Counter
from scipy.stats import mannwhitneyu 

'''
Compute O/E pos-pos contact count ratios within and between interfaces
To get random contact graph (for calculating expected contacts), shuffle node-labels in exo-specific, endo-specific, and mimicked separately 
'''

def label_pos_pos_contacts(all_exo_or_all_endo_edges, exo_spec_or_endo_spec_nodes, mimicked_nodes, exo_or_endo):
	'''
	Categorize pos-pos contacts based on whether they are occurring between
	[if exo_or_endo == 'exo']
		1. two exo-specific residues
		2. two mimicked residues
		3. between an exo-specific and a mimicked residue
	OR
	[if exo_or_endo == 'endo']
		1. two endo-specific residues
		2. two mimicked residues
		3. between an endo-specific and a mimicked residue

	'''

	between_label = 'btwn_' + exo_or_endo + '_and_mimicked'
	dict_of_pos_pos_contacts = {exo_or_endo: [], 'mimicked': [], between_label: []}
	for tp in all_exo_or_all_endo_edges:
		
		for residue_pair in all_exo_or_all_endo_edges[tp]:
			
			if all_exo_or_all_endo_edges[tp][residue_pair] == 'pos':
				res1 = residue_pair.split('-')[0]
				res2 = residue_pair.split('-')[1]
				
				if tp not in exo_spec_or_endo_spec_nodes: # these two residues must be mimicked
					dict_of_pos_pos_contacts['mimicked'].append(tp + "_" + residue_pair)
				
				elif tp not in mimicked_nodes: # these two residues must be exo or endo
					dict_of_pos_pos_contacts[exo_or_endo].append(tp + "_" + residue_pair)
				
				else:
					if res1 in exo_spec_or_endo_spec_nodes[tp] and res2 in exo_spec_or_endo_spec_nodes[tp]: # the two residues are in exo or endo
						dict_of_pos_pos_contacts[exo_or_endo].append(tp + "_" + residue_pair)
					
					elif res1 in mimicked_nodes[tp] and res2 in mimicked_nodes[tp]:
						dict_of_pos_pos_contacts['mimicked'].append(tp + "_" + residue_pair)
					
					else:
						dict_of_pos_pos_contacts[between_label].append(tp + "_" + residue_pair)

	dict_of_num_pos_pos_contacts = {}
	total = 0
	
	for label in dict_of_pos_pos_contacts:
		dict_of_num_pos_pos_contacts[label] = len(dict_of_pos_pos_contacts[label])
		total += len(dict_of_pos_pos_contacts[label])

	dict_of_num_pos_pos_contacts['total'] = total

	# print(dict_of_num_pos_pos_contacts)
	return dict_of_pos_pos_contacts, dict_of_num_pos_pos_contacts


def label_contacts(all_exo_or_all_endo_edges, exo_spec_or_endo_spec_nodes, mimicked_nodes, exo_or_endo):
	'''
	Categorize all contacts based on whether they are occurring between
	[if exo_or_endo == 'exo']
		1. two exo-specific residues
		2. two mimicked residues
		3. between an exo-specific and a mimicked residue
	OR
	[if exo_or_endo == 'endo']
		1. two endo-specific residues
		2. two mimicked residues
		3. between an endo-specific and a mimicked residue

	'''

	between_label = 'btwn_' + exo_or_endo + '_and_mimicked'
	dict_of_contacts = {exo_or_endo: [], 'mimicked': [], between_label: []}
	for tp in all_exo_or_all_endo_edges:
		for residue_pair in all_exo_or_all_endo_edges[tp]:
			# if all_exo_or_all_endo_edges[tp][residue_pair] == 'pos':
			res1 = residue_pair.split('-')[0]
			res2 = residue_pair.split('-')[1]
			if tp not in exo_spec_or_endo_spec_nodes: # these two residues must be mimicked
				dict_of_contacts['mimicked'].append(tp + "_" + residue_pair)
			elif tp not in mimicked_nodes: # these two residues must be exo or endo
				dict_of_contacts[exo_or_endo].append(tp + "_" + residue_pair)
			else:
				if res1 in exo_spec_or_endo_spec_nodes[tp] and res2 in exo_spec_or_endo_spec_nodes[tp]: # the two residues are in exo or endo
					dict_of_contacts[exo_or_endo].append(tp + "_" + residue_pair)
				elif res1 in mimicked_nodes[tp] and res2 in mimicked_nodes[tp]:
					dict_of_contacts['mimicked'].append(tp + "_" + residue_pair)
				else:
					dict_of_contacts[between_label].append(tp + "_" + residue_pair)

	dict_of_num_contacts = {}
	total = 0
	for label in dict_of_contacts:
		dict_of_num_contacts[label] = len(dict_of_contacts[label])
		total += len(dict_of_contacts[label])

	dict_of_num_contacts['total'] = total

	print(dict_of_num_contacts)
	# return dict_of_contacts, dict_of_num_contacts


# def calc_stats_of_pos_in_interface(dict_of_pos_pos_contacts, exo_spec_or_endo_spec_nodes, mimicked_nodes, exo_or_endo):
# 	# dict_of_pos_pos_contacts = {exo_or_endo: [], 'mimicked': [], between_label: []}
# 	# for all pos residues in an interface, calculate the fraction that are 
# 	# 1. in contact with within interface pos residues
# 	# 2. in contact with other interface pos residues
# 	between_label = 'btwn_' + exo_or_endo + '_and_mimicked'

# 	dict_of_fractions_for_exo_or_endo_specific_nodes = {'within_pos_pos_contacts':[], 'btwn_pos_pos_contacts': []}
# 	dict_of_fractions_for_mimicked_nodes = {'within_pos_pos_contacts':[], 'btwn_pos_pos_contacts': []}

# 	for tp_respair in dict_of_pos_pos_contacts[exo_or_endo]:
# 		tp = tp_respair.split('_')[0]
# 		res1 = tp_respair.split('_')[1].split('-')[0]
# 		res2 = tp_respair.split('_')[1].split('-')[1]
# 		if tp + '_' + res1 not in dict_of_fractions_for_exo_or_endo_specific_nodes['within_pos_pos_contacts']:
# 			dict_of_fractions_for_exo_or_endo_specific_nodes['within_pos_pos_contacts'].append(tp + '_' + res1)
# 		if tp + '_' + res2 not in dict_of_fractions_for_exo_or_endo_specific_nodes['within_pos_pos_contacts']:
# 			dict_of_fractions_for_exo_or_endo_specific_nodes['within_pos_pos_contacts'].append(tp + '_' + res2)

# 	for tp_respair in dict_of_pos_pos_contacts['mimicked']:
# 		tp = tp_respair.split('_')[0]
# 		res1 = tp_respair.split('_')[1].split('-')[0]
# 		res2 = tp_respair.split('_')[1].split('-')[1]
# 		if tp + '_' + res1 not in dict_of_fractions_for_mimicked_nodes['within_pos_pos_contacts']:
# 			dict_of_fractions_for_mimicked_nodes['within_pos_pos_contacts'].append(tp + '_' + res1)
# 		if tp + '_' + res2 not in dict_of_fractions_for_mimicked_nodes['within_pos_pos_contacts']:
# 			dict_of_fractions_for_mimicked_nodes['within_pos_pos_contacts'].append(tp + '_' + res2)

# 	for tp_respair in dict_of_pos_pos_contacts[between_label]:
# 		tp = tp_respair.split('_')[0]
# 		res1 = tp_respair.split('_')[1].split('-')[0]
# 		res2 = tp_respair.split('_')[1].split('-')[1]
		
# 		if res1 in exo_spec_or_endo_spec_nodes[tp] and res2 in mimicked_nodes[tp]:
# 			if tp + '_' + res1 not in dict_of_fractions_for_exo_or_endo_specific_nodes['btwn_pos_pos_contacts']:
# 				dict_of_fractions_for_exo_or_endo_specific_nodes['btwn_pos_pos_contacts'].append(tp + '_' + res1)
# 			if tp + '_' + res2 not in dict_of_fractions_for_mimicked_nodes['btwn_pos_pos_contacts']:
# 				dict_of_fractions_for_mimicked_nodes['btwn_pos_pos_contacts'].append(tp + '_' + res2)
# 		elif res1 in mimicked_nodes[tp] and res2 in exo_spec_or_endo_spec_nodes[tp]:
# 			if tp + '_' + res1 not in dict_of_fractions_for_mimicked_nodes['btwn_pos_pos_contacts']:
# 				dict_of_fractions_for_mimicked_nodes['btwn_pos_pos_contacts'].append(tp + '_' + res1)
# 			if tp + '_' + res2 not in dict_of_fractions_for_exo_or_endo_specific_nodes['btwn_pos_pos_contacts']:
# 				dict_of_fractions_for_exo_or_endo_specific_nodes['btwn_pos_pos_contacts'].append(tp + '_' + res2)
# 		else:
# 			print('Error!') 
# 	print(len(dict_of_fractions_for_exo_or_endo_specific_nodes['within_pos_pos_contacts']), len(dict_of_fractions_for_exo_or_endo_specific_nodes['btwn_pos_pos_contacts']))
# 	print(len(dict_of_fractions_for_mimicked_nodes['within_pos_pos_contacts']), len(dict_of_fractions_for_mimicked_nodes['btwn_pos_pos_contacts']))		
	
# 	len_exo_or_endo_pos_nodes = 0
# 	len_mimicked_pos_nodes = 0
# 	for i in exo_spec_or_endo_spec_nodes:
# 		for res in exo_spec_or_endo_spec_nodes[i]:
# 			if exo_spec_or_endo_spec_nodes[i][res] == 'pos':
# 				len_exo_or_endo_pos_nodes+=1
# 	for i in mimicked_nodes:
# 		for res in mimicked_nodes[i]:
# 			if mimicked_nodes[i][res] == 'pos':
# 				len_mimicked_pos_nodes+=1
# 	print(len_exo_or_endo_pos_nodes, len_mimicked_pos_nodes)


def shuffle_node_annos(nodes_anno): 
	nodes_anno_shuffled_all_prots = {}

	for tp in nodes_anno:
		orig_nodes_anno = nodes_anno[tp]
		orig_annos_counts = Counter(orig_nodes_anno[r] for r in orig_nodes_anno) # count to verify 

		# shuffling values
		temp = list(orig_nodes_anno.values())
		random.shuffle(temp)
		
		# reassigning to keys
		shuffled_nodes_anno = dict(zip(orig_nodes_anno, temp))
		shuffled_annos_counts = Counter(orig_nodes_anno[r] for r in orig_nodes_anno) # count to verify 

		# print(orig_nodes_anno, shuffled_nodes_anno)
		if orig_annos_counts!=shuffled_annos_counts: # count to verify 
			print('False')

		nodes_anno_shuffled_all_prots[tp] = shuffled_nodes_anno
	return nodes_anno_shuffled_all_prots


def shuffle_nodes_repeatedly(nodes_anno): # do this for exo-specific, endo-specific, and mimicked separately!

	num_trials = 10000
	nodes_anno_shuffled_all_prots_repeatedly = {}
	# Get null distribution of expected nodes
	for t in range(num_trials): 
		print(f'\t[Trial {t}] Shuffling node annotations . . .')
		nodes_anno_shuffled_all_prots_repeatedly[t] = shuffle_node_annos(nodes_anno)

	return nodes_anno_shuffled_all_prots_repeatedly

# For each trial shuffle node labels for exo-specific, endo-specific, and mimicked separetly.
# Feed in relabeled nodes for exo-specific & mimicked to relabel edges for all exo
# Feed in relabeled nodes for endo-specific & mimicked to relabel edges for all endo
def relabel_edges_in_graph_based_on_shuffled_annos(exo_or_endo_nodes_anno_shuffled_all_prots, mimicked_nodes_anno_shuffled_all_prots, orig_edges_anno):
	new_edges_anno = {}
	num_edges_all_orig = 0
	nodes_anno_shuffled_all_prots = {} # dict to store combined annotations for nodes (exo-specific + mimicked OR endo-specific + mimicked) to relabel all exo OR all endo
	
	print('Relabeling edges in graph based on shuffled node labels . . .')

	# Populate with new node labels for exo- OR endo-specific
	for tp in exo_or_endo_nodes_anno_shuffled_all_prots:
		if tp not in nodes_anno_shuffled_all_prots:
			nodes_anno_shuffled_all_prots[tp] = {}
		for res in exo_or_endo_nodes_anno_shuffled_all_prots[tp]:
			nodes_anno_shuffled_all_prots[tp][res] = exo_or_endo_nodes_anno_shuffled_all_prots[tp][res]
	
	# Populate with new node labels for mimicked
	for tp in mimicked_nodes_anno_shuffled_all_prots:
		if tp not in nodes_anno_shuffled_all_prots:
			nodes_anno_shuffled_all_prots[tp] = {}
		for res in mimicked_nodes_anno_shuffled_all_prots[tp]:
			nodes_anno_shuffled_all_prots[tp][res] = mimicked_nodes_anno_shuffled_all_prots[tp][res]

	# Now relabel edges based on new shuffled node labels
	for tp in orig_edges_anno:
		new_edges_anno[tp] = {}
		for pair in orig_edges_anno[tp]:
			num_edges_all_orig+=1
			node1 = pair.split('-')[0]
			node2 = pair.split('-')[1]
			# print(nodes_anno_shuffled_all_prots[tp])
			# print(tp,pair)
			# print(nodes_anno_shuffled_all_prots[tp])
			if nodes_anno_shuffled_all_prots[tp][node1] == 'pos' and nodes_anno_shuffled_all_prots[tp][node2] == 'pos': # two +vely selected residues are nearest neighbours
				new_edges_anno[tp][node1 + '-' + node2] = 'pos'
			elif nodes_anno_shuffled_all_prots[tp][node1] == 'neg' and nodes_anno_shuffled_all_prots[tp][node2] == 'neg': # two -vely selected residues are nearest neighbours
				new_edges_anno[tp][node1 + '-' + node2] = 'neg'
			elif nodes_anno_shuffled_all_prots[tp][node1] == 'pos' and nodes_anno_shuffled_all_prots[tp][node2] == 'neg': # one +vely one -vely selected residue are nns
				new_edges_anno[tp][node1 + '-' + node2] = 'pos_neg'
			elif nodes_anno_shuffled_all_prots[tp][node1] == 'neg' and nodes_anno_shuffled_all_prots[tp][node2] == 'pos': # one +vely one -vely selected residue are nns
				new_edges_anno[tp][node1 + '-' + node2] = 'pos_neg'
			elif (nodes_anno_shuffled_all_prots[tp][node1] == 'pos' and nodes_anno_shuffled_all_prots[tp][node2] == 'other') or (nodes_anno_shuffled_all_prots[tp][node1] == 'other' and nodes_anno_shuffled_all_prots[tp][node2] == 'pos'): # one +vely selected residue one other are nns
				new_edges_anno[tp][node1 + '-' + node2] = 'pos_other'
			elif (nodes_anno_shuffled_all_prots[tp][node1] == 'neg' and nodes_anno_shuffled_all_prots[tp][node2] == 'other') or (nodes_anno_shuffled_all_prots[tp][node1] == 'other' and nodes_anno_shuffled_all_prots[tp][node2] == 'neg'): # one -vely selected residue one other are nns
				new_edges_anno[tp][node1 + '-' + node2] = 'neg_other'
			else:
				new_edges_anno[tp][node1 + '-' + node2] = 'other'
	num_edges_new = 0
	for tp in new_edges_anno:
		num_edges_new += len(new_edges_anno[tp])
	# print(num_edges_all_orig, num_edges_new) # verify that we did this correctly

	return new_edges_anno
		

def get_nodes_to_use(all_exo_or_endo_nodes, exo_or_endo_specific_nodes, mim_nodes):
	# go through all exo /all endo and sort out which ones are exo_specific, and which are mimicked
	new_exo_or_endo_specific_nodes = {}
	new_mim_nodes = {}
	for tp in all_exo_or_endo_nodes:
		for res in all_exo_or_endo_nodes[tp]:
			if tp in exo_or_endo_specific_nodes:
				if res in exo_or_endo_specific_nodes[tp]:
					if tp not in new_exo_or_endo_specific_nodes:
						new_exo_or_endo_specific_nodes[tp] = {}
					
					new_exo_or_endo_specific_nodes[tp][res] = exo_or_endo_specific_nodes[tp][res]

					if exo_or_endo_specific_nodes[tp][res] != all_exo_or_endo_nodes[tp][res]: #sanity check
						print('ERROR!!!')	

				elif res in mim_nodes[tp]:
					if tp not in new_mim_nodes:
						new_mim_nodes[tp] = {}
					
					new_mim_nodes[tp][res] = mim_nodes[tp][res]

					if new_mim_nodes[tp][res] != mim_nodes[tp][res]: #sanity check
						print('ERROR!!!')
				else:
					print(f'ERROR! This residue {tp}, {res} has no category...')
			else: # tp not in exo_or_endo_specific
				if res in mim_nodes[tp]:
					if tp not in new_mim_nodes:
						new_mim_nodes[tp] = {}
					
					new_mim_nodes[tp][res] = mim_nodes[tp][res]

					if new_mim_nodes[tp][res] != mim_nodes[tp][res]: #sanity check
						print('ERROR!!!')
				else:
					print(f'ERROR! This residue {tp}, {res} has no category...')

	return new_exo_or_endo_specific_nodes, new_mim_nodes



def check(new_edges, spec_nodes, mim_nodes): #sanity check!
	for tp in new_edges:
		if tp in spec_nodes and tp in mim_nodes:
			for respairs in new_edges[tp]:
				r1 = respairs.split('-')[0]
				r2 = respairs.split('-')[1]
				if r1 not in spec_nodes[tp] and r1 not in mim_nodes[tp]:
					print('UH OH!')
				if r2 not in spec_nodes[tp] and r2 not in mim_nodes[tp]:
					print('UH OH!')


def calc_o_over_e(dict_of_observed_edge_counts, dict_of_expected_edge_counts_for_each_trial, exo_or_endo):
	dict_of_counts = {'observed': {}, 'expected': {}}
	between_label = 'btwn_' + exo_or_endo + '_and_mimicked'
	types_of_pos_contacts = [exo_or_endo, 'mimicked', between_label]
	
	for tc in types_of_pos_contacts:
		dict_of_counts['observed'][tc] = dict_of_observed_edge_counts[tc]

		dict_of_counts['expected'][tc] = {}
	
		for trial in dict_of_expected_edge_counts_for_each_trial:
			# print(trial)
			# print(dict_of_expected_edge_counts_for_each_trial[trial])
			dict_of_counts['expected'][tc][trial] = dict_of_expected_edge_counts_for_each_trial[trial][tc]
	

	dict_of_o_over_e = {exo_or_endo: {}, 'mimicked': {}, between_label: {}}

	for tc in types_of_pos_contacts:
		for trial in dict_of_counts['expected'][tc]:
			o_over_e = dict_of_counts['observed'][tc]/dict_of_counts['expected'][tc][trial]
			dict_of_o_over_e[tc][trial] = o_over_e

	# print(exo_or_endo, np.mean(list(dict_of_o_over_e[exo_or_endo].values())))
	# print('mimicked', np.mean(list(dict_of_o_over_e['mimicked'].values())))
	# print(between_label, np.mean(list(dict_of_o_over_e[between_label].values())))

	return dict_of_o_over_e

def calc_pvals(comparison_types, dict_of_o_over_e):
	for ct in comparison_types:
		data1 = list(dict_of_o_over_e[ct[0]].values())
		data2 = list(dict_of_o_over_e[ct[1]].values())
		res = mannwhitneyu(data1, data2, alternative='two-sided')
		
		# Output the result
		print(f'{ct:}')
		print(f"\tMann Whitney U statistic: {res.statistic}")
		print(f"\tP-value: {res.pvalue:.20f}")


def main():
	script_dir = osp.dirname(__file__)
	data_dir = osp.join(script_dir, '..', 'data', 'contact_graphs')
	outdir = osp.join(script_dir, '..', 'data', 'contact_graphs', 'interrel', 'oe_of_pos_pos_within_and_btwn_interfaces')
	os.makedirs(outdir, exist_ok=True)

	# Load nodes and edges
	all_exo_nodes = load_json(osp.join(data_dir, 'all_exo_nodes_anno.json'))
	all_endo_nodes = load_json(osp.join(data_dir, 'all_endo_nodes_anno.json'))

	exo_specific_nodes = load_json(osp.join(data_dir, 'exo_specific_nodes_anno.json'))
	endo_specific_nodes = load_json(osp.join(data_dir, 'endo_specific_nodes_anno.json'))

	mimicked_nodes = load_json(osp.join(data_dir, 'mimicked_nodes_anno.json'))

	all_exo_edges = load_json(osp.join(data_dir, 'all_exo_edges_anno_split_other.json'))
	all_endo_edges = load_json(osp.join(data_dir, 'all_endo_edges_anno_split_other.json'))

	
	# Split all_exo 
	new_exo_specific_nodes, new_mim_exo_nodes = get_nodes_to_use(all_exo_nodes, exo_specific_nodes, mimicked_nodes)
	check(all_exo_edges, new_exo_specific_nodes, new_mim_exo_nodes)

	# Split all_endo
	new_endo_specific_nodes, new_mim_endo_nodes = get_nodes_to_use(all_endo_nodes, endo_specific_nodes, mimicked_nodes)
	check(all_endo_edges, new_endo_specific_nodes, new_mim_endo_nodes)


	# Number of pos-pos contacts within all_exo/all_endo that occur between the same (exo/endo) or two different interfaces (exo/endo with mim)
	all_exo_dict_of_pos_pos_contacts, all_exo_dict_of_num_pos_pos_contacts = label_pos_pos_contacts(all_exo_edges, new_exo_specific_nodes, new_mim_exo_nodes, 'exo')

	all_endo_dict_of_pos_pos_contacts, all_endo_dict_of_num_pos_pos_contacts = label_pos_pos_contacts(all_endo_edges, new_endo_specific_nodes, new_mim_endo_nodes, 'endo')

	# print(all_exo_dict_of_num_pos_pos_contacts, all_endo_dict_of_num_pos_pos_contacts)
	

	'''Now generate random contact graphs, where we shuffle node labels repeatedly (separately for each interface)
	
		1. For each trial, relabel the edges for all exo/all endo after separately shuffling exo-specific, endo-specific, and mimicked. 
		
		2. Then label each pos-pos contact as either within exo-specific/within endo-specific OR within mimicked, OR between exo-specific and mimicked/endo-specific and mimicked
	'''
	
	## All exo ##
	
	all_exo_dict_of_num_pos_pos_contacts_random_trials = {}
	all_exo_outfi_random_trials = osp.join(outdir, 'all_exo_random_trials_pos_pos_counts.json')
	
	if not osp.exists(all_exo_outfi_random_trials):
		print('\nRestype: All exo')	
		
		# Shuffle nodes 
		exo_specific_nodes_anno_shuffled_repeatedly = shuffle_nodes_repeatedly(new_exo_specific_nodes)
		mimicked_nodes_anno_shuffled_repeatedly = shuffle_nodes_repeatedly(new_mim_exo_nodes)
		
		for t in mimicked_nodes_anno_shuffled_repeatedly:
			#1. go through the trials and get a new relabeled set of all_exo/all_endo edges for each trial
			all_exo_new_edges_anno_for_trial = relabel_edges_in_graph_based_on_shuffled_annos(exo_specific_nodes_anno_shuffled_repeatedly[t], mimicked_nodes_anno_shuffled_repeatedly[t], all_exo_edges)
			
			#2. Get number of pos-pos contacts in each label type for each randomization trial
			_, all_exo_dict_of_num_pos_pos_contacts_random_trials[t] = label_pos_pos_contacts(all_exo_new_edges_anno_for_trial, exo_specific_nodes, mimicked_nodes, 'exo')
		
		dump_json(all_exo_dict_of_num_pos_pos_contacts_random_trials, all_exo_outfi_random_trials)
	else:
		all_exo_dict_of_num_pos_pos_contacts_random_trials = load_json(all_exo_outfi_random_trials)

	## All endo ##
	
	all_endo_dict_of_num_pos_pos_contacts_random_trials = {}
	all_endo_outfi_random_trials = osp.join(outdir, 'all_endo_random_trials_pos_pos_counts.json')
	
	if not osp.exists(all_endo_outfi_random_trials):
		print('\nRestype: All endo')
		
		# Shuffle nodes 
		endo_specific_nodes_anno_shuffled_repeatedly = shuffle_nodes_repeatedly(new_endo_specific_nodes)
		mimicked_nodes_anno_shuffled_repeatedly = shuffle_nodes_repeatedly(new_mim_endo_nodes)

		for t in mimicked_nodes_anno_shuffled_repeatedly:
			#1. go through the trials and get a new relabeled set of all_exo/all_endo edges for each trial
			all_endo_new_edges_anno_for_trial = relabel_edges_in_graph_based_on_shuffled_annos(endo_specific_nodes_anno_shuffled_repeatedly[t], mimicked_nodes_anno_shuffled_repeatedly[t], all_endo_edges)
			
			#2. Get number of pos-pos contacts in each label type for each randomization trial
			_, all_endo_dict_of_num_pos_pos_contacts_random_trials[t] = label_pos_pos_contacts(all_endo_new_edges_anno_for_trial, endo_specific_nodes, mimicked_nodes, 'endo')

		dump_json(all_endo_dict_of_num_pos_pos_contacts_random_trials, all_endo_outfi_random_trials)
	else:
		all_endo_dict_of_num_pos_pos_contacts_random_trials = load_json(all_endo_outfi_random_trials)


	# Get O/E for contact counts
	all_exo_dict_of_o_over_e = calc_o_over_e(all_exo_dict_of_num_pos_pos_contacts, all_exo_dict_of_num_pos_pos_contacts_random_trials, 'exo')
	all_endo_dict_of_o_over_e = calc_o_over_e(all_endo_dict_of_num_pos_pos_contacts, all_endo_dict_of_num_pos_pos_contacts_random_trials, 'endo')
	
	all_exo_out_oe = osp.join(outdir, 'all_exo_pospos_btwn_and_within_interface_oe_counts.json')
	all_endo_out_oe = osp.join(outdir, 'all_endo_pospos_btwn_and_within_interface_oe_counts.json')

	if not osp.exists(all_exo_out_oe) or not osp.exists(all_endo_out_oe):
		dump_json(all_exo_dict_of_o_over_e, all_exo_out_oe)
		dump_json(all_endo_dict_of_o_over_e, all_endo_out_oe)

	# Calc p-vals
	all_exo_comparison_types = [('exo', 'mimicked'), ('mimicked', 'btwn_exo_and_mimicked'), ('exo', 'btwn_exo_and_mimicked')]
	calc_pvals(all_exo_comparison_types, all_exo_dict_of_o_over_e)  # all significant

	all_endo_comparison_types = [('endo', 'mimicked'), ('mimicked', 'btwn_endo_and_mimicked'), ('endo', 'btwn_endo_and_mimicked')]
	calc_pvals(all_endo_comparison_types, all_endo_dict_of_o_over_e)  # all significant


if __name__ == '__main__':
	main()