#!/bin/bash

# Step 1: Compute all pairwise distances between residues based on the highest resolution PDB structure that these residues map to
python3 get_distances_between_residues.py

# Step 2: Construct interface-specific subgraphs, label nodes (residues) and edges (residue-residue contacts)
python3 construct_sub_int_cg.py 

# Step 3: Construct whole interface contact graph, label nodes (residues) and edges (residue-residue contacts)
python3 construct_whole_interface_cg.py 

# Step 4: Compute O/E contact count ratios for different interfaces and sub-interfaces, based on 10,000 intra-protein node-label randomization trials 
python3 oe_diff_interfaces.py 

# Step 5: Compute O/E pos-pos contact count ratios between distinct interfaces, based on 10,000 intra-protein node-label randomization trials (where individual interfaces are shuffled separately) 
python3 oe_pp_cross_int.py 
