"""
Record the inter-object proximity of the atoms to compute the matrix of distances between the lipid and the protein atoms, with atom types marking the rows and the columns of the matrix -> for later evaluation of the tresholds
"""

from utils import pairwise_dist
from pymol import cmd
from pathlib import Path
import argparse
import numpy as np
import yaml
from tqdm import tqdm
import math
import os 

import pandas as pd

from utils import unite_selections, init_global_logger, add_local_logger

# #read arguments from the command line
# parser = argparse.ArgumentParser()
# parser.add_argument('--preprocessed_data_path', type=str, default='data/pdbs_preprocessed/', help="Path to the folder with the preprocessed pdb files with hetatoms and protein separated.")
# parser.add_argument('--output_path', type=str, default='results/clash_counts.txt', help="Path where to store the file with the number of potential clashes between the lipid and the protein atoms.")
# parser.add_argument('--atom_radius_type', type=str, default='cov', help="Type of the atom radius to use for the distance calculations. Possible values are 'cov' for covalent radius and 'vdw' for van der Waals radius.")
# parser.add_argument('--max_dist', type=int, default=6, help="Maximum distance between the lipid and the protein atoms to consider them for potential clash evaluation.")
# parser.add_argument('--slack_distance', type=float, default=0, help="Slack to be given to the treshold, to have manual control over it.")
# parser.add_argument('--store_distance_matrices', action='store_true', help="If the flag is provided, the distance matrices will be stored for each protein-lipid complex.")
# parser.add_argument('--atom_radii', default='conifg/atom_radii.yml', help="Path to .yml file containing atom radii.")
# args = parser.parse_args()

# #define constants, considering the arguments
# PREPROCESSED_DATA_PATH = Path(args.preprocessed_data_path)
# OUTPUT_DATA_PATH = Path(args.output_path)
# MAX_DIST = args.max_dist
# STORE_DISTANCE_MATRICES = args.store_distance_matrices
# STORE_DISTANCE_MATRICES = False
# RADIUS_TYPE = args.atom_radius_type.lower()
# TRESHOLD_SLACK = args.slack_distance
# ATOM_RADII = args.atom_radii

# OUTPUT_DATA_PATH = OUTPUT_DATA_PATH.with_stem(OUTPUT_DATA_PATH.stem+f'_{args.atom_radius_type}_slack_{TRESHOLD_SLACK}_max_dist_{MAX_DIST}')
# LIPID_SUFFIX = '_lipid'
# PROTEIN_SUFFIX = '_protein'
# POCKET_SUFFIX = '_pocket'

# LIPID_KEYWORD = 'lip'
# POCKET_KEYWORD = 'pocket'

# #setting up the logger]
# logger = init_global_logger()
# split_char = '/' if os.name != 'nt' else '\\'

# local_logger_filename = f'logs/{__file__.__repr__().split(split_char)[-1].split(".")[0]}.log'
# logger = add_local_logger(logger, local_logger_filename)

# logger.info(f'{__file__.__repr__().split("/")[-1][:-1]} started, with arguments: {args.__repr__()[10:-1]}')
# logger.info(f'Detailed logging to {local_logger_filename}')

with open('config/atom_radii.yml') as f:
    atom_radii = yaml.load(f, Loader=yaml.FullLoader)['cov']

def _pairwise_dist(sel1, sel2, max_dist, sidechain="N"):
	"""
	usage: pairwise_dist sel1, sel2, max_dist, [output=S/P/N, [sidechain=N/Y, [show=Y/N]]]
	sel1 and sel2 can be any to pre-existing or newly defined selections
	max_dist: maximum distance in Angstrom between atoms in the two selections
	--optional settings:
	sidechain: limits (Y) results to sidechain atoms (default N)
	show: shows (Y) individual distances in pymol menu (default=N)
	"""
	cmd.delete("dist*")
	extra=""
	if sidechain=="Y":
		extra=" and not name c+o+n"
	
	#builds models
	m1 = cmd.get_model(sel2+" around "+str(max_dist)+" and "+sel1+extra)
	m1o = cmd.get_object_list(sel1)
	m2 = cmd.get_model(sel1+" around "+str(max_dist)+" and "+sel2+extra)
	m2o = cmd.get_object_list(sel2)
 
	#controlers-1
	if len(m1o)==0: 
		print("warning, '"+sel1+extra+"' does not contain any atoms.")
		return None, None, None
	if len(m2o)==0: 
		print("warning, '"+sel2+extra+"' does not contain any atoms.")
		return None, None, None
	
	#measures distances
	distances= np.zeros((len(m1.atom),len(m2.atom)))
	row_atom_ids = [(m1o[0],m1.atom[r].chain,m1.atom[r].resn,m1.atom[r].resi,m1.atom[r].name,m1.atom[r].symbol) for r in range(len(m1.atom))]
	col_atom_ids = [(m2o[0],m2.atom[c].chain,m2.atom[c].resn,m2.atom[c].resi,m2.atom[c].name,m2.atom[c].symbol) for c in range(len(m2.atom))]

		
	for c1 in range(len(m1.atom)):
		for c2 in range(len(m2.atom)):
			distance=math.sqrt(sum(map(lambda f: (f[0]-f[1])**2, zip(m1.atom[c1].coord,m2.atom[c2].coord))))
			distances[c1,c2]=distance
	cmd.deselect()
	return distances, row_atom_ids, col_atom_ids



def eval_binding(pdb, selection, atom_radii, slack_distance=0, max_dist=10, lipid_keyword='lipid', pocket_keyword='pocket'):
    """
    Evaluate the binding between the lipid and the protein.
    """
    if Path(pdb).exists:
        cmd.load(pdb)
    else:
        cmd.fetch(pdb)
    # insertion = f'{selection.insertion}' if selection.insertion else '' # TODO
    resn = f'resn {selection.residue_name}' if selection.residue_name else ''
    chain = f'chain {selection.chain_id}' if selection.chain_id else ''
    resi = f'resi {selection.residue_number}' if selection.residue_number else ''
    # selection = 'and'.join([resn, chain, resi])
    select = ''
    if resn:
        select += resn
    if len(select) > 1:
        select += 'and '
    if chain:
        select += chain
    if len(select) > 1:
        select += 'and '
    if resi:
        select += resi

    # print(select)
    cmd.create(lipid_keyword, select)
    
    # cmd.create(pocket_keyword, f'{lipid_keyword} around {max(atom_radii.values()) + float(slack_distance)} and polymer')
    cmd.create(pocket_keyword, f'{lipid_keyword} around {max_dist} and polymer')
    cmd.h_add()

    distances, row_atom_ids, col_atom_ids = _pairwise_dist(lipid_keyword, pocket_keyword, max_dist)
    cmd.reinitialize()


    tresholds = np.zeros(distances.shape)
    for i, row in enumerate(row_atom_ids):
        for j, col in enumerate(col_atom_ids):
            row_atom_type = row[5].capitalize() #''.join(c for c in row[4] if not c.isnumeric()).capitalize()
            col_atom_type = col[5].capitalize() #''.join(c for c in col[4] if not c.isnumeric()).capitalize()

            row_radius = atom_radii[row_atom_type] if row_atom_type in atom_radii else None
            col_radius = atom_radii[col_atom_type] if col_atom_type in atom_radii else None
            
            if row_radius and col_radius:
                tresholds[i,j] = row_radius + atom_radii[col_atom_type] + slack_distance
            else:
                tresholds[i,j] = np.nan

            if not row_radius:
                print(f'Pdb {pdb}: Unknown atom type {row_atom_type}.')
            if not col_radius:
                print(f'Pdb {pdb}: Unknown atom type {col_atom_type}.')
    
    num_potential_clashes = np.sum(distances < tresholds)
    return num_potential_clashes

df = pd.read_csv('results/simple_val_prediction_br.csv')

res = []

# def eval_binding(pdb, selection, atom_radii, slack_distance=0, max_dist=10, lipid_keyword='lipid', pocket_keyword='pocket'):
# 
# for i, row in tqdm(df.iterrows(), total=len(df)):
for i, row in df.iterrows():
    # try:
    if i in [92, 93]: 
        res.append(np.nan)
        continue
    try:
        print(i, end=' ')

        r = eval_binding(row.pdb_path, row, atom_radii=atom_radii)
        print(r)

    except Exception as e:
        print(e)
        with open('errors.txt', 'a') as f:
            f.write(f"complex: {row.complex_name}, error: {e}\n")
        r = np.nan
    finally:
        res.append(r)    

df['clashes'] = res
df.to_csv('results/simple_val_prediction_clashes.csv', index=False)