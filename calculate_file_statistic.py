from biopandas.pdb import PandasPdb
from pathlib import Path
import argparse
from tqdm import tqdm
from pprint import pprint
from numpy import nan
from utils import init_global_logger, add_local_logger
# from pymol_util import calculate

import os
import json
import pandas as pd
import requests
import yaml

#setting up the logger


# read arguments from the command line
parser = argparse.ArgumentParser()
parser.add_argument('--data_path', type=str, default='data/pdbs/', help="Path to the folder with the preprocessed pdb files with hetatoms and protein separated.")
parser.add_argument('--lipid_paths', type=str, default='config/lipids_from_database.txt', help="Path to file where first three letters of each row are lipid codes.")
parser.add_argument('--result_path', type=str, default='results/all_data.csv', help="Path to resulting JSON file.")
parser.add_argument('--resn_inchi', type=str, default='config/lipids_from_database.txt')
parser.add_argument('--resn_to_smiles', type=str, default='config/resn_to_smiles.json')

def analyze_pdb(bdf: PandasPdb):

    # TODO: FIND MULTIPLE FRAMES IN ONE LIPID, TER DOES NOT WORK
    cols = ['residue_name' , 'insertion', 'chain_id', 'residue_number']
    filtered_df = bdf.df['HETATM'].groupby(cols).filter(lambda x: len(x) > 1)
    
    row_to_num_atoms = dict(filtered_df[cols].value_counts())

    # ATTEMPT TO REMOVE FRAMES USING TER
    s = set(filter(lambda x: x.startswith('MODEL'), bdf.pdb_text.split('\n')))
    s = set(map(lambda x: x.strip(), s))
    is_sequence = len(s) > 1

    filtered_df = filtered_df.groupby(cols).first().reset_index()[cols]
    filtered_df = filtered_df.loc[filtered_df['residue_name'] != 'HOH']

    # can be better
    filtered_df['num_atoms'] = [row_to_num_atoms[tuple(x)] for x in filtered_df[cols].values]
    filtered_df['is_sequence'] = is_sequence
    return filtered_df


def calculate_statistic(pdbs: list):
    all_data = []
    for i, pdb in tqdm(enumerate(pdbs)):
    # for i, pdb in enumerate(pdbs):
        if Path(pdb).exists:
            bdf = PandasPdb().read_pdb(str(pdb))
            pdb_code = pdb.stem
            pdb_path = pdb
        else:
            bdf = PandasPdb().fetch_pdb(str(pdb))
            pdb_code = pdb
            pdb_path = None
        unique_ligands = analyze_pdb(bdf)
        unique_ligands['pdb'] = pdb_code
        unique_ligands['pdb_path'] = pdb_path
        all_data.append(unique_ligands)

    return pd.concat(all_data)


def main():
    args = parser.parse_args()

    logger = init_global_logger()
    split_char = '/' if os.name != 'nt' else '\\'
    local_logger_filename = f'logs/{__file__.__repr__().split(split_char)[-1].split(".")[0]}.log'
    logger = add_local_logger(logger, local_logger_filename)

    logger.info(f'{__file__.__repr__().split("/")[-1][:-1]} started, with arguments: {args.__repr__()[10:-1]}')
    logger.info(f'Detailed logging to {local_logger_filename}')

    DATA_PATH = Path(args.data_path)
    RESULT_FILE = Path(args.result_path)
    RESN_TO_SMILES = Path(args.resn_to_smiles)

    resn_to_inchikey = dict()
    with open(Path(args.resn_inchi), 'r') as f:
        for line in f.readlines():
            resn, inchikey = line.split('\t')
            resn_to_inchikey[resn] = json.loads(inchikey)['input']

    # if RESN_TO_SMILES.exists():
    with open(RESN_TO_SMILES, 'r') as f:
        resn_to_smiles = json.load(f)

    data = calculate_statistic([x for x in DATA_PATH.iterdir() if x.suffix == '.pdb'])
    data['inchi_key'] = data['residue_name'].apply(lambda x: resn_to_inchikey[x] if x in resn_to_inchikey else None)
    data['smiles'] = data['residue_name'].apply(lambda x: resn_to_smiles[x] if x in resn_to_smiles else None)
    data['is_lipid'] = data['smiles'].apply(lambda x: x is not None)
    data.to_csv(RESULT_FILE, index=False)


if __name__ == '__main__':
    main()