import argparse
import pandas as pd
import numpy as np
from spyrmsd.rmsd import hrmsd
from biopandas.pdb import PandasPdb

from tqdm import tqdm


from pathlib import Path

def my_rmsd(pdb, selection, sdf, mean=False, remove_water=True):
    # Read PDB file using PDBParser
    parser = PandasPdb()
    if Path(pdb).exists():
        pdb_df = PandasPdb().read_pdb(pdb)
    else:
        pdb_df = PandasPdb().fetch_pdb(pdb)
    hetatm_df = pdb_df.df['HETATM']
    hetatm_df = hetatm_df[hetatm_df['residue_name'] == selection]
    if remove_water:
        hetatm_df = hetatm_df[hetatm_df['element_symbol'] != 'H']
    ground_truth_coords = hetatm_df[['x_coord', 'y_coord', 'z_coord']].values
    ground_truth_atoms = hetatm_df['element_symbol'] # TODO check

    # Read SDF file using Pandas
    
    with open(sdf) as f:
        lines = f.readlines()[3:]
    num_atoms = int(lines[0].split()[0].strip())

    atom_coords = [l.split()[:4] for l in lines[1:num_atoms+1]]
    ligand_df = pd.DataFrame(atom_coords, columns=['x_coord', 'y_coord', 'z_coord', 'atom'],)
    ligand_df = ligand_df.astype({'x_coord': float, 'y_coord': float, 'z_coord':float})
    if remove_water:
        ligand_df = ligand_df[ligand_df['atom'] != 'H']
    prediction_coords = ligand_df[['x_coord', 'y_coord', 'z_coord']].values
    prediction_atoms = ligand_df['atom']
    if mean:
        return np.linalg.norm(prediction_coords.mean(axis=0) - ground_truth_coords.mean(axis=0)) 

    return hrmsd(ground_truth_coords, prediction_coords, ground_truth_atoms, prediction_atoms)


def get_rmsd(df):
    prediction_base_path = Path('results/simple_val_sdf')
    # df = pd.read_csv('results/simple_val_prediction.csv')
    # sv_df = pd.read_csv('data/simple_val_diffdock_sdf.csv')
    # sv_df['pdb'] = sv_df['protein_path'].map(lambda x: Path(x).stem)
    # sv_df
    # df = df.merge(sv_df, on='pdb', how='right')
    rmsds = []
    for i, row in tqdm(df.head(1).iterrows(), total=len(df)):
        # if i>10:
        #     break
        # try:
        filename = Path(f'rank1_confidence{row.confidence_rank_1}.sdf')
        # print(row.keys())
        # print(row.complex_name)
        errors = []
        try:
            rmsd = my_rmsd(row.pdb_path, row.residue_name, prediction_base_path / Path(str(row.complex_name)) / filename, mean=False)
        except Exception as e:
            # errors.append(f"complex: {row.complex_name}, error: {e}")
            with open('errors.txt', 'a') as f:
                f.write(f"complex: {row.complex_name}, error: {e}\n")
            # print(f'row: {i}\n{e}')
            rmsd = np.nan
        finally:
            rmsds.append(rmsd)
    # print(errors)
        # rmsds.append(rmsd)
    # print(rmsds)
    df['rmsd'] = rmsds

def main():
    # parser = argparse.ArgumentParser(description='Load CSV file from path')
    # parser.add_argument('file_path', type=str, help='Path to the CSV file')
    # args = parser.parse_args()
    csv_file = Path('results/simple_val_prediction.csv')

    df = pd.read_csv(csv_file)
    # print(df.columns)
    get_rmsd(df)
    # print(df)
    # # for i, row in df.iterrows():
    # for i, row in tqdm(df.iterrows(), total=len(df)):
    #     # if i>10:
    #     #     break
    #     try:
            
    #         rmsd = my_rmsd(row.pdb_path, row.residue_name, prediction_base_path / Path(str(row.complex_name)) / Path('rank1.sdf'), mean=True)
    #     except Exception as e:
    #         print(f'row: {i}\n{e}')
    #         rmsd = np.nan
    #     finally:
    #         rmsds.append(rmsd)
    # df['rmsd'] = rmsds

    # # df['pdb'] = df['protein_path'].split('/')[-1].split('.')[0]
    


if __name__ == '__main__':
    main()
