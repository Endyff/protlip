from pymol import cmd
from pathlib import Path
import pandas as pd
import numpy as np
from tqdm import tqdm


def get_buried_ratio(pdb: [str, Path], selection, max_dist: float = 6, dot_solvent: int = 4, dot_density: int = 2):
    """
    Evaluate the binding between the lipid and the protein.
    """

    complex = 'complex'
    ligand = 'ligand'
    pocket = 'pocket'

    if Path(pdb).exists:
        cmd.load(pdb)
    else:
        cmd.fetch(pdb)
    # insertion = f'{selection.insertion}' if selection.insertion else '' # TODO
    resn = f'resn {selection.residue_name}' if selection.residue_name else ''
    chain = f'chain {selection.chain_id}' if selection.chain_id else ''
    resi = f'resi {selection.residue_number}' if selection.residue_number else ''
    
    # need this to work for some reason
    cmd.flag(flag='ignore', selection='all', action='clear')
    cmd.flag(flag='ignore', selection='solvent', action='set') 

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


    cmd.create(ligand, select)
    cmd.create(pocket, f'{ligand} around {max_dist} and polymer')
    cmd.create(complex, f'{ligand} or {pocket}')
    cmd.set('dot_solvent', dot_solvent)
    cmd.set('dot_density', dot_density)


    cmd.h_add()

    ligand_area = cmd.get_area(selection=ligand)
    pocket_area = cmd.get_area(selection=pocket)
    complex_area = cmd.get_area(selection=complex)

    cmd.reinitialize()
    if complex_area == 0:
        print('WARNING: Complex area is 0. Check the selection.')
        return np.nan

    buried_ratio = (ligand_area + pocket_area - complex_area) / (2*ligand_area) if pocket_area else 0

    return buried_ratio


# def get_buried_ratio(pdb: [str, Path], selection: [dict, pd.Row], max_dist: float = 6, dot_solvent: int = 4, dot_density: int = 2):
# 
df = pd.read_csv('results/simple_val_prediction.csv')
res = []
for i, row in tqdm(df.iterrows(), total=len(df)):
# for i, row in df.iterrows():
    # try:
    # print(i, end=' ')
    if i in [92, 93]: 
        res.append(np.nan)
        continue
    r = get_buried_ratio(row.pdb_path, row, max_dist=6, dot_solvent=4, dot_density=2)
    # print(r)

    # except Exception as e:
    #     print(e)
    #     with open('errors.txt', 'a') as f:
    #         f.write(f"complex: {row.complex_name}, error: {e}\n")
    #     r = np.nan
    # finally:
    res.append(r)    
df['buried_ratio'] = res
df.to_csv('results/simple_val_prediction_br.csv', index=False)
