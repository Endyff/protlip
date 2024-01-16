"""
This scripts downloads the pdb files from the list of pdb codes and ligand codes in 'data/pdb_lipid_complexes.txt' and saves them in 'data/pdbs/' folder.
"""

from pathlib import Path
import argparse
from Bio.PDB import PDBList
import tqdm
from utils import read_molecule_codes, OutputLogger, init_global_logger, add_local_logger
import logging
import os

# read argument from command line as a string coresponding to the path of the file with default value
# parser = argparse.ArgumentParser()
# parser.add_argument('--pdb_list_path', type=str, default='config/pdbs.txt', help="Path to the file with the list of pdb codes and ligand codes")
# parser.add_argument('--output_dir', type=str, default='data/pdbs_puress/', help="Path to the directory where to store the downloaded pdb files.")
# args = parser.parse_args()


# #setting up loggers
# logger = init_global_logger()
split_char = '/' if os.name != 'nt' else '\\'
# local_logger_filename = f'logs/{__file__.__repr__().split(split_char)[-1].split(".")[0]}.log'
# logger = add_local_logger(logger, local_logger_filename)

# logger.info(f'{__file__.__repr__().split("/")[-1][:-1]} started, with arguments: {args.__repr__()[10:-1]}')
# logger.info(f'Detailed logging to {local_logger_filename}')


# LIPID_COMPLEX_PDBS = Path(args.pdb_list_path)
# OUTPUT_DIR = Path(args.output_dir)

def download_pdb_files(pdb_list_path, output_dir, logger):
	with open(pdb_list_path, 'r') as f:
		molecule_codes = set(line.strip() for line in f.readlines())

	pdbl = PDBList()
	download_successful = 0
	download_errors = 0
	with OutputLogger(logger=logger, level='DEBUG') as log:

		for pdb_code in tqdm.tqdm(molecule_codes):
			if Path(f'{output_dir}/{pdb_code}.pdb').exists():
				continue
			for attmpt in range(3):
				path_name = pdbl.retrieve_pdb_file(pdb_code, pdir=output_dir, file_format='pdb', overwrite=False)
				path_name = Path(path_name)
				if path_name.exists():
					path_name.rename(Path(f'{output_dir}/{pdb_code}.pdb'))
					download_successful += 1
					logger.debug(f'File {pdb_code} downloaded successfully.')
					break
				else:
					download_errors += 1
					logger.debug(f'File {pdb_code} download failed. Attempt {attmpt+1}/3.')


	logger.info(f'{__file__.__repr__().split(split_char)[-1]} finished. Downloaded {download_successful} files. Encountered {download_errors} failed downloads.')



def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('--pdb_list_path', type=str, default='config/pdbs.txt', help="Path to the file with the list of pdb codes and ligand codes")
	parser.add_argument('--output_dir', type=str, default='data/pdbs/', help="Path to the directory where to store the downloaded pdb files.")
	args = parser.parse_args()

	logger = init_global_logger()
	split_char = '/' if os.name != 'nt' else '\\'
	local_logger_filename = f'logs/{__file__.__repr__().split(split_char)[-1].split(".")[0]}.log'
	logger = add_local_logger(logger, local_logger_filename)

	logger.info(f'{__file__.__repr__().split("/")[-1][:-1]} started, with arguments: {args.__repr__()[10:-1]}')
	logger.info(f'Detailed logging to {local_logger_filename}')

	download_pdb_files(args.pdb_list_path, args.output_dir, logger)
	
if __name__ == '__main__':
	main()