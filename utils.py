from Bio.PDB.PDBIO import Select
import numpy as np
import logging
import contextlib

def read_molecule_codes(lipid_complex_pdbs_path):
	molecule_codes = []
	with open(lipid_complex_pdbs_path) as f:
		for line in f:
			_line = line.strip().split()
			pdb_code, ligand_code = _line[0], _line[1]
			molecule_codes.append((pdb_code,ligand_code))
	return molecule_codes

class Lipid_select(Select):
	def __init__(self, lipid_code):
		super(Lipid_select, self).__init__()
		self.lipid_code = lipid_code

	def accept_residue(self, residue):
		if residue.get_resname() == self.lipid_code:
			return 1
		else:
			return 0
		
from pymol import cmd, math

#TODO: cmd.select("b",  "prot and resn ILE and resi 21 and name CG2")

def unite_selections(to_select: list):
	"""
	to_select: list - this should contain a list of tuples of the form (obj, chain, resn, resi, name)
	"""
	
	def atom_selection_command(obj, chain, resn, resi, name):
		sele = f"{obj} and chain {chain} and resn {resn} and resi {resi} and name {name}"
		return sele
	
	commands_to_unite = []
	for obj, chain, resn, resi, name in to_select:
		sele = atom_selection_command(obj, chain, resn, resi, name)
		commands_to_unite.append(sele)

	union_command = " or ".join(commands_to_unite)
	return union_command

# below function is derived from https://pymolwiki.org/index.php/Pairwise_distances
def pairwise_dist(sel1, sel2, max_dist, sidechain="N"):
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

class OutputLogger:
	def __init__(self, logger, level="INFO"):
		self.logger = logger
		self.name = logger.name
		self.level = getattr(logging, level)
		self.handlers = self.logger.handlers
		self._redirector = contextlib.redirect_stdout(self)
		
	def write(self, msg):
		if msg and not msg.isspace():
			self.logger.log(self.level, msg)

	def flush(self): pass

	def __enter__(self):
		self._redirector.__enter__()
		return self

	def __exit__(self, exc_type, exc_value, traceback):
		# let contextlib do any exception handling here
		self._redirector.__exit__(exc_type, exc_value, traceback)

def init_global_logger(filepath = 'logs/global.log', name = 'root', logger_level = logging.DEBUG, handler_level = logging.INFO, format = '%(asctime)s %(levelname)s %(message)s'):
	logger = logging.getLogger(name)
	logger.setLevel(logger_level)
	glob_handler = logging.FileHandler(filepath)
	glob_handler.setLevel(handler_level)

	logform_glob = logging.Formatter(format)
	glob_handler.setFormatter(logform_glob)
	logger.addHandler(glob_handler)

	return logger

def add_local_logger(logger, filename, handler_level = logging.DEBUG, format = '%(asctime)s %(message)s'):
	loc_handler = logging.FileHandler(filename)
	loc_handler.setLevel(handler_level)
	logform_loc = logging.Formatter(format)
	loc_handler.setFormatter(logform_loc)

	logger.addHandler(loc_handler)
	return logger