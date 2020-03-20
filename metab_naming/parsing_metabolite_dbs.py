#!/usr/bin/env python
# Author: Ethan D. Evans, 2018 [eevans@mit.edu]

# FOR hmdb - file: hmdb_metabolites.xml
# FOR Metacyc - file: compounds.dat
# FOR ChEBI - file: ChEBI_complete.sdf
# for lipidMAP - file: structures.sdf
# This is just to clean-up the datafiles and save them for later compoud lookup

import argparse
import pickle

class Metabolite:
	def __init__(self):
		self.mass = None
		self.name = None
		self.formula = None
		self.inchi = None
		self.inchikey = None
		self.smiles = None
		self.charge = None
		self.category = None
		self.sub_class = None
		self.major_class = None
		self.super_class = None
		self.hmdb = None
		self.direct_parent = None
		self.synonyms = []
		self.specimen_location = []
		self.accession = []

	def to_dict(self):
		return {'mass':self.mass, 'location':self.specimen_location, 'formula':self.formula, 'inchikey':self.inchikey, 'hmdb_accession':self.accession,
		        'charge':self.charge, 'category':self.category, 'inchi':self.inchi, 'smiles':self.smiles, 'hmdb':self.hmdb, 'direct_parent':self.direct_parent,
		        'sub_class':self.sub_class, 'class':self.major_class, 'super_class':self.super_class, 'synonyms':self.synonyms}

def parse_hmdb(db_file):
	compound_dict = {}
	category = False
	charge = False
	with open(db_file) as f:
		location = False
		for line in f:
			line = line.strip()
			if line == '<metabolite>':
				metabolite = Metabolite()
			if '<name>' in line:
				if metabolite.name == None:
					metabolite.name = line[6:-7]
			if '<monisotopic_molecular_weight>' in line:
				metabolite.mass = float(line[30:-31])
			if '<chemical_formula>' in line:
				metabolite.formula = line[18:-19]
			if '<smiles>' in line:
				metabolite.smiles = line[8:-9]
			if '<inchi>' in line:
				metabolite.inchi = line[13:-8]
			if '<inchikey>' in line:
				metabolite.inchikey = line[10:-11]
			if '<kind>formal_charge</kind>' in line:
				charge = True
			if charge:
				if '<value>' in line:
					metabolite.charge = int(line[7:-8])
					charge = False
			if '<synonym>' in line:
				metabolite.synonyms.append(line[9:-10])
			if '<accession>' in line:
				metabolite.accession.append(line[11:-12])
			if '<direct_parent>' in line:
				metabolite.direct_parent = line[15:-16]
			if '<sub_class>' in line:
				metabolite.sub_class = line[11:-12]
			if '<class>' in line:
				metabolite.major_class = line[7:-8]
			if '<super_class>' in line:
				metabolite.super_class = line[13:-14]
			if '<biospecimen_locations>' in line:
				location = True
			if '</biospecimen_locations>' in line:
				location = False
			if location:
				if '<biospecimen>' in line:
					metabolite.specimen_location.append(line[13:-14])

			if line == '</metabolite>':
				if metabolite.mass == None or float(metabolite.mass) == 0.0:
					# print(metabolite.mass)
					pass
				else:
					compound_dict[metabolite.name] = metabolite.to_dict()
	p_file = open('./hmdb_metabolites.pkl','wb')
	pickle.dump(compound_dict, p_file)
	p_file.close()

def parse_lipidmap(db_file):
	compound_dict = {}
	with open(db_file) as f:
		name = False
		formula = False
		mass = False
		category= False
		major_class = False
		sub_class = False
		count = 0
		inchi_key = False
		inchi = False
		smiles = False
		for line in f:
			line = line.strip()
			if '> <LM_ID>' in line:
				count += 1
				metabolite = Metabolite()

			if name:
				metabolite.name = line
				name = False
			if line == '> <NAME>':
				name = True

			if mass:
				metabolite.mass = float(line)
				mass = False
			if line == '> <EXACT_MASS>':
				mass = True

			if formula:
				metabolite.formula = line
				formula = False
			if line == '> <FORMULA>':
				formula = True

			if category:
				metabolite.category = line
				category = False
			if line == '> <CATEGORY>':
				category = True

			if major_class:
				metabolite.major_class = line
				major_class = False
			if line == '> <MAIN_CLASS>':
				major_class = True

			if sub_class:
				metabolite.sub_class = line
				sub_class = False
			if line == '> <SUB_CLASS>':
				sub_class = True

			if inchi_key:
				metabolite.inchikey = line
				inchi_key = False
			if line == '> <INCHI_KEY>':
				inchi_key = True

			if inchi:
				metabolite.inchi = line[6:]
				inchi = False
			if line == '> <INCHI>':
				inchi = True

			if smiles:
				metabolite.smiles = line
				smiles = False
			if line == '> <SMILES>':
				smiles = True

			if line == '$$$$':
				if metabolite.mass == None or float(metabolite.mass) == 0.0:
					# print(metabolite.to_dict())
					pass
				else:
					compound_dict[metabolite.name] = metabolite.to_dict()
	p_file = open('./lipidmap_metabolites.pkl','wb')
	pickle.dump(compound_dict, p_file)
	p_file.close()

def parse_chebi(db_file):
	compound_dict = {}
	with open(db_file) as f:
		name = False
		formula = False
		mass = False
		charge = False
		inchi_key = False
		inchi = False
		smiles = False
		hmdb = False
		for line in f:
			line = line.strip()
			if '> <ChEBI ID>' in line:
				metabolite = Metabolite()

			if name:
				metabolite.name = line
				name = False
			if line == '> <ChEBI Name>':
				name = True

			if mass:
				metabolite.mass = float(line)
				mass = False
			if line == '> <Monoisotopic Mass>':
				mass = True

			if formula:
				metabolite.formula = line
				formula = False
			if line == '> <Formulae>':
				formula = True

			if charge:
				metabolite.charge = int(line)
				charge = False
			if line == '> <Charge>':
				charge = True

			if inchi_key:
				metabolite.inchikey = line
				inchi_key = False
			if line == '> <InChIKey>':
				inchi_key = True

			if inchi:
				metabolite.inchi = line[6:]
				inchi = False
			if line == '> <InChI>':
				inchi = True

			if smiles:
				metabolite.smiles = line
				smiles = False
			if line == '> <SMILES>':
				smiles = True

			if hmdb:
				metabolite.hmdb = line
				hmdb = False
			if line =='> <HMDB Database Links>':
				hmdb = True

			if line == '$$$$':
				if metabolite.mass == None or float(metabolite.mass) == 0.0:
					# print(metabolite.mass)
					pass
				else:
					compound_dict[metabolite.name] = metabolite.to_dict()
	p_file = open('./chebi_metabolites.pkl','wb')
	pickle.dump(compound_dict, p_file)
	p_file.close()
	
def parse_metacyc(db_file):
	compound_dict = {}
	with open(db_file, encoding="utf8", errors='ignore') as f:
	# with open(db_file) as f:
		for line in f:
			line = line.strip()
			if line[0] == '#':
				continue
			if 'UNIQUE-ID' in line:
				metabolite = Metabolite()
				metabolite.charge = 0
				metabolite.formula = ''
				metabolite.category = []
			if 'MONOISOTOPIC-MW' in line:
				metabolite.mass = float(line.split()[2])
			if 'TYPES' in line:
				if line[8:] == 'Compounds':
					continue
				metabolite.category.append(line[8:])
			if 'COMMON-NAME' in line:
				metabolite.name = line[14:]
			if 'ATOM-CHARGES' in line:
				metabolite.charge += int(line.split()[3][:-1])
			if 'CHEMICAL-FORMULA' in line:
				metabolite.formula += line.split()[2][1:]
				metabolite.formula += line.split()[3][:-1]
			if 'SMILES' in line:
				metabolite.smiles = line[9:]
			if 'INCHI - InChI' in line:
				metabolite.inchi = line[14:]
			if 'INCHI-KEY' in line:
				metabolite.inchikey = line[21:]
			if line == '//':
				if metabolite.mass == None or float(metabolite.mass) == 0.0:
					# print(metabolite.mass)
					pass
				else:
					compound_dict[metabolite.name] = metabolite.to_dict()
					# print(compound_dict[metabolite.name])
	p_file = open('./metacyc_metabolites.pkl','wb')
	pickle.dump(compound_dict, p_file)
	p_file.close()

def parse_db(db, db_file):
	if db == 'metacyc':
		parse_metacyc(db_file)
	elif db == 'hmdb':
		parse_hmdb(db_file)
	elif db =='lipidmap':
		parse_lipidmap(db_file)
	else:
		parse_chebi(db_file)

def main():
	'''
	Define the files you want to parse and make into common look up format. 
	# FOR hmdb - file: hmdb_metabolites.xml
	# FOR Metacyc - file: compounds.dat
	# FOR ChEBI - file: ChEBI_complete.sdf
	# for lipidMAP - file: structures.sdf
	# This is just to clean-up the datafiles and save them for later compoud lookup
	'''
	parser = argparse.ArgumentParser()
	parser.add_argument('-d', help='path + database file you want to parse')
	parser.add_argument('-t', help='which database are you using? one of: "metacyc", "hmdb", "chebi", "lipidmap"')
	args = parser.parse_args()
	db_file = args.d
	db = args.t

	parse_db(db, db_file)		

if __name__ == '__main__':
	main()

