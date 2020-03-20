#!/usr/bin/env python
# Author: Ethan D. Evans, 2019 [eevans@mit.edu]

import pandas as pd
import numpy as np
import pickle
import os
from multiprocessing import Pool
import argparse 
import math
import ast
from rdkit import Chem 
import re
import string

# A few assumptions this program takes:
# 1) we are only considering +1 / -1 / 'neutral' (natively +/-1) ions. not trying to deal with multiple charge states for each feature
#		This may be added in sometime in the future
# 2) There are only certain adducts we use for comparison (pos: +H > +Na > +K ~ +NH4 and none, neg: -H > +Cl ~ none > -H-H2O,)
#		If there is one you want support for please ask for it
# 3) you need to define the RT window to be considered in a peak family
# 4) the voting can use input...but follows the idea that certain adducts are more plausible (see tanking for 2), internal to this its ppm error sorted 


def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('-t', help='ppm tolerance to limit metabolites, default=5', default=5)
	parser.add_argument('-r', help='RT tie to define what should be considers part of the same peak grouping, default=10.0', default=5.0)
	parser.add_argument('-f', help='specific file to run on or a directory of ONLY these files', default='./all/')
	parser.add_argument('-p', help='pickle file used for adding in mz/rt features if needed / performing other operations', default='none')
	parser.add_argument('-n', help='is this for the meta-analysis or some other analysis? like single or a few files', default='false')
	parser.add_argument('-s', help='seperator', default='\t')
	parser.add_argument('-o', help='path to where you want output file, best with single file use', default='./')
	args = parser.parse_args()
	meta_analysis = args.n
	seperator = args.s
	out_path = args.o
	if meta_analysis == 'false':
		meta_analysis = False
	else:
		meta_analysis = True
	ppm = float(args.t)
	rt_window = float(args.r)
	if args.p == 'none':
		models = 'none'
	else:
		models = pickle.load(open(args.p, 'rb'))
	initial = False
	if args.f[-3:] == 'csv':
		initial = True
		files = {args.f: pd.read_csv(args.f, sep=seperator)}
	elif args.f[-3:] == 'pkl':
		files = pickle.load(open(args.f, 'rb'))
	else:
		initial = True
		files = os.listdir(args.f)
		files = {fi:pd.read_csv(os.path.join(args.f,fi), sep=seperator) for fi in files}
	return ppm, rt_window, models, files, initial, meta_analysis, out_path

def get_db():
	'''
	Read in the assocaited database file
	'''
	meta_db = pickle.load(open('./metacyc_metabolites.pkl', 'rb'))
	chebi_db = pickle.load(open('./chebi_metabolites.pkl', 'rb'))
	hmdb_db = pickle.load(open('./hmdb_metabolites.pkl', 'rb'))
	lipidmap_db = pickle.load(open('./lipidmap_metabolites.pkl', 'rb'))
	return meta_db, chebi_db, hmdb_db, lipidmap_db

def add_info(column, db):
	'''
	This is the function that takes a column of data (say H_chebi) and looking in that db for all the compounds and adds the associated info
	'''
	new_column = []
	for entry in column:
		newline = []
		if entry == '[]':
			new_column.append('[]')
			continue
		for ele in ast.literal_eval(entry):
			smiles = 'None'
			formula = 'None'
			try:
				formula = db[ele[:-7]]['formula']
			except:
				try:
					formula = db[ele[:-8]]['formula']
				except:
					pass
			formula = make_common_formula(formula)
			try:
				smiles = db[ele[:-7]]['smiles']
			except:
				try:
					smiles = db[ele[:-8]]['smiles']
				except:
					pass
			try:
				charge = db[ele[:-7]]['charge']
			except:
				try:
					charge = db[ele[:-8]]['charge']
				except:
					pass
			if type(charge) != int:
				try:
					m = Chem.MolFromSmiles(smiles)
					charge = Chem.rdmolops.GetFormalCharge(m) 
				except:
					charge = 0 # best guess and normally probably good
					pass
			try:
				newline.append([ele[:-7],ele[-5:],formula,smiles, charge]) #,db[ele[:-7]]['inchi'],db[ele[:-7]]['inchikey']))
			except: 
				pass
		new_column.append(str(newline))
	return new_column

def add_inchi_smiles_formula(files):
	''' 
	basically take all the assignment features and get associated info like the inchi char, smile string and charge for later use
	CHARGE is super important for removing features which make no sense as a +1 or -1 or neutral since this is an assumption of the program
	Inputs:
		files - dictionary mapping the name of the dataset to the output of the mapping_mz_to_metablites.py program
	output;
		files - but not with added info at all the correct spots!
  	'''

	meta_db, chebi_db, hmdb_db, lipidmap_db = get_db()
	header_to_db = {'lipid':lipidmap_db, 'chebi':chebi_db, 'hmdb':hmdb_db, 'meta':meta_db}
	for fi in files:
		for header in list(files[fi]):
			split_header = header.split('_')
			if split_header[-1] not in ['lipid', 'chebi', 'hmdb', 'meta']:
				continue
			files[fi][header] = add_info(list(files[fi][header].values), header_to_db[split_header[-1]])
	return files

def extract_formula(metab):
	headers = ['H_meta', 'Na_meta', 'K_meta', 'NH4_meta', 'none_meta', 'C13_H_meta',
			   'H_chebi', 'Na_chebi', 'K_chebi', 'NH4_chebi', 'none_chebi','C13_H_chebi',
			   'H_lipid', 'Na_lipid', 'K_lipid', 'NH4_lipid', 'none_lipid','C13_H_lipid',
			   'H_hmdb', 'Na_hmdb', 'K_hmdb', 'NH4_hmdb', 'none_hmdb','C13_H_hmdb',
			   'M-H_meta', 'Cl_meta', 'M-H+H2O_meta','C13-H_meta',
			   'M-H_chebi', 'Cl_chebi', 'M-H+H2O_chebi','C13_H_chebi',
			   'M-H_lipid', 'Cl_lipid', 'M-H+H2O_lipid','C13_H_lipid',
			   'M-H_hmdb', 'Cl_hmdb', 'M-H+H2O_hmdb', 'C13_H_hmdb']
	formulas = []
	for col in metab:
		if col not in headers:
			continue
		try:
			data = metab[col].values[0]
			if data == '[]':
				continue
			data = ast.literal_eval(data)
			for d in data:
				if d[2] not in formulas:
					formulas.append(d[2])
		except:
			pass
	return formulas

def make_common_formula(orig_formula):
	elements = ['Ac','Mt','Ag','N','Al','Na','Am','Nb','Ar','Nd','As','Ne','At','Au','B','Ba','Be',
				'O','Og','Bi','P','Br','Pa','C','Pb','Ca','Pd','Cd','Pm','Ce','Pr','Cl','Pt',
				'Cm','Pu','Ra','Rb','Cr','Re','Rf','Cu','Rg','Rh','Db','Rn','Ds','Ru','Dy','S','Er','Sb',
				'Es','Eu','Se','F','Sg','Fe','Si','Fl','Sm','Fm','Fr','Sr','Ga','T','Gd','Ta','Ge','Tb','H','Tc',
				'He','Te','Th','Hg','Ti','Tl','Tm','I','Ts','U','Ir','V','K','W','Kr','Xe','La','Y','Li',
				'Yb','Lr','Zn','Lu','Zr','Lv','Mc','Md','Mg','Mn','Mo', 'R', 'X'] # removed speciic elements to avoid upper to lower case conversions
	if orig_formula[0] == '(':
		return orig_formula


	chem_f = orig_formula.split('.')
	if len(chem_f) > 1:
		combined = []
		for part in chem_f:
			if 'In' in part:
				print(part, chem_f)
			new_form = []
			sub_chem_f = re.split('(\d+)',part)
			if sub_chem_f[-1] == '':
				sub_chem_f = sub_chem_f[:-1]
			if sub_chem_f[0] == '':
				sub_chem_f = sub_chem_f[1:]
			try:
				r_ele = int(sub_chem_f[-1])
			except:
				sub_chem_f.append('1')
			try:
				l_ele = int(sub_chem_f[0])
				for ele_sub in sub_chem_f[1:]:
					try:
						if int(ele_sub):
							new_form.append(str(l_ele*int(ele_sub)))
					except:
						if ele_sub in elements:
							new_form.append(ele_sub) 
						elif len(ele_sub) == 2:
							new_ele =ele_sub[0]+ele_sub[1].lower()
							if new_ele in elements:
								new_form.append(new_ele)
						elif ele_sub == []:
							pass
						else:
							try:
								sub = re.findall('[A-Z][^A-Z]*',ele_sub)
								for e_sub in sub[:-1]:
									if e_sub in elements:
										new_form.append(e_sub)
										new_form.append(str(l_ele*1))
									else:
										pass
								new_form.append(sub[-1])
							except:
								pass
			except:
				combined.append(sub_chem_f)

		combined = [item for sublist in combined for item in sublist]
		formula = ''.join(combined)

	else:
		formula = chem_f[0]	

	chem_f = re.split('(\d+)',formula)
	if chem_f[-1] == '':
		chem_f = chem_f[:-1]
	try:
		l_ele = int(chem_f[-1])
	except:
		chem_f.append('1')
	new_form = []
	for ele in chem_f:
		try:
			if int(ele):
				new_form.append(ele)	
		except:	
			if ele in elements:
				new_form.append(ele)
			elif ele == 'NO':
				new_form.append(ele[0])
				new_form.append(str(1))
				new_form.append(ele[1])
			elif len(ele) == 2:
				new_ele = ele[0]+ele[1].lower()
				if new_ele in elements:
					new_form.append(new_ele)
			elif ele == []:
				pass
			else:
				try:
					sub = re.findall('[A-Z][^A-Z]*', ele)
					for ele_e in sub[:-1]:
						if ele_e in elements:
							new_form.append(ele_e)
							new_form.append('1')
						else:
							pass
					new_form.append(sub[-1])
				except:
					pass
	return ''.join(new_form)

def find_metab_support(sig_metab_formulas, support_metab_formulas):
	count = 0
	for f_1 in support_metab_formulas:
		if f_1 in sig_metab_formulas:
			count += 1 
	if count != 0:
		return True
	else:
		return False

def get_adduct_data(add_data, mask):
	combined_data = []
	for ele, mask_lab in zip(add_data,mask):
		if ele == '[]':
			continue
		ele = ast.literal_eval(ele)
		for sec_ele in ele:
			combined_data.append((sec_ele[1], sec_ele[0], sec_ele[2], mask_lab, sec_ele[4]))
	return sorted(combined_data, key=lambda x: x[0])

def determine_polarity(col_headers):
	if 'Na_meta' in col_headers:
		return 'positive'
	else:
		return 'negative'

def guess_metabolite(most_imp, sec_imp, third_imp, mode):
	def parse_imp_list(imp, mode):
		if mode == 'positive':
			headers = ['H_', 'none_', 'Na_', 'K_', 'NH4_', 'C13_H_','2C13_H_', '3C13_H_', 'M_ACN_H_', '2M_H_', '2M_Na_', '2M_ACN_H_']
		else:
			headers = ['M-H_', 'none_', 'Cl_', 'M-H-H2O_', 'C13-H_','2C13-H_', '3C13-H_', '2M-H_', 'M-2H_Na_', 'M-2H_K_']

		summary = {}
		compound_list = []
		for ele in imp:
			if 'R' in ele[2]:
				continue 
			if 'Pr' in ele[2]:
				continue
			if mode == 'positive' and int(ele[4]) < 0:
				continue
			if mode == 'negative' and int(ele[4]) > 0:
				continue
			if mode == 'positive' and int(ele[4]) > 1:
				continue
			if mode == 'negative' and int(ele[4]) < -1:
				continue
			if mode == 'positive' and int(ele[4]) == 1:
				if 'none_' not in ele[3]:
					continue
			if mode == 'negative' and int(ele[4]) == -1:
				if 'none_' not in ele[3]:
					continue
			if (int(ele[4]) == 0 and 'none_' in ele[3]):
				continue

			if ele[2] not in summary:
				compound_list.append((ele[2], ele[0]))
			if ele[2] in summary:
				summary[ele[2]].append((ele[0], ele[1], ele[3], ele[4]))
			else:
				summary[ele[2]] = [(ele[0], ele[1], ele[3], ele[4])]
		compound_list = sorted(compound_list, key=lambda x: x[1])
		compound_list = [ele[0] for ele in compound_list]
		return summary, compound_list
	summary_most, form_list_most = parse_imp_list(most_imp, mode)
	summary_sec, form_list_sec = parse_imp_list(sec_imp,mode)
	summary_third, form_list_third = parse_imp_list(third_imp,mode)
	combined_form = form_list_most + form_list_sec + form_list_third
	combined_summary = {**summary_most, **summary_sec, **summary_third}
	best_guess_list = [(ele,combined_summary[ele]) for ele in combined_form]
	return best_guess_list

def vote_positive(family, mode):
	# prioritze H, then none, then Na, then NH4/K, go by ppm first ...for neg its -H, none, Cl, other
	col_mask = []

	for col in family:
		if col in ['H_meta', 'H_chebi', 'H_hmdb','H_lipid']:
			col_mask.append(col)
	proton_add = family[col_mask]

	proton_add = get_adduct_data(proton_add.values[0],col_mask)
	col_mask = []
	secondary = ['none_' , 'Na_' , 'M_ACN_H_' , '2M_H_','K_','NH4_', '2M_Na_','2M_ACN_H_']
	for col in family:
		if any(name in col for name in secondary):
			col_mask.append(col)
	sec_tier = family[col_mask]
	sec_tier = get_adduct_data(sec_tier.values[0],col_mask)
	col_mask = []
	third = ['C13_H_' , '2C13_H_' ,'3C13_H_' ]
	for col in family:
		if any(name in col for name in third):
			col_mask.append(col)
	rest_add = family[col_mask]
	rest_add = get_adduct_data(rest_add.values[0],col_mask)
	return guess_metabolite(proton_add, sec_tier, rest_add, mode)

def vote_negative(family,mode):
	col_mask = []
	primary = ['M-H_meta', 'M-H_chebi', 'M-H_hmdb','M-H_lipid']
	for col in family:
		if any(name in col for name in primary):
			col_mask.append(col)
	proton_add = family[col_mask]
	proton_add = get_adduct_data(proton_add.values[0],col_mask)
	col_mask = []
	secondary = ['Cl_', 'none_', 'M-H-H2O_']
	for col in family:
		if any(name in col for name in secondary):
			col_mask.append(col)
	sec_tier = family[col_mask]
	sec_tier = get_adduct_data(sec_tier.values[0],col_mask)
	col_mask = []
	third = ['M-2H_Na_','2M-H_', 'M-2H_K_', 'C13-H_', '2C13-H_', '3C13-H_']
	for col in family:
		if any(name in col for name in third):
			col_mask.append(col)
	rest_add = family[col_mask]
	rest_add = get_adduct_data(rest_add.values[0],col_mask)
	return guess_metabolite(proton_add, sec_tier, rest_add, mode)

def vote_without_evidence(family):
	mode = determine_polarity(list(family))
	if mode == 'positive':
		metabolite_guess = vote_positive(family, mode)
	else:
		metabolite_guess = vote_negative(family,mode)
	return metabolite_guess

def ds_voting_protocol(sig_metabs, data, mode, rt, ds, out_path='./meta_analysis/'):
	def drop_col(df):
		try:
			df.drop(columns=['level_0'], inplace=True)
		except:
			pass
		try:
			df.drop(columns=['mz.1'], inplace=True)
		except:
			pass
		try:
			df.drop(columns=['Unnamed: 0'], inplace=True)
		except:
			pass
		try:
			df.drop(columns=['index'], inplace=True)
		except:
			pass
		return df 

	rt_label = ['rt', 'retention_time', 'Retention Time', 'rtmed','RT', 'row retention time', 'retention index', 'ri', 'retention time']
	mask_neighbors = np.asarray([False for i in range(data[2].shape[0])])
	header_to_keep = [header for header in list(sig_metabs) if 'meta' not in header and 'hmdb' not in header and 'lipid' not in header and 'chebi' not in header]
	rt_col = False
	all_meta_data = []
	all_feature_data = []
	feature_type = []
	for col_header in list(sig_metabs): # can you jsut make this a list of sig feat and model feat?!?!?!?!
		if col_header in rt_label:
			# once you find the right column, get the rt times for both the sig metabs and the all metabs
			sig_metabs_rt_data = list(sig_metabs[col_header].values) # list of all SIGNIFICANT (or model when added) FEATURES to look at but their RT time
			all_rt_data = list(data[2][col_header].values)  # list of all the mz features - BUT THEIR RT time
			# go over all data of peak/rt, if within rt window and molecular formulas help id a target feature, keep it
			for j,rt_time in enumerate(sig_metabs_rt_data):
				sig_metab_neighbors = []
				sig_metab_formulas = extract_formula(sig_metabs.iloc[[j]])

				for i,ele in enumerate(all_rt_data):
					if rt_time - rt <= ele <= rt_time + rt:
						support_formulas = extract_formula(data[2].iloc[[i]])
						support = find_metab_support(sig_metab_formulas, support_formulas)
						# MAKE A NEW LIST OF SIG INDICIES WITH SUPPORT....
						# once you have iterated over the all data, if support, vote with it, otherwise just vote without support
						if support:
							sig_metab_neighbors.append(i)
							mask_neighbors[i] = True

				if len(sig_metab_neighbors) != 0:
					family = data[2].iloc[sig_metab_neighbors]
					if family.shape[0] > 1:
						meta_data = sig_metabs[header_to_keep].iloc[[j]]
						vote_feature_data = vote_without_evidence(family.loc[[j]])
						all_feature_data.append(vote_feature_data)	
						all_meta_data.append(meta_data)		
						feature_type.append('primary')
					else:
						meta_data = sig_metabs[header_to_keep].iloc[[j]]
						vote_feature_data = vote_without_evidence(family)
						all_feature_data.append(vote_feature_data)	
						all_meta_data.append(meta_data)		
						feature_type.append('primary')
				if len(sig_metab_neighbors) == 0:
					meta_data = sig_metabs[header_to_keep].iloc[[j]]
					all_meta_data.append(meta_data)
					all_feature_data.append([])	
					feature_type.append('none')
			rt_col = True
			break
	# if you happen to just want all the sig peak and their associates....		
	# now some datasets lack the rt info...so just do single, without evidence compound mz look up....
	if not rt_col:
		# this is for the datasets where there is no RT information!
		for i in list(sig_metabs.index):
			family = sig_metabs.iloc[[i]]		
			meta_data = sig_metabs[header_to_keep].iloc[[i]]
			metabolite_vote_list = vote_without_evidence(family)
			all_meta_data.append(meta_data)
			all_feature_data.append(metabolite_vote_list)
			feature_type.append('primary')
	if len(all_meta_data) > 0:
		meta_data_df = pd.concat(all_meta_data, axis=0, sort=False)
		meta_data_df = drop_col(meta_data_df)
		meta_data_df = meta_data_df.reset_index()
		feature_df = pd.DataFrame(all_feature_data)
		columns_headers = ['vote_{}'.format(i) for i in range(feature_df.shape[1])]
		feature_df.columns = columns_headers
		combined_all_df = pd.concat([meta_data_df, feature_df], axis=1)
		combined_all_df = drop_col(combined_all_df)
		combined_all_df.insert(loc=0, column='support_type', value=feature_type)
		combined_all_df.to_csv(out_path+"{}_{}.csv".format(mode, ds))

def non_metaanalysis_voting(files, rt_metadata, out_path):
	for study in files:
		name = study.split('/')
		if name[-1][-3:] == 'csv':
			name = 'voted_metabs_'+name[-1][:-4]
		for mode in ['p_values', 'model_feat']:
			if mode in list(files[study]):
				print(mode)
				if mode == 'p_values':
					spec_metabs = files[study][mode] < 0.05
					ds_voting_protocol(spec_metabs, files[study], mode, rt_metadata, study)
				if mode == 'model_feat':
					spec_metabs = files[study][mode] != 0
					ds_voting_protocol(spec_metabs, files[study], mode, rt_metadata, study)
		ds_voting_protocol(files[study], [np.asarray([True for i in range(files[study].shape[0])]),[],files[study]], 'all', rt_metadata, name, out_path)

def main():
	### note 'meta_analysis' is a remnant from an alternative version that works on very different data. I have simplified it, just dont set meta_analysis to True otherwise it wont work lol
	### similarly ppm is not used anymore
	# for a one off file or directory...
	# with an output path:
	# python ./metabolite_calling.py -t 5 -r 10 -f "./underworlds/longitudinal_MIT_50000noise_metabolites_ms2.csv" -s ',' -o './underworlds/'
	# or once run:
	# python ./metabolite_calling.py -t 5 -r 10 -f "./underworlds/underworlds_charged_smiles_formula_metab_lookup.pkl" -s ',' -o './underworlds/'
	
	ppm, rt_window, models, files, initial, meta_analysis, out_path = parse_args()
	# if this is the first time running, just run on the desired file or directory and this will make a 
	# pickle of the same data but with the added info (formula / smiles) and save it... 
	if initial:
		files = add_inchi_smiles_formula(files)
		if not meta_analysis:
			pickle.dump(files, open('./underworlds/underworlds_charged_smiles_formula_metab_lookup.pkl', 'wb'))

	study_rt_metadata = rt_window
	# need to make the metabs_with_masks some how...
	non_metaanalysis_voting(files, study_rt_metadata, out_path) 

if __name__ == '__main__':
	main()
