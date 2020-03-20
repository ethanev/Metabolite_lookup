#!/usr/bin/env python
# Author: Ethan D. Evans, 2019 [eevans@mit.edu]

import pandas as pd
import numpy as np
import argparse
import ast
import random
import os
import pickle

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('-f', help='csv files from the metabolite_calling.py or csv with list of mzs of iterest')
	parser.add_argument('-n', help='is this the first time you are this', default=False)
	parser.add_argument('-c', help='which column is the mass in the csv?', default=0)
	parser.add_argument('-r', help='which column is the rt in the csv ', default=0)
	parser.add_argument('-m', help='path to metfrag outputs', default='./msms_out/')
	parser.add_argument('-v', help='is this with voted data or just simple mz csv, this assumes it is, enter anything to alter this', default='true')
	parser.add_argument('-i', help='what ion mode, positive or negative?', default='negative')
	parser.add_argument('-o', help='outfile path and name', default='./you_should_fill_this_out.csv')
	args = parser.parse_args()
	voted = args.v
	mode = args.i
	first = args.n
	out_file = args.o
	if first == 'True':
		first = True
	mass_column = int(args.c)
	rt_column = int(args.r)
	mass_file = args.f
	metfrag_res_path = args.m
	return mass_file, mass_column, rt_column, metfrag_res_path, voted, mode, first, out_file

def add_inchi(df):
	''' 
	Inputs:
		df - with entries that need inchi or whatever
	output;
		df - with inchi or whatever added
  	'''
	meta_db, chebi_db, hmdb_db, lipidmap_db = get_db()
	header_to_db = {'lipid':lipidmap_db, 'chebi':chebi_db, 'hmdb':hmdb_db, 'meta':meta_db}
	for header in list(df)[4:]:
		for i,ele in enumerate(df[header]):
			if type(ele) == float:
				continue
			if ele == str(('', [])):
				continue
			else:
				ele = ast.literal_eval(ele)
				comps = ele[1]
				comp_list_with_inchi = []
				for comp in comps:
					comp_inchi = []
					for db in [meta_db, chebi_db, hmdb_db, lipidmap_db]:
						try:
							inchi = db[comp[1]]['inchi']
							if inchi[:5] == 'NCHI ':
								inchi = inchi[7:]
							if inchi[:3] == '1S/':
								inchi = 'InChI='+inchi
							if inchi in comp_inchi:
								pass
							else:
								comp_inchi.append(inchi)
						except:
							pass
					comp_list_with_inchi.append((comp[0], comp[1],comp[2],comp[3],comp_inchi))
				entry = (ele[0],comp_list_with_inchi)
				df.loc[i,header] = str(entry)
	return df

def get_db():
	'''
	Read in the assocaited database file
	'''
	meta_db = pickle.load(open('./metacyc_metabolites.pkl', 'rb'))
	chebi_db = pickle.load(open('./chebi_metabolites.pkl', 'rb'))
	hmdb_db = pickle.load(open('./hmdb_metabolites.pkl', 'rb'))
	lipidmap_db = pickle.load(open('./lipidmap_metabolites.pkl', 'rb'))
	return meta_db, chebi_db, hmdb_db, lipidmap_db

def parse_metfrag_output(metfrag_files, metfrag_res_path):
	'''
	input:
		metfrag_files: just a list of all the metfrag output files
		metfrag_res_path: the path to the metfrag output files
	Output:
		metfrag_results: a dictionary that maps an mz/rt tuple to a list of dictionaries for the different metfrag results 
	'''
	metfrag_results = {}
	# adducts = {'H':-1.007276, 'Na': -22.989218, 'K': -38.963158, 'M-H-H2O': 19.01839,
	# 		'NH4': -18.033823, 'M-H': 1.007276, 'Cl': -34.969402, 'none':0.0, 
	# 		'C13-H':1.007276, 'C13_H':-1.007276, '2C13_H':-1.007276, '2C13-H':1.007276, 
	# 		'3C13_H':-1.007276, '3C13-H':1.007276, 'M_ACN_H':-42.034276, '2M-H':1.007276, 
	# 		'2M_H': -1.007276, '2M_Na':-22.989218, '2M_ACN_H':-42.034276, 'M-2H_Na':-20.974666, 
	# 		'M-2H_K':-36.948606 } 
	adducts = {-1.007:'H', 1.007: 'M-H', -22.989:'Na', -38.963:'K', 19.018:'M-H-H2O',
			-18.034:'NH4', -34.969:'Cl',0.0:'none',-42.034:'M_ACN_H', -20.975:'M-2H_Na',
			-36.949:'M-2H_K'}
	for f in metfrag_files:
		if f[-3:] != 'csv':
			continue 
		f_components = f.split('_')
		repeat_num = f_components[0]
		vote_num = f_components[1]+'_'+f_components[2]
		orig_mass = f_components[8]
		neut_mass = f_components[5]
		mass_diff = float(neut_mass) - float(orig_mass)
		if abs(mass_diff) > 50:
			neut_mass = 2*float(neut_mass)
			mass_diff = float(neut_mass) - float(orig_mass)
		mass_add = None
		for add in adducts:
			if add-0.001 <= mass_diff <= add+0.001:
				mass_add = adducts[add]
		if mass_add == None:
			print(neut_mass, orig_mass)
		rt = f_components[10]
		db = f_components[-1][:-4]

		metfrag_res = pd.read_csv(os.path.join(metfrag_res_path,f))
		if metfrag_res.shape[0] == 0:
			continue
		names = ['Score', 'Name', 'MolecularFormula', 'Identifier', 'SMILES', 'InChIKey','InChI', 'CompoundName', 'NoExplPeaks','NumberPeaksUsed']
		names_to_use = []
		for name in metfrag_res:
			if name in names:
				names_to_use.append(name)
		metfrag_res = metfrag_res[names_to_use]
		to_replace = {'CompoundName':'Name'}
		metfrag_res = metfrag_res.rename(columns=to_replace)
		mz_rt_info = {'results':metfrag_res,'repeat':repeat_num,'vote':vote_num,'db':db, 'adduct':mass_add}
		if (orig_mass, rt) not in metfrag_results:
			metfrag_results[(orig_mass,rt)] = [mz_rt_info]
		else:
			metfrag_results[(orig_mass,rt)].append(mz_rt_info)
	return metfrag_results

def reduce_redundancies(list_of_redun_metabs):
	list_of_redun_metabs = [(ele1.lower(),ele2.lower()) for ele1,ele2 in list_of_redun_metabs]
	grouped = {}
	seen = []
	for same_metabs in list_of_redun_metabs:
		if same_metabs[0] not in grouped and same_metabs[0] not in seen and same_metabs[1] not in grouped and same_metabs[1] not in seen:
			grouped[same_metabs[0]] = [same_metabs[1]]
			seen.append(same_metabs[0])
			seen.append(same_metabs[1])
		elif same_metabs[0] in grouped:
			grouped[same_metabs[0]].append(same_metabs[1])
		elif same_metabs[1] in grouped:
			grouped[same_metabs[1]].append(same_metabs[0])
		else:
			pass
	return grouped.keys()

def min_metfrag_other(metfrag_guess):
	min_guess = {}
	for vote in metfrag_guess:
		all_guesses = []
		number = None
		i = 0
		for group in metfrag_guess[vote]:
			if i == 0:
				all_guesses.extend(group[1])
				i = 1
			else:
				all_guesses.extend(group[0][1])
		seen = []
		for metab in all_guesses:
			if metab not in seen:
				if metab[1] != 0.0:
					seen.append(metab)
		seen_score = {}
		for metab in seen:
			if metab[0] == 'nan':
				continue
			if metab[0] not in seen_score:
				seen_score[metab[0]] = metab[1]
			else:
				if metab[1] > seen_score[metab[0]]:
					seen_score[metab[0]] = metab[1]
		seen_new = []
		for k,v in seen_score.items():
			seen_new.append((k,v))
		min_guess[vote] = sorted(seen_new, reverse=True, key=lambda x:x[1])
	votes = min_guess.keys()
	votes = [ele.split('_') for ele in votes]
	vote_order = sorted([int(vote[1]) for vote in votes])
	votes = ['vote_'+str(vote) for vote in vote_order]
	seen_molecules = {}
	all_seen = []
	for vote in votes:
		seen_molecules[vote] = []
		for chem in min_guess[vote]:
			if type(chem[0]) == float:
				continue 
			if chem[0].lower() not in all_seen:
				all_seen.append(chem[0].lower())
				seen_molecules[vote].append(chem)
	return seen_molecules

def min_matched_results(matched_results):
	min_matched = {}
	for vote in matched_results:
		all_guesses = []
		number = None
		for group in matched_results[vote]:
			all_guesses.extend(group[1])
		seen = []
		for metab in all_guesses:
			if metab not in seen:
				if metab[1] != 0.0:
					seen.append(metab)
		seen_score = {}
		for metab in seen:
			if metab[0] not in seen_score:
				seen_score[metab[0]] = metab[1]
			else:
				if metab[1] > seen_score[metab[0]]:
					seen_score[metab[0]] = metab[1]
		seen_new = []
		for k,v in seen_score.items():
			seen_new.append((k,v))
		min_matched[vote] = sorted(seen_new, reverse=True, key=lambda x:x[1])
	return min_matched

def common_molecules(metfrag_data, voted_data):
	best_guess = []
	try:
		voted_data = ast.literal_eval(voted_data)
	except:
		return best_guess, [(ele, score) for ele,score in zip(metfrag_data['Name'],metfrag_data['Score'])]
	voted_metabs_names = [ele[1].lower() for ele in voted_data[1]]
	multi_voted_inchi = [ele[-1] for ele in voted_data[1]]
	common_inchi_name = []
	same_inchi = {}
	metfrag_other_guess = []
	# other parts of ele: SMILES, InChI? maybe try these if cannot find any common metabs by name 
	for ele,score, num_peaks in zip(metfrag_data['Name'],metfrag_data['Score'],metfrag_data['NoExplPeaks']):
		if num_peaks < 1:
			continue
		if str(ele).lower() in voted_metabs_names:
			best_guess.append((str(ele).lower(), score))
		else:
			metfrag_other_guess.append((str(ele).lower(),score))
	for mfrag_name, inchi_met,score, num_peaks in zip(metfrag_data['Name'],metfrag_data['InChI'],metfrag_data['Score'],metfrag_data['NoExplPeaks']):
		if num_peaks < 1:
			continue
		inchi_met_mod = '/'.join(inchi_met.split('/')[:4])
		for j,inchi_vote in enumerate(multi_voted_inchi):
			inchi_vote = ['/'.join(inchi.split('/')[:4]) for inchi in inchi_vote]
			if inchi_met_mod in inchi_vote:
				common_inchi_name.append((mfrag_name, inchi_met,score))
				best_guess.append((str(mfrag_name).lower(), score))
			else:
				if (str(mfrag_name).lower(),score) not in best_guess and (str(mfrag_name).lower(),score) not in metfrag_other_guess: # and str(mfrag_name).lower() not in voted_metabs_names:
					metfrag_other_guess.append((str(mfrag_name).lower(),score))	
	for i_1 in common_inchi_name:
		if i_1[1] not in same_inchi:
			same_inchi[i_1[1]] = [(i_1[0],i_1[2])]
		else:
			same_inchi[i_1[1]].append((i_1[0],i_1[2]))
	for inchi in same_inchi:
		same_inchi[inchi] = sorted(same_inchi[inchi], key=lambda x:x[1], reverse=True)
	reduced_names = [same_inchi[inchi][0] for inchi in same_inchi]
	for z, ele in enumerate(reduced_names):
		if type(ele[0]) == float:
			del reduced_names[z]
	best_guess = [(ele[0].lower(),ele[1]) for ele in reduced_names]
	metfrag_other_guess = [guess for guess in metfrag_other_guess if guess not in best_guess]

	return best_guess, metfrag_other_guess

def add_results_series(data_row, matched,not_matched):
	data = data_row.to_dict()
	all_votes = list(data_row.index)[4:]
	new_order = list(data_row.index)[:4]
	in_matched = [vote for vote in matched]
	in_only_met = [vote for vote in not_matched]
	not_in_matched = [vote for vote in all_votes if vote not in in_matched]
	not_in_only_met = [vote for vote in all_votes if vote not in in_only_met]
	for vote in matched:
		if matched[vote] != []:
			data['matched_'+vote] = str(matched[vote])
		else:
			data['matched_'+vote] = np.nan
	for vote in not_matched:
		if not_matched[vote] != []:
			data['only_metfrag_'+vote] = str(not_matched[vote])
		else:
			data['only_metfrag_'+vote] = np.nan
	for vote in not_in_matched:
		data['matched_'+vote] = np.nan
	for vote in not_in_only_met:
		data['only_metfrag_'+vote] = np.nan
	for vote in all_votes:
		new_order.append('matched_'+vote)
		new_order.append(vote)
		new_order.append('only_metfrag_'+vote)
	data = pd.DataFrame(data, index=[0])
	data = data[new_order]
	return data

def add_results_series_no_metfrag(data_row): #, matched,not_matched):
	data = data_row.to_dict()
	all_votes = list(data_row.index)[4:]
	new_order = list(data_row.index)[:4]
	for vote in all_votes:
		new_order.append('matched_'+vote)
		new_order.append(vote)
		new_order.append('only_metfrag_'+vote)
	for vote in all_votes:
		data['matched_'+vote] = np.nan
		data['only_metfrag_'+vote] = np.nan
	data = pd.DataFrame(data, index=[0])
	data = data[new_order]
	return data

def redunce_dict_red(matched, not_matched):
	new_not_match = {ele:[] for ele in not_matched.keys()}
	matched_list = []
	for k,v in matched.items():
		matched_list.extend(v)
	for k,v in not_matched.items():
		if len(v) > 0:
			for metab in v:
				if metab not in matched_list:
					new_not_match[k].append(metab)
	return matched, new_not_match

def best_metab_guess(df, mode):
	def make_best_guesses(new_best, guesses, no_match=False, mz_data=False):
		for ele in sorted(guesses.keys()):
			if ele in new_best and new_best[ele] != []:
				continue
			if len(guesses[ele]) != 0 and not mz_data:
				best = sorted(guesses[ele], reverse=True, key=lambda x:x[1])
				if len(best) > 1:
					top_guess = [ele for ele in best if ele[1] == best[0][1]]
				else:
					top_guess = best
				if no_match:
					top_guess = [(ele[0]+'*', ele[1]) for ele in top_guess]
				#### start here use this to come up with the new column of 'best' metab guesses to add to DF and be done
				new_best[ele] = top_guess
			else:
				new_best[ele] = []
		for priority in new_best:
			if new_best[priority] != []:
				return new_best	
		for ele in sorted(guesses.keys()):
			if len(guesses[ele]) != 0 and mz_data:					
				best = guesses[ele]
				top_guess = []
				for spec in best:
					if type(spec) == str:
						formula = spec
						continue
					for m in spec:
						if float(m[0]) <= float(best[1][0][0]):
							top_guess.append((m[1]+'**',m[0],m[2], formula))
				new_best[ele] = top_guess
		return new_best
	
	best_guess = []
	possible_col = ['vote_{}'.format(i) for i in range(int((df.shape[1]-4)/3))]
	for i, row in df.iterrows():
		if i % 200 == 0:
			print('on best metab: ', i)
		importance_order = []
		row_dict = row.to_dict()
		mz_lookup_data = []
		for possible in possible_col:
			try:
				vote_data = row_dict[possible]
			except:
				pass
			if type(vote_data) == float:
				continue
			if vote_data == str(('',[])):
				continue
			try:
				vote_data =	ast.literal_eval(vote_data)
			except:
				vote_data = vote_data
			adduct = vote_data[1][0][2]
			importance_order.append(adduct)
			mz_lookup_data.append(vote_data)
		if mode == 'negative':
			importance = []
			order = []
			guesses = {'0':[],'1':[],'2':[],'3':[]}
			guesses_mz_lookup = {'0':[],'1':[],'2':[],'3':[]}
			guesses_no_match = {'0':[],'1':[],'2':[],'3':[]}
			for adduct in importance_order:
				if adduct[:4] == 'M-H_':
					importance.append(0)
				elif 'Cl' in adduct or 'none' in adduct:
					importance.append(1)
				elif '2M' in adduct or 'Na' in adduct or 'H2O' in adduct:
					importance.append(2)
				elif 'K' in adduct:
					importance.append(3)
				else:
					importance.append(3)
			for i in range(4):
				for j, num in enumerate(importance):
					if i == num:
						if type(row_dict['matched_vote_{}'.format(j)]) != float:
							guesses[str(i)].extend(ast.literal_eval(row_dict['matched_vote_{}'.format(j)]))
						if type(row_dict['vote_{}'.format(j)]) != float:
							guesses_mz_lookup[str(i)].extend(ast.literal_eval(row_dict['vote_{}'.format(j)]))
						if type(row_dict['only_metfrag_vote_{}'.format(j)]) != float:
							guesses_no_match[str(i)].extend(ast.literal_eval(row_dict['only_metfrag_vote_{}'.format(j)]))
			new_best = {}
			new_best = make_best_guesses(new_best,guesses)
			new_best = make_best_guesses(new_best,guesses_no_match, no_match=True)
			new_best = make_best_guesses(new_best,guesses_mz_lookup, mz_data=True)	
		else:
			importance = []
			order = []
			guesses = {'0':[],'1':[],'2':[],'3':[]}
			guesses_mz_lookup = {'0':[],'1':[],'2':[],'3':[]}
			guesses_no_match = {'0':[],'1':[],'2':[],'3':[]}
			for adduct in importance_order:
				if adduct[:2] == 'H_':
					importance.append(0)
				elif 'Na' in adduct or 'none' in adduct:
					importance.append(1)
				elif 'ACN' in adduct or 'K' in adduct or 'NH4' in adduct:
					importance.append(2)
				else:
					importance.append(3)
			for i in range(4):
				for j, num in enumerate(importance):
					if i == num:
						if type(row_dict['matched_vote_{}'.format(j)]) != float:
							guesses[str(i)].extend(ast.literal_eval(row_dict['matched_vote_{}'.format(j)]))
						if type(row_dict['vote_{}'.format(j)]) != float:
							guesses_mz_lookup[str(i)].extend(ast.literal_eval(row_dict['vote_{}'.format(j)]))
						if type(row_dict['only_metfrag_vote_{}'.format(j)]) != float:
							guesses_no_match[str(i)].extend(ast.literal_eval(row_dict['only_metfrag_vote_{}'.format(j)]))
			new_best = {}
			new_best = make_best_guesses(new_best,guesses)
			new_best = make_best_guesses(new_best,guesses_no_match, no_match=True)
			new_best = make_best_guesses(new_best,guesses_mz_lookup, mz_data=True)

		added = 0
		for ele in sorted(new_best.keys()):
			if new_best[ele] == []:
				continue
			else:
				added = 1
				best_guess.append(new_best[ele])
				break
		if added == 0:
			best_guess.append(np.nan)
	df.insert(4, 'metfrag_matched_best_guess', best_guess)
	return df

def compare_results(voted_data, mass_column, rt_column, meta_frag_data):
	'''
	compare the results from the voted data to the metfrag results. 
		Inputs:
			voted_data: df the is input to program with voted order of primary mz lookups
			mass_column: from input, which column the mz data is from
			rt_column: from input, which column the rt data is from
			meta_frag_data: the dictionary of metfrag results
		Outputs;
			new_df: df from combining, has both the data when metfrag available and not. 
	'''
	### make a dict that maps the mz,rt tuple to the row index
	voted_data_dict = {}
	support_i = []
	for i, row in voted_data.iterrows():
		mass = round(row.iloc[mass_column],5)
		rt = round(row.iloc[rt_column],5)
		voted_data_dict[('{:.5f}'.format(mass),'{:.5f}'.format(rt))] = i
	### now iterate over all the metfrag results
	new_df = []
	z = 0
	seen_k = []
	for k in meta_frag_data:
		if k not in seen_k:
			seen_k.append(k)
		z+=1
		if z % 50 == 0:
			print('comparing: ',z)
		adducts_seen = {}
		### each mz,rt tuple will map to a list of both votes but the 1-3 random but different metfrags on the same vote (diff 2nd ms peaks)
		### ok need to keep track of the 0-3 metfrag runs per vote and the votes 
		seen_metfrag_summary = {}
		other_metfrag_summary = {}
		for single_guess in meta_frag_data[k]:
			if single_guess['vote'] not in adducts_seen:
				adducts_seen[single_guess['vote']] = [single_guess['adduct']]
			else:
				adducts_seen[single_guess['vote']].append(single_guess['adduct'])
			vote = single_guess['vote']
			row_index = voted_data_dict[k]
			### these are the two major things to be compared.
			add_dict = {}
			inds = voted_data.loc[row_index].index[4:]
			seen_nan = False
			for ind in inds:
				try:
					data = ast.literal_eval(voted_data.loc[row_index][ind])
					if '_'.join(data[1][0][2].split('_')[:-1]) not in add_dict:
						add_dict['_'.join(data[1][0][2].split('_')[:-1])] = [ind]
					else:
						add_dict['_'.join(data[1][0][2].split('_')[:-1])].append(ind)
				except:
					if not seen_nan:
						if single_guess['adduct'] not in add_dict:
							add_dict[single_guess['adduct']] = [ind]
						else:
							add_dict[single_guess['adduct']].append(ind)
						seen_nan = True
					else:
						pass
			votes = add_dict[single_guess['adduct']]
			for v in votes:
				vote_data = voted_data.loc[row_index][v]
				single_guess_metfrag_res = single_guess['results']
				### do the comparison
				metfrag_w_votes, other_guesses = common_molecules(single_guess_metfrag_res, vote_data)
				if v not in seen_metfrag_summary:
					seen_metfrag_summary[v] = [(single_guess['repeat'],metfrag_w_votes)]
					other_metfrag_summary[v] = [(single_guess['repeat'],other_guesses)]
				else:
					seen_metfrag_summary[v].append((single_guess['repeat'],metfrag_w_votes))
					other_metfrag_summary[v].append([(single_guess['repeat'],other_guesses)])
		min_metfrag_other_summary = min_metfrag_other(other_metfrag_summary)
		seen_metfrag_summary = min_matched_results(seen_metfrag_summary)
		seen_metfrag_summary, min_metfrag_other_summary = redunce_dict_red(seen_metfrag_summary, min_metfrag_other_summary)
		#### ok now i have the reduced list...need to make into a nice dataframe. 
		new_row = add_results_series(voted_data.iloc[row_index,:], seen_metfrag_summary,min_metfrag_other_summary)
		new_df.append(new_row)
	print('number of mz/rt tuples for ms/ms', len(seen_k))
	print('size of the actual voted data: ',voted_data.shape)
	for ele in voted_data_dict:
		if ele not in seen_k:
			new_df.append(add_results_series_no_metfrag(voted_data.iloc[voted_data_dict[ele],:]))
	new_df = pd.concat(new_df)
	new_df = new_df.reset_index()
	new_df = new_df.drop(columns='index')
	print('size of the new df: ',new_df.shape)
	return new_df

def main():
	# inputs are:
	# 1) the CSV output from the metabolite calling script (ie voted candidates)
	# 2) the path to the output of the MetF rag launching program 'parse_mgf_comb_feat_prep_metfrag.py'

	## for the first time its run: 
	# python combine_metfrag_w_votes.py -f './all_voted_metabs_longitudinal_MIT_50000noise_metabolites_ms2.csv' -c 2 -r 3 -m './ms-ms2/msms_out/' -v 'true' -n 'True' -o './all_named_metabs.csv'
	# or
	# python combine_metfrag_w_votes.py -f './pos/all_voted_metabs_mapped_all_pos.csv' -c 2 -r 3 -m './pos/msms_out/' -v 'true' -n 'True' -o './MTBLS264_all_pos.csv'
	## for the next times....
	# python combine_metfrag_w_votes.py -f './all_voted_metabs_longitudinal_MIT_50000noise_metabolites_ms2.csv' -c 2 -r 3 -m './ms-ms2/msms_out/' -v 'true'

	## for a first time with a positive:
	# python combine_metfrag_w_votes.py -f './MTBLS264/rbc_pos/all_voted_metabs_mapped_rbc_pos.csv' -c 2 -r 3 -m './MTBLS264/rbc_pos/msms_out/' -v 'true' -n 'True' -i 'positive' -o './rbc_pos.csv'

	mass_file, mass_column, rt_column, metfrag_res_path, voted, mode, first, out_file = parse_args()

	### read in the mz data:
	study_rt_metadata = pd.read_csv(mass_file)#.iloc[:,:]
	### will be used to look up the voted compounds for smiles / inchi
	study_rt_metadata = add_inchi(study_rt_metadata)
	### get the metfrag results
	metfrag_files = os.listdir(metfrag_res_path)
	metfrag_files = metfrag_files[:]

	if voted != 'true':
		print('whomp whomp, what you are trying to do is not currently supported...sorry, talk to ethan, still working on implementing this')
		exit(0)
	### parse all the metfrag files:
	print('parsing metfrag')
	if first:
		metfrag_results = parse_metfrag_output(metfrag_files, metfrag_res_path)
		pickle.dump(metfrag_results, open('./parsed_metfrag_1.pkl','wb'))
	else:
		metfrag_results = pickle.load(open('./parsed_metfrag_1.pkl', 'rb'))

	### ok now to compare the guesses with this metfrag output...
	### scan over the masses in the csv files, truncate and look up in the dict...
	print('comparing results')
	metfrag_added_data = compare_results(study_rt_metadata, mass_column, rt_column, metfrag_results)
	print('getting best metab')
	metfrag_added_data_best_guess = best_metab_guess(metfrag_added_data, mode)
	metfrag_added_data_best_guess.to_csv(out_file)

if __name__ == '__main__':
	main()
