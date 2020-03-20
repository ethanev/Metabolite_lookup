#!/usr/bin/env python
# Author: Ethan D. Evans, 2019 [eevans@mit.edu]

import argparse
import pickle
import pandas as pd
from multiprocessing import Pool
import os
import numpy as np

def get_db():
	'''
	Read in the assocaited database file
	'''
	meta_db = pickle.load(open('./metacyc_metabolites.pkl', 'rb'))
	chebi_db = pickle.load(open('./chebi_metabolites.pkl', 'rb'))
	hmdb_db = pickle.load(open('./hmdb_metabolites.pkl', 'rb'))
	lipidmap_db = pickle.load(open('./lipidmap_metabolites.pkl', 'rb'))
	return meta_db, chebi_db, hmdb_db, lipidmap_db

def neutralize_masses(mass_data, mode):
	'''
	Will make a dictionary of modes mapping to lists of masses
	NOTE: -H means it lost a H, _H means it got the H...
	'''
	adducts = { 'H':-1.007276, 'Na': -22.989218, 'K': -38.963158, 'M-H-H2O': 19.01839,
				'NH4': -18.033823, 'M-H': 1.007276, 'Cl': -34.969402, 'none':0.0, 
				'C13-H':0.007276,'2C13-H':-0.99943353,'3C13-H':-1.992724, #using a value of 1 for a 13C since the masses actually are over a window  in the processed data and not single peak data
				'C13_H':-2.007276, '2C13_H':-3.007276,'3C13_H':-4.007276,  
				'M_ACN_H':-42.034276, '2M-H':1.007276,'2M_H': -1.007276, '2M_Na':-22.989218, 
				'2M_ACN_H':-42.034276, 'M-2H_Na':-20.974666, 'M-2H_K':-36.948606 } 
				#  highly prob neg mode are 2M-h and M-H-h2o
	if mode == 'positive':
		possible_adduct_mz = {'H':[], 'Na':[], 'K':[], 'NH4':[], 'none':[], 'C13_H':[], '2C13_H':[], '3C13_H':[], 'M_ACN_H':[], '2M_H':[], '2M_Na':[], '2M_ACN_H':[]}
	elif mode == 'negative':
		possible_adduct_mz = {'M-H':[], 'Cl':[], 'M-H-H2O':[], 'none':[], 'C13-H':[], '2C13-H':[], '3C13-H':[], '2M-H':[], 'M-2H_Na':[], 'M-2H_K':[]} # maybe add in +Formic acid?
	else:
		possible_adduct_mz = {'H':[], 'Na':[], 'K':[], 'NH4':[], 'M-H':[], 'Cl':[], 'none':[], 'C13-H':[], 'C13_H':[], '2C13_H':[], '3C13_H':[], 'M_ACN_H':[], '2M_H':[], '2M_Na':[], '2M_ACN_H':[], '2C13-H':[], '3C13-H':[], '2M-H':[], 'M-2H_Na':[], 'M-2H_K':[]}

	for mz in mass_data:
		for adduct in possible_adduct_mz.keys():
			if '2M' in adduct:
				possible_adduct_mz[adduct].append(float((float(mz)+adducts[adduct])/2))
			else:
				possible_adduct_mz[adduct].append(float(mz)+adducts[adduct])
	return possible_adduct_mz

def mz_db_lookup(mz_list, db, ppm):
	all_mz_data = []
	for mass in mz_list:
		mz_compounds = []
		for compound in db:
			ppm_err = abs((float(mass) - float(db[compound]['mass'])))/float(db[compound]['mass']) * 1000000
			if ppm_err <= ppm:
				mz_compounds.append('%s, %.3f'%(compound, ppm_err))
		all_mz_data.append(mz_compounds)
	return all_mz_data

def map_db_adducts_lookup(mz_dict, db_dict, ppm):
	'''
	input:
		mz_dict: dictionary of adduct mapping to list of masses corrected for the adduct
		db_dict: dictionary of databases
		ppm: ppm error max for compound ID
	'''
	mapped_metabolites = {}
	mapped_metabolites['mz'] = mz_dict['none']
	for db in db_dict:
		for adduct in mz_dict:
			add_db = adduct+'_'+db
			mapped_metabolites[add_db] = mz_db_lookup(mz_dict[adduct], db_dict[db], ppm)
	return pd.DataFrame(mapped_metabolites)

def parse_args():
	
	parser = argparse.ArgumentParser()
	parser.add_argument('-f', help='input file of masses - Note there should be a column header! (csv and tsv supported, pickled model/sigfeat')
	parser.add_argument('-w', help='which features to use? use: "sig", "model" or "all"', default='all')
	parser.add_argument('-p', help='ppm error acceptable in compound matching, "positive" or "negative" give you specific results \
								otherwise it will give you both positive and negative')
	# the following are less critical and more special use (ie if you just input a single csv/tsv, narrow which database etc)
	parser.add_argument('-c', help='which column number in the csv do you want to use? USE 0 INDEXIN - only needed if inputing a csv or tsv', default=0)
	parser.add_argument('-d', help='database you want to search against - defaults to all 4', default='all')
	parser.add_argument('-m', help='ion mode used in the mass spec', default='none')
	parser.add_argument('-o', help='output file name', default='./feat_to_metab.csv')
	parser.add_argument('-s', help='if using pickled data, which study? ie MTBLS17 etc', default=False)
	args = parser.parse_args()
	mz_file = str(args.f)
	csv_col = int(args.c)
	database = args.d
	ppm = float(args.p)
	study = args.s
	features_to_use = args.w
	mode = args.m
	out_file = args.o
	return mz_file, csv_col, database, ppm, mode, study, features_to_use, out_file

def main():
	### examples for say a single file:
	# python ./mapping_mz_to_metabolites.py -f './underworlds/alm_mix_plus_controls.csv' -p 10 -c 1 -m 'negative' -o './underworlds/comp_mapped_almix_controls.csv' 
	# python ./mapping_mz_to_metabolites.py -f './underworlds/post_filt_samples_blanks_controls.csv' -p 10 -c 1 -m 'negative' -o './underworlds/longitudinal_MIT_50000noise_metabolites.csv'

	mz_file, csv_col, database, ppm, mode, study, features_to_use, out_file = parse_args()
	# get the databases you want to do the look up from
	meta_db, chebi_db, hmdb_db, lipidmap_db = get_db()
	dbs = {'meta':meta_db, 'chebi':chebi_db, 'hmdb':hmdb_db, 'lipid':lipidmap_db}
	# read in the masses to look up - if you want to run on all files in the meta-analysis input a pickle of all the models above (output of feature_analysis.ipynb)
	if mz_file[-3:] == 'csv':
		data = pd.read_csv(mz_file)
		data = data.iloc[:,:]
		mz_data = data.iloc[:,csv_col]
	elif mz_file[-3:] == 'tsv':
		data = pd.read_csv(mz_file, sep='\t')
		mz_data = data.iloc[:,csv_col]
	# get neutral masses for all mz values
		# ASSUMPTION: all mz values come from singly ionized species! It's not a great assumption. 
		# gives a dictionary of adducts mapping to the associated masses
	mz_data_neutralized = neutralize_masses(mz_data, mode)
	# look up masses in each of the databases
	# multiple looks ups: all 4 databases with all adducts. 
	looked_up_metabolites = map_db_adducts_lookup(mz_data_neutralized, dbs, ppm)
	columns = list(data)
	rt_label = ['rt', 'retention_time', 'Retention Time', 'rtmed','RT', 'row retention time', 'retention index', 'ri', 'retention time']
	for col in columns[1:]:
		if col in rt_label:
			looked_up_metabolites['rt'] = data[col]
			break
	looked_up_metabolites.to_csv(out_file)


if __name__ == '__main__':
	main()
