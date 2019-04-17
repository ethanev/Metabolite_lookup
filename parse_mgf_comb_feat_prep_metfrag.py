import pandas as pd
import numpy as np
import argparse
import ast
import random
import os
from multiprocessing import Pool

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('-f', help='csv files from the metabolite_calling.py or csv with list of mzs of iterest')
	parser.add_argument('-c', help='which column is the mass in the csv?', default=0)
	parser.add_argument('-r', help='which column is the rt in the csv ', default=0)
	parser.add_argument('-m', help='mgf file to parse with secondary mass spectra')
	parser.add_argument('-v', help='is this with voted data or just simple mz csv, this assumes it is, enter anything to alter this', default='true')
	parser.add_argument('-p', help='is it positive mode? "True" for pos', default='False')
	args = parser.parse_args()
	mass_file = args.f
	voted = args.v
	is_pos = args.p
	mass_column = int(args.c)
	rt_column = int(args.r)
	mgf_file = args.m
	return mass_file, mass_column, rt_column, mgf_file, voted, is_pos

def parse_mgf(mgf_file):
	mfg_data = {}
	with open(mgf_file, 'r') as f:
		begin = False
		for line in f:
			if 'PEAK_ID' in line:
				begin = True
				new_entry = ''
				continue
			line = line.strip().split()
			if 'END' in line:
				begin = False
				mfg_data[mass][rt].append(new_entry)
				continue
			if line == []:
				continue
			if 'TITLE=msLevel' in line:
				rt = line[3][:-1]
				mass = line[9][:-1]
				if mass not in mfg_data:
					mfg_data[mass] = {rt:[]}
				else:
					if rt not in mfg_data[mass]:
						mfg_data[mass].update({rt:[]})
					else:
						pass
			if begin:
				assert len(line)>0
				new_entry += line[0]
				new_entry += '\t'
				new_entry += line[1]
				new_entry += '\n'
	return mfg_data

def write_peak_list_file(mass, rt, mgf_dict, out_path='./', mass_tol = 0.003, rt_tol=5.0):
	possible_mass_data = []
	for ele in mgf_dict:
		if float(ele)-mass_tol <= float(mass) <= float(ele)+mass_tol:
			for rt_d in mgf_dict[ele]:
				if float(rt_d)-rt_tol <= float(rt) <= float(rt_d)+rt_tol:
					possible_mass_data.append(mgf_dict[ele][rt_d])
	mz_data_file = str(mass)+'_'+str(rt)+'.txt'
	files = []
	if len(possible_mass_data) == 0:
		return 0
	#### to use only 3 max ms2's since the results tend to be the same...use the followig:
	if len(possible_mass_data) > 3:
		rc_1 = random.randint(0,len(possible_mass_data)-1)
		rc_2 = random.randint(0,len(possible_mass_data)-1)
		while rc_2 == rc_1:
			rc_2 = random.randint(0,len(possible_mass_data)-1)
		rc_3 = random.randint(0,len(possible_mass_data)-1)
		while rc_3 == rc_2 or rc_3 == rc_1:
			rc_3 = random.randint(0,len(possible_mass_data)-1)
		for j in [rc_1,rc_2,rc_3]:
			with open(out_path+str(j)+'_'+mz_data_file, 'w') as f:
				files.append(str(j)+'_'+mz_data_file)
				f.write(possible_mass_data[j][0])
	else:
		for j in range(len(possible_mass_data)):
			with open(out_path+str(j)+'_'+mz_data_file, 'w') as f:
				files.append(str(j)+'_'+mz_data_file)
				f.write(possible_mass_data[j][0])		
	# to use all possible ms2s....use the following:
	# for j in range(len(possible_mass_data)):
	# 	with open(out_path+str(j)+'_'+mz_data_file, 'w') as f:
	# 		files.append(str(j)+'_'+mz_data_file)
	# 		f.write(possible_mass_data[j][0])
	return files

def neutralize_mass(mass, mass_row, is_pos):
	### NOTE: not actually running metfrag on any of the predicted C13 peaks since these wont reliably match to anything...
	adducts = {'H':-1.007276, 'Na': -22.989218, 'K': -38.963158, 'M-H-H2O': 19.01839,
			'NH4': -18.033823, 'M-H': 1.007276, 'Cl': -34.969402, 'none':0.0, 
			'M_ACN_H':-42.034276, '2M-H':1.007276,'2M_H': -1.007276, '2M_Na':-22.989218, 
			'2M_ACN_H':-42.034276, 'M-2H_Na':-20.974666, 'M-2H_K':-36.948606 } 
			#'C13-H':1.007276, 'C13_H':-1.007276, '2C13_H':-1.007276, '2C13-H':1.007276, '3C13_H':-1.007276, '3C13-H':1.007276,
	try:	
		vote_1_data = ast.literal_eval(mass_row)
		adduct = vote_1_data[1][0][2] #[1] get the list of compounds, [0] get the first compound of the list [2] get the adduct of that comp
		if is_pos != 'True':
			adduct = adduct.split('_')[:-1]
			if len(adduct) > 1:
				adduct = '_'.join(adduct)
			else:
				adduct = adduct[0]
		if is_pos == 'True':
			adduct = '_'.join(adduct.split('_')[:-1])
		if '2M' in adduct:
			neut_mass = float((float(mass)+adducts[adduct])/2)
		else:
			neut_mass = mass + adducts[adduct]
	except:
		if is_pos == 'True':
			adduct = 'H'
		else:
			adduct = 'M-H'
		neut_mass = mass + adducts[adduct]
	if adduct == 'H' or adduct =='M-H':
		h_seen = True
	else:
		h_seen = False
	return neut_mass, h_seen

def write_prm_file(mz_data_files, neut_mass, orig_mass, rt, vote, db, is_postive, add, out_path='./msms_out'):
	all_prm_files = []
	for j, mz_data_file in enumerate(mz_data_files):
		with open('./{}_{}_orig_mz_{:.5f}_neut_mz_{:.5f}_db_{}.prm'.format(j,vote,orig_mass,neut_mass,db), 'w') as f:
			f.write('PeakListPath = {} \n'.format(mz_data_file))
			if db == 'KEGG':
				f.write('MetFragDatabaseType = KEGG \n')
			else:
				f.write('MetFragDatabaseType = LocalCSV \n')
				if db == 'hmdb':
					f.write('LocalDatabasePath = ./hmdb_2017-07-23.csv \n')
				if db == 'lipidmaps':
					f.write('LocalDatabasePath = ./lipidmaps.csv \n')
			f.write('NeutralPrecursorMass = {} \n'.format(neut_mass))
			f.write('DatabaseSearchRelativeMassDeviation = 5 \n')
			f.write('FragmentPeakMatchAbsoluteMassDeviation = 0.002 \n')
			f.write('FragmentPeakMatchRelativeMassDeviation = 5 \n')
			if is_postive == 'True':
				f.write('PrecursorIonMode = {} \n'.format(add))
			else:
				f.write('PrecursorIonMode = {} \n'.format(add))
			f.write('IsPositiveIonMode = {} \n'.format(is_postive))
			f.write('MetFragCandidateWriter = CSV \n')
			f.write('SampleName = {}_{}_neut_mz_{:.5f}_orig_mass_{:.5f}_rt_{:.5f}_db_{} \n'.format(j,vote,neut_mass,orig_mass,rt,db))
			f.write('ResultsPath = {} \n'.format(out_path))
			f.write('MaximumTreeDepth = 2 \n')
			f.write('MetFragPreProcessingCandidateFilter = UnconnectedCompoundFilter \n')
			f.write('MetFragPostProcessingCandidateFilter = InChIKeyFilter \n')
		all_prm_files.append('{}_{}_orig_mz_{:.5f}_neut_mz_{:.5f}_db_{}.prm'.format(j,vote,orig_mass,neut_mass,db))
	return all_prm_files

def write_bash(file_list, out_path='./'):
	with open(out_path+'launch_metfrags.sh', 'w') as f:
		for fi in file_list:
			f.write('java -jar MetFrag2.4.5-CL.jar {} \n'.format(fi))

def build_parallel_metfrag(file_list):
	return [('java -jar ./MetFrag2.4.5-CL.jar ./{}'.format(fi)) for fi in file_list]

def run_parallel(met_frag_arg):
	os.system(met_frag_arg)

def main():
	##### NOTE: run this in the same folder as the MetFrag program and the databases!!! #####
	# inputs are:
	# 1) the CSV output from the metabolite calling scri
	# 2) the mgf file with all the secondardy mass spectra for the dataset, this is parsed first
	### Example command line:  FOR NEGATIVE
	# python parse_mgf_comb_feat_prep_metfrag.py -f './all_voted_metabs_longitudinal_MIT_50000noise_metabolites_ms2.csv' -c 2 -r 3 -m './MIT_long_50000noise_MS-MS.mgf' -v 'true'
	#### FOR POSTIVE:
	# python parse_mgf_comb_feat_prep_metfrag.py -f './all_voted_metabs_mapped_plasma_pos.csv' -c 2 -r 3 -m './MTBLS264_plasma_pos.mgf' -v 'true' -p 'True'

	mass_file, mass_column, rt_column, mgf_file, voted, is_pos = parse_args()
	pool = Pool(os.cpu_count()-1)
	# extract the ms2 data from the mgf file
	mgf_dict = parse_mgf(mgf_file)

	# read in the mz data:
	study_rt_metadata = pd.read_csv(mass_file).iloc[:,:]
	# get just the mass data: 
	mass_data = study_rt_metadata.iloc[:,mass_column]
	rt_data = list(study_rt_metadata.iloc[:,rt_column])
	metab_types = list(study_rt_metadata.iloc[:,1])
	file_list = []
	
	# ok now basically for each interesting mass, look to see if it has associated ms2 data in mgf_dict
	for i, (mass,rt,metab_type) in enumerate(zip(mass_data,rt_data,metab_types)):
		# if 'support' in metab_type:
		# 	continue
		### write the data.txt file (just be a print out of the mgf_dict - tab spacing between mz / intensity)
		sec_ms_file_names = write_peak_list_file(mass, rt, mgf_dict)
		### basically if there are no secondary ms for the given mz/rt move on to the next...
		if sec_ms_file_names == 0:
			continue
		### if this data didn't come from my voting script sorry lol....run it first
		if voted != 'true':
			print('what you are trying to do is not currently supported...sorry, talk to ethan')
			exit(0)
			# need to do something special since all i will have are mz's
		# neutralize the mass based on the results of the voting...
		mass_row = study_rt_metadata.loc[i]
		votes = list(mass_row.index)[4:]
		no_masses = False
		h_seen = False
		for vote in votes:
			v = mass_row[vote]
			if is_pos == 'True':
				add = '1'
			else:
				add = '-1'
			try:
				tuple_row = ast.literal_eval(v)
				# print(tuple_row[1][0][2])
				if 'C13' in tuple_row[1][0][2]:
					continue
				elif 'H-H2O' in tuple_row[1][0][2]:
					continue
				elif '2M' in tuple_row[1][0][2]:
					continue
				elif 'M-2H' in tuple_row[1][0][2]:
					continue
				elif 'Cl' in tuple_row[1][0][2]:
					add = '35'
				elif '-H_' in tuple_row[1][0][2]:
					add = '-1'
				elif 'none' in tuple_row[1][0][2]:
					add = '0'
				elif 'ACN_H' in tuple_row[1][0][2]:
					add = '42'
				elif 'M_H' in tuple_row[1][0][2]:
					add = '1'
				elif 'NH4' in tuple_row[1][0][2]:
					add = '18'
				elif 'Na' in tuple_row[1][0][2]:
					add = '23'
				elif 'K' in tuple_row[1][0][2]:
					add = '39'	
				# print(tuple_row[1][0][2], add)							
			except:
				pass
			# print(v, 'here', add)
			# ok use the following for those there they does even have a 1st vote
			if pd.isna(v):
				# print(v, 'breaking after h+ / H-')
				if h_seen:
					# print('already saw a H+ or H-, not duplicating')
					break
				# this means no compounds at all...we will just run one H+ or H- 
				no_masses = True
			if v == str(('',[])):
				# print(v, vote, votes[-1])
				## This means there is probably some compounds somewhere...we will wait to run the H+ or H- till all votes are in and we dont see it
				if vote == votes[-1]:
					if not h_seen:
						# print('last vote!')
						neut_mass, h_seen = neutralize_mass(mass, v, is_pos)
				else:
					continue

			neut_mass, h_seen = neutralize_mass(mass, v, is_pos)
			##### write the param.prm file 
			for db in ['KEGG', 'hmdb', 'lipidmaps']:
				files = write_prm_file(sec_ms_file_names, neut_mass, mass, rt, vote, db, is_pos, add)
				file_list += files

			if no_masses:
				break
	###### OK HERE IS WHERE YOU NEED TO MULTIPROCESS...
	###### LAUNCH METFRAGS!
	arg_list = build_parallel_metfrag(file_list)
	pool.map(run_parallel,arg_list)

		##### or to just write out....
		#### write_bash(file_list)		

if __name__ == '__main__':
	main()