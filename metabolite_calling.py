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
		# files = {args.f: pd.read_csv(args.f, sep='\t')}
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

def skip_ds():
	skip_ds = ['m_oxylipin_chronic_hep_b', 'm_GC_nmfi_and_bsi_diagnosis_v2_maf',
			###### the following are good after running the MTBLS352 parsing function:
           # 'DEMO_neg-norm-metaboAnalystInput_0_NGT', 'DEMO_neg-norm-metaboAnalystInput_2_Pre-DM',
           # 'DEMO_neg-norm-metaboAnalystInput_1_T2D', 'DEMO_pos-norm-metaboAnalystInput_0_NGT', 
           # 'DEMO_pos-norm-metaboAnalystInput_2_Pre-DM', 'DEMO_pos-norm-metaboAnalystInput_1_T2D', 
           'm_CER_mass_spectrometry_v4_3_CS', 'm_CER_mass_spectrometry_v4_0_NS', 
           'm_CER_mass_spectrometry_v4_2_FS', 'm_CER_mass_spectrometry_v4_1_COPD', 'm_EICO_mass_spectrometry_v4_3_CS',
           'm_EICO_mass_spectrometry_v4_3_CS','m_EICO_mass_spectrometry_v4_0_NS', 'm_EICO_mass_spectrometry_v4_2_FS',
           'm_EICO_mass_spectrometry_v4_1_COPD','m_SHOT_mass_spectrometry_v4_3_CS','m_SHOT_mass_spectrometry_v4_0_NS',
           'm_SHOT_mass_spectrometry_v4_2_FS', 'm_SHOT_mass_spectrometry_v4_1_COPD', 'm_TAG_mass_spectrometry_v4_3_CS',
           'm_TAG_mass_spectrometry_v4_0_NS', 'm_TAG_mass_spectrometry_v4_2_FS', 'm_TAG_mass_spectrometry_v4_1_COPD', 
           'm_typhoid_carriage_metabolite_profiling_mass_spectrometry_v2_maf', 'AN000452_0_Healthy',
           'AN000452_1_CRC', 'AN000452_2_Polyp', 'AN000525_1_MCD', 'AN000525_2_FSGS', 'AN000525_0_Control', 'AN000526_0_Control', 'AN000526_2_FSGS', 'AN000526_1_MCD', 
           'AN000580', 'AN000581', 'AN000582', 'AN000583', 'AN000705', 'AN000706', 'AN000931', 'AN000076_0_CN', 
           'AN000077_0_CN', 'AN000078_0_CN', 'AN000079_0_CN', 'AN000076_2_MCI', 'AN000077_2_MCI', 'AN000078_2_MCI', 'AN000079_2_MCI',
           'AN000076_1_AD', 'AN000077_1_AD', 'AN000078_1_AD', 'AN000079_1_AD',
            ###### from the one-v-one datasets
           'm_CER_mass_spectrometry_v4_COPD_FS', 'm_CER_mass_spectrometry_v4_COPD_CS', 'm_CER_mass_spectrometry_v4_COPD_NS', 'm_CER_mass_spectrometry_v4_FS_CS', 
           'm_CER_mass_spectrometry_v4_FS_NS', 'm_CER_mass_spectrometry_v4_CS_NS', 'm_EICO_mass_spectrometry_v4_COPD_FS', 'm_EICO_mass_spectrometry_v4_COPD_CS',
           'm_EICO_mass_spectrometry_v4_COPD_NS', 'm_EICO_mass_spectrometry_v4_FS_CS', 'm_EICO_mass_spectrometry_v4_FS_NS', 'm_EICO_mass_spectrometry_v4_CS_NS',
           'm_SHOT_mass_spectrometry_v4_COPD_FS', 'm_SHOT_mass_spectrometry_v4_COPD_CS', 'm_SHOT_mass_spectrometry_v4_COPD_NS', 'm_SHOT_mass_spectrometry_v4_FS_CS',
           'm_SHOT_mass_spectrometry_v4_FS_NS', 'm_SHOT_mass_spectrometry_v4_CS_NS', 'm_TAG_mass_spectrometry_v4_COPD_FS', 'm_TAG_mass_spectrometry_v4_COPD_CS',
           'm_TAG_mass_spectrometry_v4_COPD_NS', 'm_TAG_mass_spectrometry_v4_FS_CS', 'm_TAG_mass_spectrometry_v4_FS_NS', 'm_TAG_mass_spectrometry_v4_CS_NS',
           'AN000076_CN_MCI', 'AN000076_CN_AD', 'AN000076_MCI_AD', 'AN000077_CN_MCI', 'AN000077_CN_AD', 'AN000077_MCI_AD', 'AN000078_CN_MCI', 'AN000078_CN_AD', 
           'AN000078_MCI_AD', 'AN000079_CN_MCI', 'AN000079_CN_AD', 'AN000079_MCI_AD', 'AN000452_Healthy_CRC', 'AN000452_Healthy_Polyp', 'AN000452_CRC_Polyp',
           'AN000525_MCD_FSGS', 'AN000525_MCD_Control', 'AN000525_FSGS_Control', 'AN000526_MCD_FSGS', 'AN000526_MCD_Control', 'AN000526_FSGS_Control']
	low_res = ['IPO_aligned_MTBLS105_qMS', 'IPO_aligned_MTBLS105_SIM-MS', 'IPO_aligned_MTBLS315_mzData', 'AN000100', 'AN000101',
	           'XCMS-Report-annotated-SingleClass-GCTOF.', 'XCMS-Report-annotated-SingleClass-GCTOF.plasma', 'AN000615', 'IPO_aligned_ST000381_pos', 'AN000618',
	           'AN000603_plasma', 'AN000603_serum', 'AN000620_plasma', 'AN000620_serum', 'IPO_aligned_ST000385_adc2_plasma',
	           'IPO_aligned_ST000385_adc2_serum', 'IPO_aligned_ST000385_adc1_plasma', 'IPO_aligned_ST000385_adc1_serum',
	           'AN000625_GC', 'IPO_aligned_ST000388_GC', 'AN000628_plasma', 'AN000628_serum', 'IPO_aligned_ST000392_plasma',
	           'IPO_aligned_ST000392_serum', 'AN000633', 'IPO_aligned_ST000396', 'AN001390all_author', 
	           'IPO_aligned_ST000865_batch2_raw', 'IPO_aligned_ST000865_batch3_raw', 'IPO_aligned_ST000865_onebatch']
	return skip_ds + low_res 

def remove_skips(models, skip_ds):
	new_models = {}
	for k, v in models.items():
		new_models[k] = []
		for ds in v:
			if ds['data_set'] in skip_ds:
				continue
			else:
				new_models[k].append(ds)
	return new_models

def fix_mz_rt_datasets(models, files): 
	'''
	models is the pickle file for the data used to make models and all accompanying metadata
	files is a dictionary mapping the actual file name to the pd df representation of the CSV to adding the metadata to
	'''
	no_mz_rt = ['ST000918_AN001503_metabolites_all.csv', 'ST000888_AN001450_metabolites_all.csv', 'ST000763_AN001202_Healthy_PAH_metabolites_all.csv',
			    'ST000763_AN001201_Healthy_PAH_metabolites_all.csv', 'ST000608_AN000930_metabolites_all.csv', 'ST000608_AN000929_metabolites_all.csv',
			    'ST000578_AN000888_metabolites_all.csv', 'ST000578_AN000889_metabolites_all.csv', 'ST000388_AN000624_LC_metabolites_all.csv',
			    'Feng_plasmaall_author_metabolites_all.csv', 'Feng_urineall_author_metabolites_all.csv',
			    'MTBLS279_m_chronic_hep_b_NEG_metabolites_all.csv', 'MTBLS279_m_chronic_hep_b_POS_metabolites_all.csv']
	data_sets = ['AN001503', 'AN001450', 'AN001202_Healthy_PAH', 'AN001201_Healthy_PAH', 'AN000930', 'AN000929','AN000888', 'AN000889', 'AN000624_LC',
				 'plasmaall_author', 'urineall_author', 'm_chronic_hep_b_NEG', 'm_chronic_hep_b_POS']
	ds_to_csv = dict(zip(data_sets, no_mz_rt))
	mz_labels = ['mz', 'Mass-to-Charge', 'mass_to_charge', 'Mass', 'mzmed', 'm/z', 'moverz_quant',
	'quant mz', 'quantified m/z', 'row m/z', 'Quant mx', 'Quantified m/z']
	rt_labels = ['rt', 'retention_time', 'Retention Time', 'rtmed','RT', 'row retention time', 
				'retention index', 'ri', 'retention time']

	for k, v in models.items():
		for ds in v:
			if ds['data_set'] in data_sets:
				# this is the key to look up in the files dict to get the df of interest
				csv_file = ds_to_csv[ds['data_set']]
				# this is the peaks with rt/mz data to add in
				peaks_data = ds['peaks']
				# print(files[csv_file].shape, peaks_data.shape)
				for mz_header in list(peaks_data):
					if mz_header in mz_labels:
						mz_data = peaks_data[mz_header].reset_index()
						break
				for rt_header in list(peaks_data):
					if rt_header in rt_labels:
						rt_data = peaks_data[[rt_header,mz_header]].reset_index()
						break 
				try:
					files[csv_file] = pd.concat([rt_data,files[csv_file]], axis=1)
					files[csv_file].dropna(axis=1, how='all', inplace=True)
				except:
					files[csv_file] = pd.concat([mz_data,files[csv_file]], axis=1)
					files[csv_file].dropna(axis=1, how='all', inplace=True)
	return files

def share_data():
	'''
	for the sake of time, if this had the same MZ/rt features i jsut took one and processed it.
	below are the data_sets removed from the mapping but that also have unique features / model features
	'''
	repeat_mz = ['IPO_aligned_ST000763_untar_neg_Borderline Pressures_LowRisk', 'IPO_aligned_ST000763_untar_neg_Normal Pressures_LowRisk', 
			 'IPO_aligned_ST000763_untar_neg_Normal Pressures_Borderline Pressures', 'IPO_aligned_ST000763_untar_neg_PAH_LowRisk', 
			 'IPO_aligned_ST000763_untar_neg_PAH_Borderline Pressures', 'IPO_aligned_ST000763_untar_neg_PAH_Normal Pressures',
			 'IPO_aligned_ST000763_untar_neg_Healthy_LowRisk', 'IPO_aligned_ST000763_untar_neg_Healthy_Borderline Pressures',
			 'IPO_aligned_ST000763_untar_neg_Healthy_Normal Pressures',
			 'IPO_aligned_ST000763_untar_pos_Borderline Pressures_LowRisk', 'IPO_aligned_ST000763_untar_pos_Normal Pressures_LowRisk', 
			 'IPO_aligned_ST000763_untar_pos_Normal Pressures_Borderline Pressures', 'IPO_aligned_ST000763_untar_pos_PAH_LowRisk', 
			 'IPO_aligned_ST000763_untar_pos_PAH_Borderline Pressures', 'IPO_aligned_ST000763_untar_pos_PAH_Normal Pressures',
			 'IPO_aligned_ST000763_untar_pos_Healthy_LowRisk', 'IPO_aligned_ST000763_untar_pos_Healthy_Borderline Pressures',
			 'IPO_aligned_ST000763_untar_pos_Healthy_Normal Pressures',
			 'AN001202_Borderline Pressures_LowRisk', 'AN001202_Normal Pressures_LowRisk', 'AN001202_Normal Pressures_Borderline Pressures',
			 'AN001202_PAH_LowRisk', 'AN001202_PAH_Borderline Pressures', 'AN001202_PAH_Normal Pressures', 'AN001202_Healthy_LowRisk',
			 'AN001202_Healthy_Borderline Pressures', 'AN001202_Healthy_Normal Pressures', 
			 'AN001201_Borderline Pressures_LowRisk', 'AN001201_Normal Pressures_LowRisk', 'AN001201_Normal Pressures_Borderline Pressures',
			 'AN001201_PAH_LowRisk', 'AN001201_PAH_Borderline Pressures', 'AN001201_PAH_Normal Pressures', 'AN001201_Healthy_LowRisk',
			 'AN001201_Healthy_Borderline Pressures', 'AN001201_Healthy_Normal Pressures',
			 'IPO_aligned_ST000329_neg_FSGS_Control', 'IPO_aligned_ST000329_neg_MCD_Control', 'IPO_aligned_ST000329_pos_FSGS_Control',
			 'IPO_aligned_ST000329_pos_MCD_Control', 'XCMS-Report-annotated-SingleClass.27jun12_MCI_AD', 'XCMS-Report-annotated-SingleClass.27jun12_CN_AD',
			 'XCMS-Report-annotated-SingleClass.11jun12_MCI_AD', 'XCMS-Report-annotated-SingleClass.11jun12_CN_AD', 'XCMS-Report-annotated-SingleClass.04jun12_MCI_AD',
			 'XCMS-Report-annotated-SingleClass.04jun12_CN_AD', 'IPO_aligned_ST000046_20120625_pos_c18_MCI_AD', 'IPO_aligned_ST000046_20120625_pos_c18_CN_AD',
			 'IPO_aligned_ST000046_20120620_neg_c18_MCI_AD', 'IPO_aligned_ST000046_20120620_neg_c18_CN_AD', 'IPO_aligned_ST000046_20120618_pos_c18_MCI_AD',
			 'IPO_aligned_ST000046_20120618_pos_c18_CN_AD', 'IPO_aligned_ST000046_20120613_neg_hilic_MCI_AD', 'IPO_aligned_ST000046_20120613_neg_hilic_CN_AD',
			 'IPO_aligned_ST000046_20120606_neg_hilic_MCI_AD', 'IPO_aligned_ST000046_20120606_neg_hilic_CN_AD', 'IPO_aligned_ST000045_17mar_neg_II_IW',
			 'IPO_aligned_ST000045_17mar_neg_ND_IW', 'IPO_aligned_ST000045_11mar_pos_II_IW', 'IPO_aligned_ST000045_11mar_pos_ND_IW', 'IPO_aligned_ST000045_11feb_neg_II_IW',
			 'IPO_aligned_ST000045_11feb_neg_ND_IW', 'IPO_aligned_ST000045_2feb_pos_II_IW', 'IPO_aligned_ST000045_2feb_pos_ND_IW', '17March10-21-_II_IW',
			 '17March10-21-_ND_IW', '11March10-21-_II_IW', '11March10-21-_ND_IW', '11Feb10-21-r0_II_IW', '11Feb10-21-r0_ND_IW', '02Feb10-21-r0_II_IW',
			 '02Feb10-21-r0_II_IW', 'IPO_aligned_MTBLS352_neg_NGT_Pre-DM', 'IPO_aligned_MTBLS352_neg_T2D_Pre-DM', 'DEMO_pos-norm-metaboAnalystInput_NGT_Pre-DM',
			 'DEMO_pos-norm-metaboAnalystInput_T2D_Pre-DM', 'DEMO_neg-norm-metaboAnalystInput_NGT_Pre-DM', 'DEMO_neg-norm-metaboAnalystInput_T2D_Pre-DM',
			 'ulsam_IPO_reprocessedsecond_mean', 'ulsam_IPO_reprocessedfirst_mean', 'ulsam_IPO_reprocessedmean', 'ulsam_IPO_reprocessedsecond']

def clean_files(files):
	for dataset in files:
		files[dataset].drop(columns=['Unnamed: 0'], inplace=True)
		if files[dataset].iloc[:1,:4].isnull().values.any():
			# now found the first element in some datasets are NaNs and the data is shifted down 1
			# remove nans, shift data up, remove last row
			if 'MTBLS19' in dataset or 'MTBLS17' in dataset:
				files[dataset]['rt'] = files[dataset]['rt'].shift(-1)
				files[dataset]['mz'] = files[dataset]['mz'].shift(-1)
				files[dataset].drop(files[dataset].tail(1).index, inplace=True)  
	return files

def build_model_to_metab_table_dict(models, files):
	# need to make the dataset names in the models to the names of the files dict, then can mask this dict 
	# with the sig / model features. 
	# 1) get map of data_set names to name of dataset in files:
	models_to_metabs = {}
	for k,v in models.items():
		for ds in v:
			possibles = []
			# print('ON DS: ',ds['data_set'])
			match = False
			for j in files:
				if ds['peaks'].shape[0] == files[j].shape[0]:
					possibles.append(j)
			if len(possibles) == 0:
				continue
			if len(possibles) == 1:
					best = possibles[0]
			# print(possibles)
			if 'MTBLS19_IPO' in possibles[0]:
				possibles_new = [ele[8:-20] for ele in possibles]
				for poss in possibles_new:
					if poss == ds['data_set']:
						best =  possibles[0][:8]+ poss + possibles[0][-20:]
			if 'MTBLS352' in possibles[0]:
				possibles_new = [ele[9:-20] for ele in possibles]
				for poss in possibles_new:
					if poss == ds['data_set']:
						best = possibles[0][:9] + poss + possibles[0][-20:]
			if '02Feb10' in possibles[0]:
				possibles_new = [ele[9:-20] for ele in possibles]
				best = None
				for poss in possibles_new:
					if poss == ds['data_set']:
						best = possibles[0][:9] + poss + possibles[0][-20:]
				if best == None:
					best = possibles[0]
			if 'IPO_aligned_ST000392' in possibles[0]:
				possibles_new = [ele[9:-20] for ele in possibles]
				for poss in possibles_new:
					if poss == ds['data_set']:
						best = possibles[0][:9] + poss + possibles[0][-20:]
			models_to_metabs[ds['data_set']] = best
	### ok dictionary mapping data_set (ie model) to the looked up metabolites is fully functioning! 
	return models_to_metabs

def extract_sig(models):
	'''
	this makes 4 dictionaries with mapping of study name to a mask of the features for analysis - ie the stat sig features or the model features 
	'''
	stat_sig  = {}
	actual_stat_sig = {}
	model_feat = {}
	actual_model_feat = {}
	real_features = {}
	for k,v in models.items():
		for ds in v:
			# print(list(ds))
			real_features[ds['data_set']] = np.vstack((ds['pvalues'],ds['clf'].coef_)).T

			stat_sig[ds['data_set']] = ds['pvalues'] < 0.05 # this is a boolean mask for everything (ALL FEATURES)
			actual_stat_sig[ds['data_set']] = ds['pvalues'][ds['pvalues']<0.05] # this is only the stat sig values
			actual_stat_sig_model_feat = ds['clf'].coef_[stat_sig[ds['data_set']].reshape((1,-1))]
			actual_stat_sig[ds['data_set']] = np.vstack((actual_stat_sig[ds['data_set']], actual_stat_sig_model_feat)).T
			# model_feat_where_stat_sig = ds['clf'].coef_['pvalues'] < 0.05
	
			model_feat[ds['data_set']] = ds['clf'].coef_.squeeze() != 0
			actual_model_feat[ds['data_set']] = ds['clf'].coef_.squeeze()[ds['clf'].coef_.squeeze() != 0]
			actual_model_feat_stat_sig = ds['pvalues'][model_feat[ds['data_set']]].reshape((1,-1))
			actual_model_feat[ds['data_set']] = np.vstack((actual_model_feat[ds['data_set']], actual_model_feat_stat_sig)).T

	return stat_sig, model_feat, actual_stat_sig, actual_model_feat, real_features

def map_feat_to_metab_table(sig_feat, model_feat, files, mapper, actual_stat_sig,  actual_model_feat, real_feat):
	'''
	model_feat: dict mapping data_set name from models to a mask of features 
	sig_feat: dict mapping data_set name to stat sig feature mask 
	files: dict mapping arbitrary name of the looked up metabs to pd dfs of mz/rt/metabolites etc, 
	mapper: a dictionary mapping the name of the data_set to the arbitrary name of the metab lookup file it should be matched to

	returns:
		new_metab_dataset : dict mapping the data_set name to a list of 1) mask of model feature 2) sig features 3) the mz/rt/metabo name df
	'''
	new_metab_dataset = {}
	for k in mapper.keys():
		new_metab_dataset[k] = [sig_feat[k], model_feat[k], files[mapper[k]], actual_stat_sig[k], actual_model_feat[k], real_feat[k]]
	return new_metab_dataset

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
	# if 'NO' in orig_formula:
		# print('old: ', orig_formula, 'new form:', formula, 'final:', ''.join(new_form))
	return ''.join(new_form)

def find_metab_support(sig_metab_formulas, support_metab_formulas):
	count = 0
	# print(sig_metab_formulas, support_metab_formulas)
	# print('input:', sig_metab_formulas)
	# sig_metab_formulas = [make_common_formula(chem) for chem in sig_metab_formulas]
	# print('out: ',sig_metab_formulas)
	# support_metab_formulas = [make_common_formula(chem) for chem in support_metab_formulas] 
	for f_1 in support_metab_formulas:
		if f_1 in sig_metab_formulas:
			count += 1 
	if count != 0:
		return True
	else:
		return False

	# sig_formulas = extract_formula(support_metab)
	## basically take all the formulas in the sig (as a set) and see if that formula is in the support
	## if nothing in support toss it
	## if there is something, keep this metabolite, trim possible metabolites that dont fit. 

def get_adduct_data(add_data, mask):
	# print(add_data, mask)
	combined_data = []
	# print(add_data)
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
				# print('positive mode tossing because change <0: ', ele)
				continue
			if mode == 'negative' and int(ele[4]) > 0:
				# print('negative mode tossing because change >0: ', ele)
				continue
			if mode == 'positive' and int(ele[4]) > 1:
				# print('postive mode tossing because charge would be 3+ (or 2 if "none"): ', ele)
				continue
			if mode == 'negative' and int(ele[4]) < -1:
				# print('negative mode tossing because charge would be 3- (of -2 if "none"): ', ele)
				continue
			if mode == 'positive' and int(ele[4]) == 1:
				if 'none_' not in ele[3]:
					# print('positive mode tossing because charge would be 2: ', ele)
					continue
			if mode == 'negative' and int(ele[4]) == -1:
				if 'none_' not in ele[3]:
					# print('negative mode tossing because charge would be -2: ', ele)
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

def vote_with_evidence(family, metabo_num, sig_metab='none', mode='none'):
	# need to rewrite....get the data from the SIG element, then take the other metabs to dict, take the metabs and label etc 
	# that map to the mattching formula and use just this to not add additional useless data
	# currently its not getting the second or third data...start here!!!
	ms_mode = determine_polarity(list(family))
	if type(sig_metab) == int:
		sig_ps = [sig_metab]
	if  mode != 'none':
		mode_values = list(family[mode].values)
		sig_metab = [round(sig_metab[0],7)]
		mode_values = [round(p,7) for p in mode_values]
		# ok this gets the data for the proper compound (important because a family could have 2 sig family members and dont want process teh same one twice and
		# this has it match up with the j indexing which is critical for writing)
		sig_ps = [i for i, ele in enumerate(mode_values) if ele == sig_metab]
	primary_or_support = {}

	for z in range(family.shape[0]):
		if z == sig_ps[0]:
			primary_or_support[z] = 'primary_{}'.format(metabo_num)
		else:
			primary_or_support[z] = 'support_{}'.format(metabo_num)
	family_data = {}
	if ms_mode == 'positive':
		data_dict = {}
		# extract the data from each of the rows, 
		for i in range(family.shape[0]):
			family_row = family.iloc[[i]]
			metabolite_guess = vote_positive(family_row, ms_mode)
			family_data[i] = metabolite_guess
		for j in range(family.shape[0]):
			new_metabolites = []
			support = family_data[j]
			metab_2_forms = [make_common_formula(ele[0]) for ele in support]
			good_metabs_1 = []
			for metabolites_1 in family_data[sig_ps[0]]:
				if make_common_formula(metabolites_1[0]) not in metab_2_forms:
					new_metabolites.append(tuple(('',[])))
				else:
					for metabolites_2 in support:
						if make_common_formula(metabolites_2[0]) == make_common_formula(metabolites_1[0]):
							new_metabolites.append(metabolites_2)
			data_dict[j] = new_metabolites
	else:
		data_dict = {}
		for i in range(family.shape[0]):
			family_row = family.iloc[[i]]
			metabolite_guess = vote_negative(family_row,ms_mode)
			family_data[i] = metabolite_guess
		for j in range(family.shape[0]):
			new_metabolites = []
			support = family_data[j]
			metab_2_forms = [make_common_formula(ele[0]) for ele in support]
			good_metabs_1 = []
			for metabolites_1 in family_data[sig_ps[0]]:
				if make_common_formula(metabolites_1[0]) not in metab_2_forms:
					new_metabolites.append(tuple(('',[])))
				else:
					for metabolites_2 in support:
						if make_common_formula(metabolites_2[0]) == make_common_formula(metabolites_1[0]):
							new_metabolites.append(metabolites_2)
			data_dict[j] = new_metabolites
	return data_dict, primary_or_support

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
						###### 3-29-19: basically the secondary and primary thing sucked...so this removes it, no duplicates now
						# if mode == 'all':
						# 	for k in range(family.shape[0]):
						# 		if list(family[col_header].iloc[[k]]) == list(sig_metabs[col_header].iloc[[j]]):
						# 			sig_metab = k
						# 	vote_data_dict, primary_or_support = vote_with_evidence(family, j, sig_metab)
						# else:
						# 	sig_metab = sig_metabs[mode].iloc[[j]].values
						# 	vote_data_dict, primary_or_support = vote_with_evidence(family, j, sig_metab, mode)
						# for i in vote_data_dict.keys():
						# 	feature_type.append(primary_or_support[i])
						# 	all_feature_data.append(vote_data_dict[i])
						# new_header = [ele for ele in list(family) if 'meta' not in ele and 'hmdb' not in ele and 'lipid' not in ele and 'chebi' not in ele]
						# meta_data = family[new_header]
						# all_meta_data.append(meta_data)
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

def mask_file_metab_dict(metabs_with_masks, rt_metadata):
	for ds in metabs_with_masks:
		# print(ds)
		data = metabs_with_masks[ds]
		data[2][['p_values','model_feat']] = data[5]
		#### what the data is: 
		#### FIRST i change the p-values and model features to match the correct data set due to sharing. 
		#### data = [sig_feat[k], model_feat[k], files[mapper[k]], actual_stat_sig[k], actual_model_feat[k]]
		#### if just doing the significant or model features:
		# if ds != 'serum_IPO_aligned_Feng_serum_batch1':
		# 	continue
		#### get the rt time window that depends on the dataset...
		if type(rt_metadata) != float:
			rt = rt_metadata[rt_metadata.analysis==ds]['delta_time'].values[0]
			if math.isnan(rt ):
				rt = 2.0
		else:
			rt = rt_metadata

		#### mask by the significant metabolites indicies or by the model features!
		sig_metabs = data[2][data[0]].copy().reset_index() 
		model_metabs = data[2][data[1]].copy().reset_index()
		#### set the real data_set specific p-values - due to sharing of datasets for some studies
		#### NOT EVEN SURE I NEED THIS ANYMORE SINCE I JUST CHANGE ALL THE DATA above...
		sig_metabs[['p_values', 'model_feat']] = data[3]
		model_metabs[['model_feat', 'p_values']] = data[4]
		for mode, metabs in zip(['p_values', 'model_feat'], [sig_metabs,model_metabs]): # at a later point ru with 'all' in mode and 'data[2]' in the dataset
			ds_voting_protocol(metabs, data, mode, rt, ds)

def mask_file_metab_dict_parallel(metabs_with_masks, rt_metadata):
	data_map = []
	for ds in metabs_with_masks:
		data = metabs_with_masks[ds]
		data[2][['p_values','model_feat']] = data[5]
		if type(rt_metadata) != float:
			rt = rt_metadata[rt_metadata.analysis==ds]['delta_time'].values[0]
			if math.isnan(rt):
				rt = 2.0
		else:
			rt = rt_metadata
		#### mask by the significant metabolites indicies or by the model features!
		sig_metabs = data[2][data[0]].copy().reset_index() 
		model_metabs = data[2][data[1]].copy().reset_index()
		#### set the real data_set specific p-values - due to sharing of datasets for some studies
		#### NOT EVEN SURE I NEED THIS ANYMORE SINCE I JUST CHANGE ALL THE DATA above...
		sig_metabs[['p_values', 'model_feat']] = data[3]
		model_metabs[['model_feat', 'p_values']] = data[4]
		data_map.append((sig_metabs, data, 'p_values', rt, ds))
		data_map.append((model_metabs, data, 'model_feat', rt, ds))
	pool = Pool(os.cpu_count())
	# data_map = [(sig_metabs,data,mode,rt,ds)]
	pool.starmap(ds_voting_protocol, data_map)

def non_metaanalysis_voting(files, rt_metadata, out_path):
	for study in files:
		# print(list(files[study].loc[0]))
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

	#way to run this program:
	#python ./metabolite_calling.py -t 5 -r 10 -p "../30avg_log_reg_YES_bn_ds_models_and_sigfeat_NO_log_NO_standscal_YES_ovo.pkl" -n 'true'
	# or if you already made the pickle:
	# python ./metabolite_calling.py -t 5 -r 10 -p "../30avg_log_reg_YES_bn_ds_models_and_sigfeat_NO_log_NO_standscal_YES_ovo.pkl" -f "./charged_smiles_formula_metab_lookup.pkl" -n 'true'
	#### NOTE 2-22-19: reran the full this 
	## thus made a new 'changed...csv' file

	# or for a new one off file or directory:
	# with an output path
	# python ./metabolite_calling.py -t 5 -r 10 -f "./underworlds/longitudinal_MIT_50000noise_metabolites_ms2.csv" -s ',' -o './underworlds/'
	# or once run:
	# python ./metabolite_calling.py -t 5 -r 10 -f "./underworlds/underworlds_charged_smiles_formula_metab_lookup.pkl" -s ',' -o './underworlds/'


	# does two analyses: 1) based on significant feature from the model. 2) based on stat sig features
	# 'files' is a dictionary mapping file name for metabolite lookup to a pd df of that csv file!
	# 'models' is a dict mapping study to list of different datasets / models - ie the main chunck of the project that has most info
	ppm, rt_window, models, files, initial, meta_analysis,out_path = parse_args()
	# if this is the first time running, just run on the desired file or directory and this will make a 
	# pickle of the same data but with the added info (formula / smiles) and save it...
	# WHAT OUT: change the file name below to not overwrite important things! 
	if initial:
		files = add_inchi_smiles_formula(files)
		if meta_analysis:
			pickle.dump(files, open('./charged_smiles_formula_metab_lookup.pkl', 'wb'))
		if not meta_analysis:
			pickle.dump(files, open('./underworlds/underworlds_charged_smiles_formula_metab_lookup.pkl', 'wb'))

	if meta_analysis:
		study_rt_metadata = pd.read_csv('./ms_instrument_column_polarity_dataset_names.csv', sep='\t')[['Accession', 'analysis', 'delta_time']]
	# remove low res datasets of datasets without the mz/rt data...
		ds_to_remove = skip_ds()
		models = remove_skips(models, ds_to_remove)

		# go in and fix the data sets lacking the rt / mz info even if its there
		# datasets without the mz/rt columns added...
		files = fix_mz_rt_datasets(models, files)

		# clean up the dataframe to remove the excess rt/ mz columns that are empty:
		files = clean_files(files)

		#build the mapper of model names to names used for metabolite mapping 
		model_to_metab_mapper = build_model_to_metab_table_dict(models, files)

		# Extract the significant features / model features and use this to make a new DF of just these
		stat_sig, model_feat, actual_stat_sig, actual_model_feat, real_features = extract_sig(models)

		# need to apply the stat / model features to the data frames and replace the values since they are just for one dataset and not necessaryily the mapped one
		metabs_with_masks = map_feat_to_metab_table(stat_sig, model_feat, files, model_to_metab_mapper, actual_stat_sig, actual_model_feat, real_features)

		# what can happen now:
		# 1) extract the features for sig / model and simply use  - for this still need a way to get the correct metabolite name or top 3 picks
		### how do deal with metabs in different dbs that are the same but different names...maybe use smiles instead? do they all have them?
		# for each feature, get a cluster family with a user defined RT cluster
		# mask_file_metab_dict(metabs_with_masks, study_rt_metadata)

		###### When ready for full run use the parallel version:
		mask_file_metab_dict_parallel(metabs_with_masks, study_rt_metadata)
	
	else:
		study_rt_metadata = rt_window
		# need to make the metabs_with_masks some how...
		non_metaanalysis_voting(files, study_rt_metadata, out_path)
	# open analysis:
	# use clary's data of known compounds to see how often a na adduct also has a H adduct-- do we even need to worry aoubt these adducts 
	# take clary's list of annotated things, and clean things up and see how many features we can explain with known chemicals
	# are we seeing conserted changes within mz's, if there is a signal it is probably conserved. 
if __name__ == '__main__':
	main()
