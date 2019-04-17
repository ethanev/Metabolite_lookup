import argparse
import pickle
import pandas as pd
from multiprocessing import Pool
import os
import numpy as np

#not used anymore: 
def get_rt_cluster(rt,mask, mz_data, rt_data, time_window=0.1):
	if rt:
		spec_rt_list = list(rt_data[mask].values)
		rt_and_neighbor_mask, sec_ind_to_prim_ind, family_names = find_rt_neighbors(spec_rt_list, rt_data, time_window)
		return rt_and_neighbor_mask, sec_ind_to_prim_ind, family_names
		# write something to cluster and return a boolean mask as list of [primary, sec, sec, primary....]
	else:
		num_metab = mz_data[mask].shape[0]
		return mask, {}, {str(mz_data.iloc[i]):['primary_{}'.format(str(i))] for i in range(num_metab)}

def find_rt_neighbors(rt_list, rt_data, time_window):
	'''
	rt_list - list of the rt times for the primary sig / model features
	rt_data - the full dataframe of rt_values
	time_window - the delta time for which two things are part of the same cluster 
	'''
	# family_dict - maps rt time of the primary feature to family members. 
	family_dict = {}
	# secondary_to_primary - maps each features to which primary clusters it belongs to. 
	secondary_to_primary = {str(i):[] for i in list(range(rt_data.shape[0]))}
	# new_mask - a boolean mask that expant the mask from the sig/model features with their cluster family
	new_mask = np.zeros(rt_data.shape, dtype=bool)
	for j, rt in enumerate(rt_list):
		family_dict[str(rt)] = []
		for i, rt_ in enumerate(list(rt_data.values)):
			if abs(float(rt) - float(rt_)) < time_window:
				new_mask[i] = True
				secondary_to_primary[str(i)].append(j)
				family_dict[str(rt)].append((i,str(rt_)))
	secondary_to_primary = trim_dict(secondary_to_primary)
	return new_mask, secondary_to_primary, family_dict

def trim_dict(dict_with_empty):
	return {ele:dict_with_empty[ele] for ele in dict_with_empty if len(dict_with_empty[ele])!=0}

def mz_fam_from_rt_fam(mz_data, sec_ind_to_prim_ind, original_mask, mz_to_lookup):
	mz_to_fam = {str(ele):[] for ele in list(mz_to_lookup.values)}
	primary_mz = list(mz_data[original_mask].values)
	if len(mz_to_fam) == 0:
		return {}
	else:
		for i, ele in enumerate(primary_mz):
			mz_to_fam[str(ele)].append('primary_{}'.format(str(i)))
	for ele in sec_ind_to_prim_ind:
		# is this index accuratly getting the mass? 
		mz = mz_data.iloc[int(ele)]
		primary_links = sec_ind_to_prim_ind[ele]
		for link in primary_links:
			if 'primary_{}'.format(str(link)) in mz_to_fam[str(mz)]:
				pass
			else:
				mz_to_fam[str(mz)].append('secondary_{}'.format(str(link)))
	return mz_to_fam

# functions only for the meta-analysis
def get_all_ds_mz(datasets, features_to_use='all'):
	'''
	make a dictionary mapping study to list of dicts where each dict maps data set name to mz features of interest

	Basically if the dataset does not have mz data its in the skip_ds list, if its a low resolution study its in the low_res 
		NEITHER ARE USED
	'''
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
	mz_data_dict = {}
	for k,v in datasets.items():
		mz_data_dict[k] = []
		for ds in v:
			if ds['data_set'] in skip_ds:
				continue
			if ds['data_set'] in low_res:
				continue
			if ds['data_set'] in repeat_mz:
				continue
			### commented out the next line since i want the program to do lookups for all, then be able to alter clustering in a seperate program
			# mz_data,mz_to_families = get_data_from_ds(ds, features_to_use)
			# mz_data_dict[k].append([ds['data_set'],mz_data, mz_to_families])

			#### the following line are for when running on all mz values
			mz_data, extra_data = get_data_from_ds(ds, features_to_use)
			mz_data_dict[k].append([ds['data_set'], mz_data, extra_data])
	return mz_data_dict

def get_data_from_ds(dataset, features_to_use='all'):
	'''
	given a single dataset, this will return the desired features from the 'peaks' list. 
	'''
	peaks_data = dataset['peaks']
	mz_labels = ['mz', 'Mass-to-Charge', 'mass_to_charge', 'Mass', 'mzmed', 'm/z', 'moverz_quant',
	'quant mz', 'quantified m/z', 'row m/z', 'Quant mx', 'Quantified m/z']
	rt_labels = ['rt', 'retention_time', 'Retention Time', 'rtmed','RT', 'row retention time', 
				'retention index', 'ri', 'retention time']
	extra_data = pd.DataFrame({'p_values':dataset['pvalues'], 'model_feat':dataset['clf'].coef_.squeeze()})
	for mz_header in list(peaks_data):
		if mz_header in mz_labels:
			mz_data = peaks_data[mz_header]
			break
	for rt_header in list(peaks_data):
		if rt_header in rt_labels:
			rt_data = peaks_data[[rt_header,mz_header]]
			if rt_data[rt_header].isnull().any():
				pass
			else:
				extra_data = extra_data.join(rt_data)
			break
	############ THE IF...ELIF part is not used!!!###### just run on all data, process later 
	# 1) if you want to use the features that are significant, mask the peaks list by the list of sig feats
	if features_to_use == 'sig':
		mask = dataset['pvalues'] < 0.05
	# 2) if you want to use the features from the model, mask by non-zero model features
	elif features_to_use == 'model':
		mask = dataset['clf'].coef_.squeeze() != 0
	# 3) otherwise, use all the features and look all of them up.
	else:
		return mz_data, extra_data
	# rt_and_neighbor_mask, sec_ind_to_prim_ind, family_names = get_rt_cluster(rt, mask, mz_data, rt_data)
	# if rt:
		# mz_to_family = mz_fam_from_rt_fam(mz_data,sec_ind_to_prim_ind, mask, mz_data[rt_and_neighbor_mask])
	# else:
		# mz_to_family = family_names
	# print('shape just sig: ',mz_data[mask].shape, 'shape post new mask: ', mz_data[rt_and_neighbor_mask].shape, 'total shape: ', mz_data.shape)
	# return mz_data[rt_and_neighbor_mask], mz_to_family

def parse_MTBLS352(datasets):
	'''
	Goal: fix the peaks data for one specific dataset!
	Input: the full dataset object 
	'''
	for ds in datasets['MTBLS352']:
		if 'DEMO' in ds['data_set']:
			data = list(ds['peaks'])
			rt = [ele.split('_')[0] for ele in data]
			mz = [ele.split('_')[1][:-3] for ele in data]
			peaks = pd.DataFrame({'rt': rt, 'mz':mz})
			ds['peaks'] = peaks
	return datasets

def process_data_parallel(study, dataset, dbs, ppm, metadata, features_to_use='all'):
	'''
	In:
		study - the string name of the study, ie MTBLS266 etc
		dataset - list of 0) the dataset name, 1) the mz_data used for metabolite lookup and
			3) the extra data, df with p_values, model_sig_feat and rt/mz if the rt is available
		dbs - dictionary of the 4 databases used for lookup
		ppm - the ppm setting used for the lookup
		metadata - the metadata file used to get the mode of ms used
		features_to_use - defaulting to 'all' since this should just map all possible metabolites
	Out: 
		a lot of csv files for all the lookups
	'''
	out_file = './'+features_to_use+'/'+study+'_'+dataset[0]+'_metabolites_'+features_to_use+'.csv'
	mode = metadata.loc[dataset[0]]['mode']
	extra_data = dataset[2]
	mz_data_neutralized = neutralize_masses(dataset[1], mode)
	looked_up_metabolites = map_db_adducts_lookup(mz_data_neutralized, dbs, ppm)
	looked_up_metabolites = pd.concat([extra_data,looked_up_metabolites], axis=1)
	looked_up_metabolites.to_csv(out_file, sep='\t')
	print('finished {}'.format(dataset[0]))

#functions for general use
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
			# if db[compound]['mass'] == None or float(db[compound]['mass'])==0.0:
				# print(correctedompound, db[compound]['mass'])
			# print(compound, db[compound])
			# try:
			ppm_err = abs((float(mass) - float(db[compound]['mass'])))/float(db[compound]['mass']) * 1000000
			# print(mass, db[compound]['mass'], ppm_err)
			# except:
				# pass
				# print('error', compound, db[compound])
			if ppm_err <= ppm:
				# print(compound, ppm_err, db[compound])
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
	### example: 
	### in the parallel for of the metaanalysis version this writes csv files for all useful datasets to the 'all' folder
	# python ./mapping_mz_to_metabolites.py -f '../30avg_log_reg_YES_bn_ds_models_and_sigfeat_NO_log_NO_standscal_YES_ovo.pkl' -p 10 -w 'sig'
	# for say a single file:
	# python ./mapping_mz_to_metabolites.py -f './underworlds/alm_mix_plus_controls.csv' -p 10 -c 1 -m 'negative' -o './underworlds/comp_mapped_almix_controls.csv' 
	# python ./mapping_mz_to_metabolites.py -f './underworlds/post_filt_samples_blanks_controls.csv' -p 10 -c 1 -m 'negative' -o './underworlds/longitudinal_MIT_50000noise_metabolites.csv'
	mz_file, csv_col, database, ppm, mode, study, features_to_use, out_file = parse_args()
	pool = Pool(os.cpu_count())
	# get the databases you want to do the look up from
	meta_db, chebi_db, hmdb_db, lipidmap_db = get_db()
	dbs = {'meta':meta_db, 'chebi':chebi_db, 'hmdb':hmdb_db, 'lipid':lipidmap_db}
	# next line is kinda a bad assumption...
	try:
		metadata = pd.read_csv('../ms_instrument_column_polarity_dataset_names.csv', sep='\t').set_index('analysis')
	except:
		print('no metadata file available, not using')

	# read in the masses to look up - if you want to run on all files in the meta-analysis input a pickle of all the models above (output of feature_analysis.ipynb)
	if mz_file[-3:] == 'pkl':
		# use a pickle with models and sig feats and 'peaks' for full utility
		pickled_data = pickle.load(open(mz_file, 'rb'))
		pickled_data = parse_MTBLS352(pickled_data) # just fixed the data for this one study
		#get dict mapping study to list of datasets mapping ds_name to list of feat
		all_mz_data = get_all_ds_mz(pickled_data, features_to_use)		

		data_map = []
		for study in all_mz_data:
			for dataset in all_mz_data[study]:
				data_map.append((study, dataset, dbs, ppm, metadata,'all'))
		pool.starmap(process_data_parallel, data_map)
			#### for one at a time processing: 
			# for dataset in all_mz_data[study]:
				# out_file = './'+features_to_use+'/'+study+'_'+dataset[0]+'_metabolites_'+features_to_use+'.csv'
				# mode = metadata.loc[dataset[0]]['mode']
				# mz_data_neutralized = neutralize_masses(dataset[1], mode)
				# looked_up_metabolites = map_db_adducts_lookup(mz_data_neutralized, dbs, ppm)
				# looked_up_metabolites.to_csv(out_file, sep='\t')

	else:
		if mz_file[-3:] == 'csv':
			data = pd.read_csv(mz_file)
			data = data.iloc[:,:]
			mz_data = data.iloc[:,csv_col]
		elif mz_file[-3:] == 'tsv':
			data = pd.read_csv(mz_file, sep='\t')
			mz_data = data.iloc[:,csv_col]


		### based on LC type / column used (ie if its say it could be a polar metab or lipid but you ran c18 its probably lipid)
	
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