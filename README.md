# Metabolite_lookup
Project for taking raw m/z features from high resulution mass spec studies and giving putative empirical formulas with associated chemical names, ppm error for various ion adducts. The lookup happens in four databases: HMDB, Lipidmap, ChEBI and MetaCyc.
Run order: 
1) use parsing_metabolite_dbs.py 
	this will read in and process the various databases

2) use mapping_mz_to_metabolites.py
	this will perform the mapping from mz features to small molecule names in the databases. You will mostly want to use the options for csv file inputs. It is build for a meta-analysis I have been working on but I made sure its compatible with single files that you want to know what metabolites are there.
	this program makes csv files for each datasets and outputs them to a user specified directory. 
  
3) use metabolite_calling.py
	this will go throguh the csv file(s) from the last program (enter their location), cluster by defined RT for peak groupings (user set), find similar chemicals and call the best (if possible) metabolite for the sig / model features.
