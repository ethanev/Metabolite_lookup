# Metabolite_lookup - 4-18-19: can now use MS2 info with metfrag to perform higher confidence assignment
Project for taking raw m/z features from high resulution mass spec studies and giving putative empirical formulas with associated chemical names, ppm error for various ion adducts. The lookup happens in four databases: HMDB, Lipidmap, ChEBI and MetaCyc.
How to run each script should have an example or multiple commented out in the main() function. 

Run order:
0) run XCMS to get your feature table and extract the MS2 spectra as a .mgf
1) use parsing_metabolite_dbs.py 
	this will read in and process the various databases - NOTE: you need to download the appropriate databases! I just wasn't sure if I was allowed to put them up here. sorry!

2) use mapping_mz_to_metabolites.py
	this will perform the mapping from mz features to small molecule names in the databases. You will mostly want to use the options for csv file inputs. It is build for a meta-analysis I have been working on but I made sure its compatible with single files that you want to know what metabolites are there.
	this program makes csv files for each datasets and outputs them to a user specified directory. 
  
3) use metabolite_calling.py (on the output file from mapping_mz_to_metabolites.py)
	this will go throguh the csv file(s) from the last program (enter their location), cluster by defined RT for peak groupings (user set), find similar chemicals and call the best (if possible) metabolite for the sig / model features.

MS2 part
4) run parse_mgf_comb_feat_prep_metfrag.py (on output of metabolite_calling.py) - highly suggest running on a instance of cluster with plenty of CPUs and memory (it launches a lot of metfrag programs). 
	* requirements: the command line version of metfrag and hmb_2017-07-23.csv, kegg_2017-07-23.csv, lipidmaps.csv all in the same directory as this script and the data input. 
	* output requirements: make a dir called 'msms_out' prior to running ... have not added the feature to auto do this yet sorry. 

5) run combine_metfrag_w_votes.py 



	
	
