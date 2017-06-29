# -*- coding: utf-8 -*-

####################################################################################
### Author: Patrick Tran Van (patrick.tranvan@gmail.com)
###
### This script is useful for running PAML's codeml with RAxML's tree.
###
### Need these folders and files :
###
### alignments/ : contain the multiple alignments NON formatted in phylip format.
### branch_lengths/ : contain the RAxML's tree NON formatted in newick format.
### codeml.ctl : control template.
###
### At the end, a .sh file is created for VITAL-IT infrastructure. 
###
### How to use:
###
### python timema2codeml.py -s timema -i1 alignments -i2 branch_lengths -i3 codeml_timema.ctl -e <your_email> -o paml.sh
###
### python timema2codeml.py -s table_3rate -i1 output/ -o 3rate_table.txt
### python timema2codeml.py -s table_2rate -i1 output/ -o 2rate_table.txt
### python timema2codeml.py -s table_free -i1 output/ -o free_table.txt
### 
####################################################################################

import os
import re
import argparse
import sys
import glob
import os.path
from Bio import SeqIO
import shutil

def alignments_format(alignment_before, alignment_after):
	"""
	Add spaces after header (required for codeml). 
	"""
	
	alignment_before_file = open(alignment_before,"r")
	content = alignment_before_file.readlines()
	
	alignment_after_file = open(alignment_after,"a")			
	alignment_after_file.write(content[0])

	specie_row = 1
	
	while specie_row < len(content):
		alignment_after_file.write(content[specie_row].split()[0] + "   " + content[specie_row].split()[1] + "\n")
		specie_row+=1
	
	alignment_before_file.close()
	alignment_after_file.close()


def codeml_control(organism, model, sample_name, treefile, outfile, dest):
	"""
	Create a codeml's control file for each job. 
	"""

	if organism == "timema":
	
		if model == "1":	# free model
			replacements = {"alignment.phy":sample_name, "timema.tree":treefile, "paml.out":outfile, "model_number":str(1)}
		else:
			replacements = {"alignment.phy":sample_name, "timema.tree":treefile, "paml.out":outfile, "model_number":str(2)}
		
		with open("codeml_timema.ctl") as infile, open(dest, 'w') as outfile:
		    for line in infile:
		        for src, target in replacements.iteritems():
		            line = line.replace(src, target)
		        outfile.write(line)
		        	        
def model3rate(raxml_before, raxml_after):
	"""
	Format to 3rate model. 
	"""
	
	rate3_before_file = open(raxml_before,"r")
	
	# Read the newick file
	
	for raxml_tree in rate3_before_file:
			
		split_raxml = raxml_tree.split(",")
		
		complete_format = []

		# The branch length are multiplied by 3 and branch label is added (sex or asex)
		
		part0 = split_raxml[0].split(":")[0] + ":" + str(round(float(split_raxml[0].split(":")[1])*3,9)) + "#1"
		part1 = split_raxml[1].split(":")[0] + ":" + str(round(float(split_raxml[1].split(":")[1].split(")")[0])*3,9)) + "#2" + "):" + str(round(float(split_raxml[1].split(":")[2])*3,9))
		part2 = split_raxml[2].split(":")[0] + ":" + str(round(float(split_raxml[2].split(":")[1])*3,9)) + "#1"
		part3 = split_raxml[3].split(":")[0] + ":" + str(round(float(split_raxml[3].split(":")[1].split(")")[0])*3,9)) + "#2" + "):" + str(round(float(split_raxml[3].split(":")[2].split(")")[0])*3,9)) + "):" + str(round(float(split_raxml[3].split(":")[3])*3,9))
		part4 = split_raxml[4].split(":")[0] + ":" + str(round(float(split_raxml[4].split(":")[1])*3,9)) + "#1"
		part5 = split_raxml[5].split(":")[0] + ":" + str(round(float(split_raxml[5].split(":")[1].split(")")[0])*3,9)) + "#2" + "):" + str(round(float(split_raxml[5].split(":")[2].split(")")[0])*3,9)) + "):" + str(round(float(split_raxml[5].split(":")[3])*3,9))
		part6 = split_raxml[6].split(":")[0] + ":" + str(round(float(split_raxml[6].split(":")[1])*3,9)) + "#1"
		part7 = split_raxml[7].split(":")[0] + ":" + str(round(float(split_raxml[7].split(":")[1].split(")")[0])*3,9)) + "#2" + "):" + str(round(float(split_raxml[7].split(":")[2].split(")")[0])*3,9)) + "):" + str(round(float(split_raxml[7].split(":")[3])*3,9))
		part8 = split_raxml[8].split(":")[0] + ":" + str(round(float(split_raxml[8].split(":")[1])*3,9)) + "#1"
		part9 = split_raxml[9].split(":")[0] + ":" + str(round(float(split_raxml[9].split(":")[1].split(")")[0])*3,9)) + "#2):" + split_raxml[9].split(":")[2]
		
		complete_format.append(part0)
		complete_format.append(part1)
		complete_format.append(part2)
		complete_format.append(part3)
		complete_format.append(part4)
		complete_format.append(part5)
		complete_format.append(part6)
		complete_format.append(part7)
		complete_format.append(part8)
		complete_format.append(part9)
		
		new_tree = ','.join(complete_format)

		rate3_before_file = open(raxml_after,"w")		
		rate3_before_file.write(new_tree)					
		rate3_before_file.close()

	rate3_before_file.close()

def model2rate(raxml_before, raxml_after):
	"""
	Format to 2rate model. 
	"""
	
	rate2_before_file = open(raxml_before,"r")
	
	# Read the newick file
	
	for raxml_tree in rate2_before_file:
			
		split_raxml = raxml_tree.split(",")
		
		complete_format = []

		# The branch length are multiplied by 3 and branch label is added (sex + asex)
		
		part0 = split_raxml[0].split(":")[0] + ":" + str(round(float(split_raxml[0].split(":")[1])*3,9)) + "#1"
		part1 = split_raxml[1].split(":")[0] + ":" + str(round(float(split_raxml[1].split(":")[1].split(")")[0])*3,9)) + "#1" + "):" + str(round(float(split_raxml[1].split(":")[2])*3,9))
		part2 = split_raxml[2].split(":")[0] + ":" + str(round(float(split_raxml[2].split(":")[1])*3,9)) + "#1"
		part3 = split_raxml[3].split(":")[0] + ":" + str(round(float(split_raxml[3].split(":")[1].split(")")[0])*3,9)) + "#1" + "):" + str(round(float(split_raxml[3].split(":")[2].split(")")[0])*3,9)) + "):" + str(round(float(split_raxml[3].split(":")[3])*3,9))
		part4 = split_raxml[4].split(":")[0] + ":" + str(round(float(split_raxml[4].split(":")[1])*3,9)) + "#1"
		part5 = split_raxml[5].split(":")[0] + ":" + str(round(float(split_raxml[5].split(":")[1].split(")")[0])*3,9)) + "#1" + "):" + str(round(float(split_raxml[5].split(":")[2].split(")")[0])*3,9)) + "):" + str(round(float(split_raxml[5].split(":")[3])*3,9))
		part6 = split_raxml[6].split(":")[0] + ":" + str(round(float(split_raxml[6].split(":")[1])*3,9)) + "#1"
		part7 = split_raxml[7].split(":")[0] + ":" + str(round(float(split_raxml[7].split(":")[1].split(")")[0])*3,9)) + "#1" + "):" + str(round(float(split_raxml[7].split(":")[2].split(")")[0])*3,9)) + "):" + str(round(float(split_raxml[7].split(":")[3])*3,9))
		part8 = split_raxml[8].split(":")[0] + ":" + str(round(float(split_raxml[8].split(":")[1])*3,9)) + "#1"
		part9 = split_raxml[9].split(":")[0] + ":" + str(round(float(split_raxml[9].split(":")[1].split(")")[0])*3,9)) + "#1):" + split_raxml[9].split(":")[2]
		
		complete_format.append(part0)
		complete_format.append(part1)
		complete_format.append(part2)
		complete_format.append(part3)
		complete_format.append(part4)
		complete_format.append(part5)
		complete_format.append(part6)
		complete_format.append(part7)
		complete_format.append(part8)
		complete_format.append(part9)

		new_tree = ','.join(complete_format)

		rate2_before_file = open(raxml_after,"w")		
		rate2_before_file.write(new_tree)					
		rate2_before_file.close()

	rate2_before_file.close()
	

def modelFree(raxml_before, raxml_after):
	"""
	Format to free model. 
	"""
	
	free_before_file = open(raxml_before,"r")
	
	# Read the newick file
	
	for raxml_tree in free_before_file:
			
		split_raxml = raxml_tree.split(",")
		
		complete_format = []

		# The branch length are multiplied by 3
				
		part0 = split_raxml[0].split(":")[0] + ":" + str(round(float(split_raxml[0].split(":")[1])*3,9))
		part1 = split_raxml[1].split(":")[0] + ":" + str(round(float(split_raxml[1].split(":")[1].split(")")[0])*3,9)) + "):" + str(round(float(split_raxml[1].split(":")[2])*3,9))
		part2 = split_raxml[2].split(":")[0] + ":" + str(round(float(split_raxml[2].split(":")[1])*3,9))
		part3 = split_raxml[3].split(":")[0] + ":" + str(round(float(split_raxml[3].split(":")[1].split(")")[0])*3,9)) + "):" + str(round(float(split_raxml[3].split(":")[2].split(")")[0])*3,9)) + "):" + str(round(float(split_raxml[3].split(":")[3])*3,9))
		part4 = split_raxml[4].split(":")[0] + ":" + str(round(float(split_raxml[4].split(":")[1])*3,9))
		part5 = split_raxml[5].split(":")[0] + ":" + str(round(float(split_raxml[5].split(":")[1].split(")")[0])*3,9)) + "):" + str(round(float(split_raxml[5].split(":")[2].split(")")[0])*3,9)) + "):" + str(round(float(split_raxml[5].split(":")[3])*3,9))
		part6 = split_raxml[6].split(":")[0] + ":" + str(round(float(split_raxml[6].split(":")[1])*3,9))
		part7 = split_raxml[7].split(":")[0] + ":" + str(round(float(split_raxml[7].split(":")[1].split(")")[0])*3,9)) + "):" + str(round(float(split_raxml[7].split(":")[2].split(")")[0])*3,9)) + "):" + str(round(float(split_raxml[7].split(":")[3])*3,9))
		part8 = split_raxml[8].split(":")[0] + ":" + str(round(float(split_raxml[8].split(":")[1])*3,9))
		part9 = split_raxml[9].split(":")[0] + ":" + str(round(float(split_raxml[9].split(":")[1].split(")")[0])*3,9)) + "):" + split_raxml[9].split(":")[2]
		
		complete_format.append(part0)
		complete_format.append(part1)
		complete_format.append(part2)
		complete_format.append(part3)
		complete_format.append(part4)
		complete_format.append(part5)
		complete_format.append(part6)
		complete_format.append(part7)
		complete_format.append(part8)
		complete_format.append(part9)

		new_tree = ','.join(complete_format)

		free_after_file = open(raxml_after,"w")		
		free_after_file.write(new_tree)					
		free_after_file.close()

	free_before_file.close()
		
	
def timema(alignment_directory, branch_lengths, control_file, email, output_file):
	"""
	Creating script for the 3rate, 2rate and free models. 
	"""
	
	folder = "output"
	
	if os.path.isfile(output_file):
			os.remove(output_file)
	
	script = open(output_file,"a")
	script.write("#!/bin/bash\n\n")
	script.close()
			
	# Reset the result folder if existing
					
	if os.path.isdir(folder):
		shutil.rmtree(folder, ignore_errors=True)
	else:
		os.mkdir(folder)

	# Look for every alignment file and create all the files needed for codeml

	alignments_repo = os.path.join(alignment_directory, "*")
	
	alignment_file = 0
	while alignment_file < len(glob.glob(alignments_repo)):
	
		sample = sorted(glob.glob(alignments_repo))[alignment_file]	# alignments file
		sample_name = sample.split("/")[-1]
		
		sample_repo = os.path.join(folder, sample_name)
		model3rate_repo = os.path.join(sample_repo, "model3rate")
		model2rate_repo = os.path.join(sample_repo, "model2rate")
		modelfree_repo = os.path.join(sample_repo, "modelFree")
		
		# Create a directory for each sample and model
	
		os.makedirs(sample_repo)
		
		os.makedirs(model3rate_repo)
		alignments_format(sample, os.path.join(model3rate_repo, sample_name))	# alignment re-format
			
		codeml_control("timema", "2", sample_name, "RAxML_result." + sample_name, sample_name + "3rate", os.path.join(model3rate_repo, "codeml.ctl"))	# control creation		
		model3rate(os.path.join(branch_lengths, "RAxML_result." + sample_name + ".results"), os.path.join(model3rate_repo, "RAxML_result." + sample_name))
		

		os.makedirs(model2rate_repo)
		alignments_format(sample, os.path.join(model2rate_repo, sample_name))	# alignment re-format
			
		codeml_control("timema", "2", sample_name, "RAxML_result." + sample_name, sample_name + "2rate", os.path.join(model2rate_repo, "codeml.ctl"))	# control creation		
		model2rate(os.path.join(branch_lengths, "RAxML_result." + sample_name + ".results"), os.path.join(model2rate_repo, "RAxML_result." + sample_name))

		os.makedirs(modelfree_repo)
		alignments_format(sample, os.path.join(modelfree_repo, sample_name))	# alignment re-format
				
		codeml_control("timema", "1", sample_name, "RAxML_result." + sample_name, sample_name + "Free", os.path.join(modelfree_repo, "codeml.ctl"))	# control creation
		modelFree(os.path.join(branch_lengths, "RAxML_result." + sample_name + ".results"), os.path.join(modelfree_repo, "RAxML_result." + sample_name))

		
		# Writing a script for VITAL-IT.
		
		script = open(output_file,"a")
				
		script.write("echo 'cd {} && module add Phylogeny/paml/4.9a && codeml' | bsub -u {} -N -J paml_{}\n\n".format(model3rate_repo, email, sample_name + "_3rate"))
		script.write("echo 'cd {} && module add Phylogeny/paml/4.9a && codeml' | bsub -u {} -N -J paml_{}\n\n".format(model2rate_repo, email, sample_name + "_2rate"))
		script.write("echo 'cd {} && module add Phylogeny/paml/4.9a && codeml' | bsub -u {} -N -J paml_{}\n\n".format(modelfree_repo, email, sample_name + "_Free"))
							
		script.close()				
			
		alignment_file += 1


def table_3rate(output_directory, output_file):
	"""
	Final table for the 3rate model. 
	"""
		
	if os.path.isfile(output_file):
			os.remove(output_file)
	
	script = open(output_file,"a")
	script.write("OG_ID\tlnL\tw0\tw1\tw2\tdS_Tdi\tdS_Tps\tdS_Tsi\tdS_Tcm\tdS_Tms\tdS_Tce\tdS_Tge\tdS_Tpa\tdS_Tte\tdS_Tbi\tdN_Tdi\tdN_Tps\tdN_Tsi\tdN_Tcm\tdN_Tms\tdN_Tce\tdN_Tge\tdN_Tpa\tdN_Tte\tdN_Tbi\n")
	script.close()
								
	data_repo = os.path.join(output_directory, "*")
	
	sample_file = 0
	while sample_file < len(glob.glob(data_repo)):
		
		tab_line = []
		
		sample = sorted(glob.glob(data_repo))[sample_file]
		sample_name = sample.split("/")[-1].split(".")[0]		
		rate3_file_path = os.path.join(sample, "model3rate/{}.phy3rate".format(sample_name))
				
		codeml_output = open(rate3_file_path,"r")
		lnl_info = ([line.split() for line in codeml_output if re.search("lnL",line)][0])
		codeml_output.close()

		codeml_output = open(rate3_file_path,"r")
		w_info = ([line.split() for line in codeml_output if re.search("w",line)][1])
		codeml_output.close()
			
		codeml_output = open(rate3_file_path,"r")
		ds_info = ([next(codeml_output).split() for line in codeml_output if re.search("dS tree",line)][0])
		codeml_output.close()

		codeml_output = open(rate3_file_path,"r")
		dn_info = ([next(codeml_output).split() for line in codeml_output if re.search("dN tree",line)][0])
		codeml_output.close()

								
		tab_line.append(sample_name)
		tab_line.append(lnl_info[4])
		tab_line.append(w_info[4])
		tab_line.append(w_info[5])
		tab_line.append(w_info[6])
		tab_line.append(ds_info[1].split(",")[0])
		tab_line.append(ds_info[3].split(")")[0])
		tab_line.append(ds_info[6].split(",")[0])
		tab_line.append(ds_info[8].split(")")[0])
		tab_line.append(ds_info[12].split(",")[0])
		tab_line.append(ds_info[14].split(")")[0])
		tab_line.append(ds_info[18].split(",")[0])
		tab_line.append(ds_info[20].split(")")[0])
		tab_line.append(ds_info[24].split(",")[0])
		tab_line.append(ds_info[26].split(")")[0])
		tab_line.append(dn_info[1].split(",")[0])
		tab_line.append(dn_info[3].split(")")[0])
		tab_line.append(dn_info[6].split(",")[0])
		tab_line.append(dn_info[8].split(")")[0])
		tab_line.append(dn_info[12].split(",")[0])
		tab_line.append(dn_info[14].split(")")[0])
		tab_line.append(dn_info[18].split(",")[0])
		tab_line.append(dn_info[20].split(")")[0])
		tab_line.append(dn_info[24].split(",")[0])
		tab_line.append(dn_info[26].split(")")[0])
		
			 	
	 	script = open(output_file,"a")
		script.write('\t'.join(tab_line)+"\n")
		script.close()
	
		sample_file += 1


def table_2rate(output_directory, output_file):
	"""
	Final table for the 2rate model. 
	"""
		
	if os.path.isfile(output_file):
			os.remove(output_file)
	
	script = open(output_file,"a")
	script.write("OG_ID\tlnL\tw0\tw1\tdS_Tdi\tdS_Tps\tdS_Tsi\tdS_Tcm\tdS_Tms\tdS_Tce\tdS_Tge\tdS_Tpa\tdS_Tte\tdS_Tbi\tdN_Tdi\tdN_Tps\tdN_Tsi\tdN_Tcm\tdN_Tms\tdN_Tce\tdN_Tge\tdN_Tpa\tdN_Tte\tdN_Tbi\n")
	script.close()
								
	data_repo = os.path.join(output_directory, "*")
	
	sample_file = 0
	while sample_file < len(glob.glob(data_repo)):
		
		tab_line = []
		
		sample = sorted(glob.glob(data_repo))[sample_file]
		sample_name = sample.split("/")[-1].split(".")[0]		
		rate2_file_path = os.path.join(sample, "model2rate/{}.phy2rate".format(sample_name))
				
		codeml_output = open(rate2_file_path,"r")
		lnl_info = ([line.split() for line in codeml_output if re.search("lnL",line)][0])
		codeml_output.close()

		codeml_output = open(rate2_file_path,"r")
		w_info = ([line.split() for line in codeml_output if re.search("w",line)][1])
		codeml_output.close()
				
		codeml_output = open(rate2_file_path,"r")
		ds_info = ([next(codeml_output).split() for line in codeml_output if re.search("dS tree",line)][0])
		codeml_output.close()

		codeml_output = open(rate2_file_path,"r")
		dn_info = ([next(codeml_output).split() for line in codeml_output if re.search("dN tree",line)][0])
		codeml_output.close()
					
		tab_line.append(sample_name)
		tab_line.append(lnl_info[4])
		tab_line.append(w_info[4])
		tab_line.append(w_info[5])
		tab_line.append(ds_info[1].split(",")[0])
		tab_line.append(ds_info[3].split(")")[0])
		tab_line.append(ds_info[6].split(",")[0])
		tab_line.append(ds_info[8].split(")")[0])
		tab_line.append(ds_info[12].split(",")[0])
		tab_line.append(ds_info[14].split(")")[0])
		tab_line.append(ds_info[18].split(",")[0])
		tab_line.append(ds_info[20].split(")")[0])
		tab_line.append(ds_info[24].split(",")[0])
		tab_line.append(ds_info[26].split(")")[0])
		tab_line.append(dn_info[1].split(",")[0])
		tab_line.append(dn_info[3].split(")")[0])
		tab_line.append(dn_info[6].split(",")[0])
		tab_line.append(dn_info[8].split(")")[0])
		tab_line.append(dn_info[12].split(",")[0])
		tab_line.append(dn_info[14].split(")")[0])
		tab_line.append(dn_info[18].split(",")[0])
		tab_line.append(dn_info[20].split(")")[0])
		tab_line.append(dn_info[24].split(",")[0])
		tab_line.append(dn_info[26].split(")")[0])
					 	
	 	script = open(output_file,"a")
		script.write('\t'.join(tab_line)+"\n")
		script.close()

		sample_file += 1


def table_free(output_directory, output_file):
	"""
	Final table for the free model. 
	"""
		
	if os.path.isfile(output_file):
			os.remove(output_file)
	
	script = open(output_file,"a")
	script.write("OG_ID\tlnL\tdS_Tdi\tdS_Tps\tdS_Tsi\tdS_Tcm\tdS_Tms\tdS_Tce\tdS_Tge\tdS_Tpa\tdS_Tte\tdS_Tbi\tdN_Tdi\tdN_Tps\tdN_Tsi\tdN_Tcm\tdN_Tms\tdN_Tce\tdN_Tge\tdN_Tpa\tdN_Tte\tdN_Tbi\twtr_Tdi\twtr_Tps\twtr_Tsi\twtr_Tcm\twtr_Tms\twtr_Tce\twtr_Tge\twtr_Tpa\twtr_Tte\twtr_Tbi\n")
	script.close()
								
	data_repo = os.path.join(output_directory, "*")
	
	sample_file = 0
	while sample_file < len(glob.glob(data_repo)):
		
		tab_line = []
		
		sample = sorted(glob.glob(data_repo))[sample_file]
		sample_name = sample.split("/")[-1].split(".")[0]		
		free_file_path = os.path.join(sample, "modelFree/{}.phyFree".format(sample_name))
				
		codeml_output = open(free_file_path,"r")
		lnl_info = ([line.split() for line in codeml_output if re.search("lnL",line)][0])
		codeml_output.close()

		codeml_output = open(free_file_path,"r")
		ds_info = ([next(codeml_output).split() for line in codeml_output if re.search("dS tree",line)][0])
		codeml_output.close()

		codeml_output = open(free_file_path,"r")
		dn_info = ([next(codeml_output).split() for line in codeml_output if re.search("dN tree",line)][0])
		codeml_output.close()

		codeml_output = open(free_file_path,"r")
		w_info = ([next(codeml_output).split() for line in codeml_output if re.search("w ratios as labels for TreeView",line)][0])
		codeml_output.close()
		
		index_w_info = [1,4,9,12,19,22,29,32,39,42]
						
		tab_line.append(sample_name)
		tab_line.append(lnl_info[4])
		tab_line.append(ds_info[1].split(",")[0])
		tab_line.append(ds_info[3].split(")")[0])
		tab_line.append(ds_info[6].split(",")[0])
		tab_line.append(ds_info[8].split(")")[0])
		tab_line.append(ds_info[12].split(",")[0])
		tab_line.append(ds_info[14].split(")")[0])
		tab_line.append(ds_info[18].split(",")[0])
		tab_line.append(ds_info[20].split(")")[0])
		tab_line.append(ds_info[24].split(",")[0])
		tab_line.append(ds_info[26].split(")")[0])
		tab_line.append(dn_info[1].split(",")[0])
		tab_line.append(dn_info[3].split(")")[0])
		tab_line.append(dn_info[6].split(",")[0])
		tab_line.append(dn_info[8].split(")")[0])
		tab_line.append(dn_info[12].split(",")[0])
		tab_line.append(dn_info[14].split(")")[0])
		tab_line.append(dn_info[18].split(",")[0])
		tab_line.append(dn_info[20].split(")")[0])
		tab_line.append(dn_info[24].split(",")[0])
		tab_line.append(dn_info[26].split(")")[0])
		
		for index in index_w_info:
			tab_line.append(w_info[index].split("#")[-1])
			 	
	 	script = open(output_file,"a")
		script.write('\t'.join(tab_line)+"\n")
		script.close()
	
		sample_file += 1


##############
# Main script
						
def main(argv):
	
	mod=[]
	mod.append('\n%(prog)s -s timema -i1 <alignment_directory> -i2 <branch_lengths> -i3 <control_file> -e <email> -o <output_file>')
	mod.append('%(prog)s -s table_3rate -i1 <output_directory> -o <output_file>')
	mod.append('%(prog)s -s table_2rate -i1 <output_directory> -o <output_file>')
	mod.append('%(prog)s -s table_free -i1 <output_directory> -o <output_file>')	
		
	parser = argparse.ArgumentParser(prog = 'timema2codeml.py',
                                 usage = "\n".join(mod))

	parser.add_argument('-s', action='store', dest='step_value',
	                    help='Step')
	                                                 	
	parser.add_argument('-i1', action='store', dest='input_value',
	                    help='Input 1')

	parser.add_argument('-i2', action='store', dest='input2_value',
	                    help='Input 2')

	parser.add_argument('-i3', action='store', dest='input3_value',
	                    help='Input 3')

	parser.add_argument('-e', action='store', dest='email_value',
	                    help='Email')
       	
	parser.add_argument('-o', action='store', dest='output_value',
	                    help='Output')
	                    	                    	
	parser.add_argument('--version', action='version', version='%(prog)s 1.0')

	results = parser.parse_args()
	
	# Run PAML codeml in Vital-IT cluster.
		
	if results.step_value == "timema" and results.input_value and results.input2_value and results.input3_value and results.email_value and results.output_value:
		timema(results.input_value, results.input2_value, results.input3_value, results.email_value, results.output_value)

	# Table output.

	if results.step_value == "table_3rate" and results.input_value and results.output_value:
		table_3rate(results.input_value, results.output_value)

	if results.step_value == "table_2rate" and results.input_value and results.output_value:
		table_2rate(results.input_value, results.output_value)

	if results.step_value == "table_free" and results.input_value and results.output_value:
		table_free(results.input_value, results.output_value)
													
if __name__ == "__main__":
	main(sys.argv[1:])
				
