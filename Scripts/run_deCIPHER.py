#!/usr/bin/env python

#############################################################################################
# run_snakFile.py
#
# Copyright (C) 2021 Fatou Seck Thiam , Mohammad Salma & Emira Cherif
#
# This script is used to run the pipeline steps 
#
# usage:
# python3 run_snakFile_v2.py [-h] --step STEP --params PARAMS --samples SAMPLES --threads INT
#############################################################################################
#############################################################################################

import os, sys, pwd, argparse # glob, re,
from datetime import datetime
import subprocess
from pathlib import Path
import pandas as pd
import numpy as np
import io, glob
import seaborn as sns
import matplotlib.pyplot as plt

###############################################################################
############# Some global variables ###########################################
## The list of steps
steps_list=["pycoQC","pip_core","m_r_p","all"]

## Config file exemple
config_file_exemple="""

			#barcode	sample_name
			barcode01	sample_1
			barcode02	sample_2
			barcode03	sample_3

			"Note: the columns are separated by '\t'."\n """

############# Functions's section #############################################

def getSnakeDir():
	
	""" Read the argv and extract the path to this script to use it later to 
	call the 'covid.smk' file.
    Warning: the feature is expected to have sankemake.smk file in the same 
    directory (i.e. no check
    is done).
    
    Return the path as a string """

	tool_dir="{}".format(sys.argv[0])
	directorys=tool_dir.split("/")[0:-1]
	tool_dir="/".join(directorys)
	return (tool_dir)


def printInfo():

	""" Print on screen some usefull information like : 
	the date of the run, the user who launched this script, the working 
	directory, the dirctory of this script. 
	
	Return NULL """
	
	print("\n -------------------------------",'\n',
	"This job has been started on :",datetime.today().strftime('%d-%m-%Y'),
	"at :", datetime.today().strftime('%H:%M:%S'),'\n',
	"-------------------------------",'\n',
	"Launched by:",pwd.getpwuid( os.getuid())[4],
	", with this username:",pwd.getpwuid( os.getuid())[0],'\n',
	"-------------------------------\n Working directory is:")

	subprocess.call(["pwd"])
	print("-------------------------------",'\n',
		"You call this script form:", '\n',
		getSnakeDir(),'\n')



def readConfigFile(samplesf, paramsf):

	""" Read the config file.
    
    Return sample name, params """

	mySampleDict={}
	# Open the config file
	with open(samplesf, 'r') as infile:
			for line in infile:
				if not line.startswith('#'):
					try:
						# Split file lines by 'tab'
						columns = line.split("\t")
						# Check if the line dont start by '#'
						if (len(columns)!=2):
							print("ERROR: The ",samplesf," is not compatible. There are :",str(len(columns)),
								"Columns but 2 columns are expected"'\n', 
						config_file_exemple)
							sys.exit(1)
						else:
							if columns[0] not in mySampleDict:
								mySampleDict.update({columns[0].strip(" "):columns[1].strip("\n")})
					except IndexError: 
						print("ERROR: The ",samplesf," is not compatible. Please respect this structure",'\n', 
						config_file_exemple)
						sys.exit(1)

	
	subprocess.call(["mkdir","-p","Step1_usedConfigs"])
	subprocess.call(["cp",samplesf,"Step1_usedConfigs/config_samplename.tsv"])
	if (paramsf):
		myParamDict={}
		subprocess.call(["mkdir","-p","Step1_usedConfigs"])
		subprocess.call(["cp",paramsf,"Step1_usedConfigs/config.txt"])
		with open("Step1_usedConfigs/config.txt","r") as myConfig:
			for line in myConfig:
				if not line.startswith('#'):
					try:
						line=line.replace(" ", "").replace("\"", "").replace("'", "").replace("\t", "").strip("\n")
						if line!="\n":
							#print(line)
							columns=line.split("=")
							if columns[0] not in myParamDict:
								myParamDict.update({columns[0]:columns[1].strip("\n")})
					except:
						pass
	return(mySampleDict, myParamDict)
	


####################################################	
def runPipeline(stepf,coref, mySampleDict, myParamDict): 
	
	""" Run the called script. This function need 4 parameters :
	step : the name of the step to run.
	threads : the number of core.
	mySampleDict : Dict for file name in fast5/fastq directories and samples name.
	myParamDict : Dict for parameters. 
	
	Return NULL """
	
	
	# Get snakefile path
	pathToSnake=getSnakeDir()
	###---------------------
	# Run pycoQC
	if (stepf==steps_list[0]):
		p = os.popen('pycoQC -f '+str(myParamDict['input_sequence_summary'])+' -o pycoQC_report.html') 
		print(p.read())
	
	###---------------------
	# Run snakemake
	if (stepf==steps_list[1]):
		print("Running snakemake file ...")
		subprocess.call(["mkdir","-p","DagFiles"])
		p = os.popen('snakemake -s '+pathToSnake+'/deCIPHER.smk --cores '+str(coref) +' --dag | dot -Tpdf > DagFiles/dag_start_'+str(datetime.today().strftime('at_%H-%M-%S_on_%d-%m-%Y'))+'.pdf') 
		print(p.read())

		p = os.popen('snakemake -s '+pathToSnake+'/deCIPHER.smk --cores '+str(coref) +' --use-conda') 
		print(p.read())

		p = os.popen('snakemake --report report.html -s '+pathToSnake+'/deCIPHER.smk')
		print(p.read())

		p = os.popen('snakemake -s '+pathToSnake+'/deCIPHER.smk --cores '+str(coref) +' --dag | dot -Tpdf > DagFiles/dag_end_'+str(datetime.today().strftime('at_%H-%M-%S_on_%d-%m-%Y'))+'.pdf') 
		print(p.read())


	###---------------------
	# Run mafft, raxmlHPC & pangolin
	if (stepf==steps_list[2]):
		
		print("Running step2 ...")
		subprocess.call(["mkdir","-p","Step8_consensus_fasta"])
		p = os.popen('cat Step3_artic_medaka_result/*.consensus.fasta > Step8_consensus_fasta/all_fasta.fasta') 
		print(p.read())
		
		p = os.popen('mafft --thread 4 --threadtb 5 --threadit 0 --reorder --auto Step8_consensus_fasta/all_fasta.fasta > Step8_consensus_fasta/all_alignment.fasta') 
		print(p.read())
		
		p = os.popen('raxmlHPC -m GTRGAMMA -p 12345 -s Step8_consensus_fasta/all_alignment.fasta -n '+str(myParamDict['name'])+' -f a -x 1000 -N 100') 
		print(p.read())
		
		p = os.popen('pangolin --panGUIlin Step8_consensus_fasta/all_fasta.fasta --threads '+str(coref)) 
		print(p.read())

		inputfiles=""
		for key, val in mySampleDict.items():
			inputfiles+="Step3_artic_medaka_result/"+val+"_"+key+".trimmed.rg.sorted.bam "
		
		print(inputfiles)
		p = os.popen("plotCoverage --bamfiles "+inputfiles+"-o Summary/All.pdf --smartLabels --plotWidth 25 --plotHeight 7 --outRawCounts Summary/All.txt --outCoverageMetrics Summary/All.txt --plotFileFormat pdf -p "+str(coref))
		print(p.read())


		####
		samples_dir = sorted(glob.glob("Summary/*_snake1.txt"))
		listSample=[]
		for samp in samples_dir:
			listSample.append(samp)
		print(listSample)
		
		df_all="sample,file,format,type,num_seqs,sum_len,min_len,avg_len,max_len\n"
		for file in listSample:
			with open(file, 'r') as infile:
				nameSample=file.split("_snake1")[0].split("/")[-1]
				
				df = ""
				for line in infile:
					newLine=""
					if (not line.startswith('#')):
						if (not line.startswith('--')):
							if (not line.startswith('file')):
								columns =line.split(" ")
								for col in columns:
									if col != "":
										if col=="-":
											col="before_filtring"
										if col.startswith('Step2'):
											col ="after_filtring"
										
										newLine+= col+","
								
								print(nameSample+","+newLine.strip(",")+"\n")
								df+=nameSample+","+newLine.strip(",")+"\n"
			df_all+=df
		
		###
		data = io.StringIO(df_all)
		df2 = pd.read_csv(data, sep=",")


		# set plot style: grey grid in the background:
		sns.set(style="darkgrid")

		# Set the figure size
		plt.figure(figsize=(20, 15))

		# grouped barplot
		sns.barplot(x="sample", y="num_seqs", hue="file", data=df2, ci=None);
		plt.xticks(rotation=45)
		plt.ylabel("number of sequences")
		plt.xlabel("Samples")
		plt.title("Summary")
		plt.savefig("Summary/stat.pdf")

	###---------------------
	# Run covid.smk file
	if (stepf==steps_list[3]):
		print("Running all the pipeline : pycoQC --> pip_core --> mafft, raxmlHPC and pangolin ;)")
		runPipeline("pycoQC",coref,mySampleDict,myParamDict) 
		runPipeline("pip_core",coref,mySampleDict,myParamDict)
		runPipeline("m_r_p",coref,mySampleDict,myParamDict)


			

####################################################
####################################################


def main():
	# Get some infos
	printInfo()

	##############################################################################
	
	parser = argparse.ArgumentParser(description='run_snakFile_v2.py')
	parser.add_argument('--step', help='Select step: pycoQC, pip_core, m_r_p or all. pycoQC: runs pycoQC tool. pip_core: runs deCIPHER.smk. m_r_p: runs mafft, raxmlHPC and pangolin', required=True)
	parser.add_argument('--params', help='Config file containing some parameters to run the pipeline.', required=True)
	parser.add_argument('--samples', help='Config file to associate barcodes with sample names.', required=True)
	parser.add_argument('--threads','-t', type=int , default= 4, help='The number of cores to run the pipeline')
	args = parser.parse_args()

	mySampleDict, myParamDict = readConfigFile(args.samples,args.params)
	
	if args.step in steps_list:
		runPipeline(args.step,args.threads, mySampleDict, myParamDict) 
	else:
		print("ERROR: ",args.step," is not supported. Please select one these:",steps_list)
		sys.exit(1)

	####################################################
	print(" Done ! ")


if __name__ == '__main__':
	main()
sys.exit(0) 
