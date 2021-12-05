#!/usr/bin/env python

#############################################################################################
# run_snakFile.py
#
# Copyright (C) 2021 Mohammad Salma, Fatou Seck Thiam & Emira Cherif
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
import time
import resource
import tracemalloc
import colorama
from colorama import Fore, Back, Style
import snakemake


colorama.init(autoreset=True)
###############################################################################
############# Some global variables ###########################################

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

def checkfiles(myParamDict):
	filesToCheck=["input_fastq","input_fast5"] 
	dataOK = True
	for file in filesToCheck:
		if file in myParamDict.keys():
			if not os.path.isdir(myParamDict[file]):				
				dataOK= False
				print(Fore.RED +"ERROR: ",myParamDict[file],Fore.RED +": No such directory.")
				sys.exit(1)
	
	if "input_sequence_summary" in myParamDict.keys():
		if not os.path.isfile(myParamDict["input_sequence_summary"]):
			dataOK= False
			print(Fore.RED +"ERROR: ",myParamDict["input_sequence_summary"],Fore.RED +": No such file.")
			sys.exit(1)

	return dataOK
		






def getCondaPath():
	commandConda="echo $CONDA_PREFIX"
	cmd=subprocess.Popen(commandConda, shell=True, stdout=subprocess.PIPE)
	cmd_data=""
	for line in cmd.stdout:
		line= line.decode(encoding="utf-8", errors="ignore")
		cmd_data+=line
	print(cmd_data)
	return cmd_data.split("/envs/ontdecipher")[0]

def printInfo():

	""" Print on screen some usefull information like :
	the date of the run, the user who launched this script, the working
	directory, the dirctory of this script.

	Return NULL """

	print(Fore.GREEN +"\n -------------------------------",'\n',
	Fore.GREEN +"This job has been started on :",datetime.today().strftime('%d-%m-%Y'),
	Fore.GREEN +"at :", datetime.today().strftime('%H:%M:%S'),'\n',
	Fore.GREEN +"-------------------------------",'\n',
	Fore.GREEN +"Launched by:",pwd.getpwuid( os.getuid())[4],
	Fore.GREEN +", with this username:",pwd.getpwuid( os.getuid())[0],'\n',
	Fore.GREEN +"-------------------------------\n Working directory is:")

	subprocess.call(["pwd"])
	print(Fore.YELLOW +"-------------------------------",'\n',
		Fore.GREEN +"You call this script form:", '\n',
		getSnakeDir(),'\n')



def readConfigFile( paramsf, threads):

	""" Read the config file.

    Return sample name, params """

	getSnakeDir_path= '"'+str(getSnakeDir())+'"'
	conda_path='"'+str(getCondaPath())+'"'

	subprocess.call(["mkdir","-p","Step1_usedConfigs"])
	if (paramsf):
		myParamDict={}
		subprocess.call(["mkdir","-p","Step1_usedConfigs"])
		subprocess.call(["cp",paramsf,"Step1_usedConfigs/config.txt"])
		cmd="""echo '\nscripts_path="""+getSnakeDir_path+"""' >> Step1_usedConfigs/config.txt"""
		subprocess.Popen(cmd,shell=True)
		cmd="""echo '\nconda_path="""+conda_path+"""' >> Step1_usedConfigs/config.txt"""
		subprocess.Popen(cmd,shell=True)
		cmd="""echo '\ncpu="""+str(threads)+"""' >> Step1_usedConfigs/config.txt"""
		subprocess.Popen(cmd,shell=True)

		subprocess.call(["sleep","3"])

		with open("Step1_usedConfigs/config.txt","r") as myConfig:
			for line in myConfig:
				if not line.startswith('#'):
					try:
						line=line.replace(" ", "").replace("\"", "").replace("'", "").replace("\t", "").strip("\n")
						if line!="\n":
							print(line)
							columns=line.split("=")
							if columns[0] not in myParamDict:
								myParamDict.update({columns[0]:columns[1].strip("\n")})
					except:
						pass
	#print(myParamDict)
	return(myParamDict)

def runCustomStep(databases, fastaCsonsensus, coref, myParamDict, output):

	""" Run the called script. This function need 4 parameters :
	step : the name of the step to run.
	threads : the number of core.
	mySampleDict : Dict for file name in fast5/fastq directories and samples name.
	myParamDict : Dict for parameters.

	Return NULL """


	# Get snakefile path
	pathToSnake=getSnakeDir()
	###---------------------

	if "maxambiguous" in myParamDict.keys():
			maxambiguous="--maxambiguous "+str(myParamDict['maxambiguous'])+" "
	else:
		maxambiguous="--maxambiguous 0.3 "

	if "number_distinct_starting_trees" in myParamDict.keys():
		number_distinct_starting_trees="-N "+str(myParamDict['number_distinct_starting_trees'])+" "
		print(number_distinct_starting_trees)
	else:
		number_distinct_starting_trees="-N 100 "


	print(Fore.GREEN +"Generate All fasta consensus file ...")
	subprocess.call(["mkdir","-p","Step9_consensus_fasta"])

	print(Fore.GREEN +"Run mafft ...")
	## mafft --add sequences.fasta --reorder ~/Documents/IE-2018/ToolDevelop/deCIPHER/ONTdeCIPHER/Scripts/artic-ncov2019_data/primer_schemes/nCoV-2019/V3/nCoV-2019.reference.fasta > output.fasta
	
	p = os.popen('mafft --add '+databases+' --thread '+str(coref) +' --reorder '+str(myParamDict['scripts_path'])+'/artic-ncov2019_data/primer_schemes/'+str(myParamDict['primers'])+'/'+str(myParamDict['primers'].split('/')[0])+'.reference.fasta >  Step9_consensus_fasta/'+str(output)+'_reference.fasta  2>> Logs/'+str(output)+'_mafft.log')
	print(p.read())

	p = os.popen('mafft --thread '+str(coref) +' --6merpair '+str(maxambiguous)+' --addfragments '+fastaCsonsensus+' Step9_consensus_fasta/'+str(output)+'_reference.fasta >  Step9_consensus_fasta/'+str(output)+'_all_alignment.fasta 2>> Logs/'+str(output)+'_mafft.log')
	print(p.read())

	print(Fore.GREEN +"Run raxmlHPC ...")
	p = os.popen('raxmlHPC -T '+str(coref) +' -m GTRGAMMA -p 12345 -s Step9_consensus_fasta/'+str(output)+'_all_alignment.fasta -n '+str(output)+' -f a -x 1000 '+str(number_distinct_starting_trees)+' >> Logs/'+str(output)+'_raxml.log 2>&1')
	print(p.read())

	# ####
	if os.path.isfile("RAxML_bestTree."+str(output)):

		print(Fore.GREEN +"Plot trees ...")
		cmd="source "+str(myParamDict["conda_path"])+"/etc/profile.d/conda.sh && conda activate ete3 && python3 "+pathToSnake+"/plot_tree.py --file RAxML_bestTree."+str(output)+" --output RAxML_bestTree_"+str(output)+" --format svg && python3 "+pathToSnake+"/plot_tree.py --file RAxML_bipartitions."+str(output)+" --output RAxML_bipartitions_"+str(output)+" --format svg && conda deactivate"
		subprocess.run(cmd, shell=True, executable='/bin/bash')

		cmd="source "+str(myParamDict["conda_path"])+"/etc/profile.d/conda.sh && conda activate ete3 && python3 "+pathToSnake+"/plot_tree.py --file RAxML_bestTree."+str(output)+" --output RAxML_bestTree_"+str(output)+" --format pdf && python3 "+pathToSnake+"/plot_tree.py --file RAxML_bipartitions."+str(output)+" --output RAxML_bipartitions_"+str(output)+" --format pdf && conda deactivate"
		subprocess.run(cmd, shell=True, executable='/bin/bash')


####################################################
####################################################


def main():
	# Get some infos
	tracemalloc.start()

	time_start = time.perf_counter()
	printInfo()

	##############################################################################

	parser = argparse.ArgumentParser(description='run_ONTdeCIPHER.py')
	parser.add_argument('--databases', help='A user databases fasta file.', required=True)
	parser.add_argument('--fastaCsonsensus', help='Step9_consensus_fasta/all_fasta.fasta or any fasta consensus file.', required=True)
	parser.add_argument('--params', help='Config file containing some parameters to run the pipeline.', required=True)
	parser.add_argument('--threads','-t', type=int , default= 4, help='The number of cores to run the pipeline')
	parser.add_argument('--output', type=str , default= 'user_custom', help='The number of cores to run the pipeline')
	args = parser.parse_args()

	myParamDict = readConfigFile(args.params,args.threads)

	
	if checkfiles(myParamDict):
		runCustomStep(args.databases, args.fastaCsonsensus, args.threads, myParamDict,args.output)

	####################################################
	
	print("###################")
	print(Fore.GREEN +"Completed at:",str(datetime.today().strftime('%H:%M:%S on %d-%m-%Y')))
	time_elapsed = (time.perf_counter() - time_start)
	dur=time.strftime('%H hour %M min %S sec', time.gmtime(time_elapsed))
	print(Fore.BLUE +"Duration: {} ".format(dur))

	# current, peak = tracemalloc.get_traced_memory()
	# print(Fore.BLUE + f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")

	tracemalloc.stop()

	print(" Done ! ")
	


if __name__ == '__main__':
	main()
sys.exit(0)
