#!/usr/bin/env python

#############################################################################################
# run_ONTdeCIPHER.py
#
# Copyright (C) 2021 Mohammad Salma, Fatou Seck Thiam & Emira Cherif
#
# This script is used to run the pipeline steps
#
# usage:
# python3 run_ONTdeCIPHER.py [-h] --step STEP --params PARAMS --samples SAMPLES
#                          [--threads THREADS]
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

def checkfiles(myParamDict,mySampleDict):
	filesToCheck=["input_fastq","input_fast5"] 
	dataOK = True
	for file in filesToCheck:
		if file in myParamDict.keys():
			if not os.path.isdir(myParamDict[file]):				
				dataOK= False
				print(Fore.RED +"ERROR: ",myParamDict[file],Fore.RED +": No such directory.")
				sys.exit(1)
			else:
				for sample in mySampleDict.keys():
					if not os.path.isdir(myParamDict[file]+"/"+sample):
						dataOK= False
						print(Fore.RED +"ERROR: ",myParamDict[file]+"/"+sample,Fore.RED +": No such file.")
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



def readConfigFile(samplesf, paramsf, threads):

	""" Read the config file.

    Return sample name, params """

	getSnakeDir_path= '"'+str(getSnakeDir())+'"'
	conda_path='"'+str(getCondaPath())+'"'
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
		cmd="""echo '\nscripts_path="""+getSnakeDir_path+"""' >> Step1_usedConfigs/config.txt"""
		subprocess.Popen(cmd,shell=True)
		# cmd="""echo '\nconda_path="""+conda_path+"""' >> Step1_usedConfigs/config.txt"""
		# subprocess.Popen(cmd,shell=True)
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
						if line.startswith('usher'):
							myParamDict.update({"usher":"usher"})
					except:
						pass
		if "conda_path" not in myParamDict.keys():
			print(Fore.YELLOW +"User doesn't provide the path to etc/profile.d/conda.sh, ontdecipher will try to detect it!")
			cmd="""echo '\nconda_path="""+conda_path+"""' >> Step1_usedConfigs/config.txt"""
			subprocess.Popen(cmd,shell=True)
		else: 
			if os.path.isfile(str(myParamDict["conda_path"])+"/etc/profile.d/conda.sh"):
				pass
			else:
				print(Fore.RED +"ERROR: ",str(myParamDict["conda_path"])+"/etc/profile.d/conda.sh",Fore.RED +": No such file.")
				sys.exit(1)


	#print(myParamDict)
	return(mySampleDict, myParamDict)
####################################################
def generateStats(mySampleDict,myParamDict,coref):
	print(Fore.GREEN +"Generate stats ...")
	inputfiles=""
	for key, val in mySampleDict.items():
		inputfiles+="Step3_artic_medaka_result/"+val+"_"+key+".trimmed.rg.sorted.bam "

	print(inputfiles)
	p = os.popen("plotCoverage --bamfiles "+inputfiles+"-o Summary/All.pdf --smartLabels --plotWidth 25 --plotHeight 7 --outRawCounts Summary/All.txt --outCoverageMetrics Summary/All.txt --plotFileFormat pdf -p "+str(coref)+ " >> Logs/all_samples_plotCoverage.log 2>&1")
	print(p.read())

	samples_dir = sorted(glob.glob("Summary/*_snake_Q.txt"))
	if "input_fast5" in myParamDict.keys():
		samples_dir = sorted(glob.glob("Summary/*_snake_Q5.txt"))

	listSample=[]
	for samp in samples_dir:
		listSample.append(samp)
	
	#print(listSample)

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

							#print(nameSample+","+newLine.strip(",")+"\n")
							df+=nameSample+","+newLine.strip(",")+"\n"
		df_all+=df

	###
	data = io.StringIO(df_all)
	df2 = pd.read_csv(data, sep=",")
	df2['sample'] = df2['sample'].str.replace("_snake_Q.txt","")
	df2['sample'] = df2['sample'].str.replace("_snake_QC.txt","")

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
	plt.savefig("Summary/stat.pdf", bbox_inches='tight')



	gb = df2.groupby("file")
	listTilte=["After","Before"]
	i=0
	for x in gb.groups:

		df3= pd.DataFrame(gb.get_group(x),columns=["sample","min_len","avg_len","max_len"])
		df3['sample'] = df3['sample'].str.replace("_snake_Q.txt","")
		df3['sample'] = df3['sample'].str.replace("_snake_QC.txt","")
		print (df3)
		df3 = df3.set_index('sample')
		plt.figure(figsize=(70, 75))
		df3.plot.barh(stacked=True, title= listTilte[i]);
		plt.xticks(rotation=45, fontsize=7)
		plt.yticks(rotation=45, fontsize=7)
		plt.xlabel("number of sequences")
		plt.ylabel("Samples")
		plt.legend(borderaxespad=0, fontsize=7)
		i+=1
		plt.savefig("Summary/summary_"+str(listTilte[i-1])+".pdf", bbox_inches='tight')


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
		cmd='source '+str(myParamDict['conda_path'])+'/etc/profile.d/conda.sh && conda activate pycoqc && pycoQC -f '+str(myParamDict['input_sequence_summary'])+' -o pycoQC_report.html && conda deactivate'
		subprocess.Popen(cmd, shell=True, executable='/bin/bash')


	###---------------------
	# Run snakemake
	if (stepf==steps_list[1]):
		print("Running snakemake file ...")
		subprocess.call(["mkdir","-p","DagFiles"])

		snakeFileTarget="ontdeCIPHER.smk"
		if "input_fast5" in myParamDict.keys():
			snakeFileTarget="ontdeCIPHER_Q5.smk"

		# p = os.popen('snakemake -s '+pathToSnake+'/'+snakeFileTarget+' --cores '+str(coref) +' --dag | dot -Tpdf > DagFiles/dag_start_'+str(datetime.today().strftime('at_%H-%M-%S_on_%d-%m-%Y'))+'.pdf')
		# print(p.read())

		cmd='snakemake -s '+pathToSnake+'/'+snakeFileTarget+' --cores '+str(coref) +' --dag | dot -Tpdf > DagFiles/dag_start_'+str(datetime.today().strftime('at_%H-%M-%S_on_%d-%m-%Y'))+'.pdf'
		subprocess.Popen(cmd, shell=True, executable='/bin/bash')

		p = os.popen('snakemake -s '+pathToSnake+'/'+snakeFileTarget+' --cores '+str(coref))
		print(p.read())

		snakemake_version=str(snakemake.__version__)

		if snakemake_version == '3.13.3':
			cmd='source '+str(myParamDict['conda_path'])+'/etc/profile.d/conda.sh && conda activate pangolin && snakemake --report report_'+str(datetime.today().strftime('at_%H-%M-%S_on_%d-%m-%Y'))+'.html -s '+pathToSnake+'/'+snakeFileTarget
			subprocess.Popen(cmd, shell=True, executable='/bin/bash')
		else:	
			p = os.popen('snakemake --report report_'+str(datetime.today().strftime('at_%H-%M-%S_on_%d-%m-%Y'))+'.html -s '+pathToSnake+'/'+snakeFileTarget)
			print(p.read())

		p = os.popen('snakemake -s '+pathToSnake+'/'+snakeFileTarget+' --cores '+str(coref) +' --dag | dot -Tpdf > DagFiles/dag_end_'+str(datetime.today().strftime('at_%H-%M-%S_on_%d-%m-%Y'))+'.pdf')
		print(p.read())


	###---------------------
	# Run mafft, raxmlHPC & pangolin
	if (stepf==steps_list[2]):
		## generate plots
		generateStats(mySampleDict,myParamDict,coref)


		if "usher" in myParamDict.keys():
			usher="--usher "
		else:
			usher=""

		if "maxambiguous" in myParamDict.keys():
			maxambiguous="--maxambiguous "+str(myParamDict['maxambiguous'])+" "
		else:
			maxambiguous="--maxambiguous 0.3 "

		if "number_distinct_starting_trees" in myParamDict.keys():
			number_distinct_starting_trees="-N "+str(myParamDict['number_distinct_starting_trees'])+" "
			print(number_distinct_starting_trees)
		else:
			number_distinct_starting_trees="-N 1000 "

		if "max-ambig" in myParamDict.keys():
			max_ambig="--max-ambig "+str(myParamDict['max-ambig'])+" "
		else:
			max_ambig="--max-ambig 0.3 "


		if "min-length" in myParamDict.keys():
			min_length="--min-length "+str(myParamDict['min-length'])+" "
		else:
			min_length="--min-length 25000 "


		print(Fore.GREEN +"Generate All fasta consensus file ...")
		subprocess.call(["mkdir","-p","Step9_consensus_fasta"])
		subprocess.call(["rm","-f"," Step9_consensus_fasta/all_fasta.fasta"])
		p = os.popen('cat Step3_artic_medaka_result/*.consensus.fasta > Step9_consensus_fasta/all_fasta.fasta')
		print(p.read())
		print(Fore.GREEN +"Run mafft ...")
		p = os.popen('mafft --thread '+str(coref) +' --6merpair '+str(maxambiguous)+'--addfragments Step9_consensus_fasta/all_fasta.fasta '+str(myParamDict['scripts_path'])+'/artic-ncov2019_data/primer_schemes/'+str(myParamDict['primers'])+'/'+str(myParamDict['primers'].split('/')[0])+'.reference.fasta >  Step9_consensus_fasta/all_alignment.fasta')
		print(p.read())


		print(Fore.GREEN +"Run raxmlHPC ...")
		p = os.popen('raxmlHPC -T '+str(coref) +' -m GTRGAMMA -p 12345 -s Step9_consensus_fasta/all_alignment.fasta -n '+str(myParamDict['name'])+' -f a -x 1000 '+str(number_distinct_starting_trees))
		print(p.read())

		print(Fore.GREEN +"Run pangolin ...")
		cmd='source '+str(myParamDict['conda_path'])+'/etc/profile.d/conda.sh && conda activate pangolin && pangolin '+str(max_ambig)+str(min_length)+usher+'Step9_consensus_fasta/all_fasta.fasta --threads '+str(coref)+' && conda deactivate'
		subprocess.Popen(cmd, shell=True, executable='/bin/bash')
		#print(p.read())

		# ####
		if os.path.isfile("RAxML_bestTree."+str(myParamDict["name"])):

			print(Fore.GREEN +"Plot trees ...")
			cmd="source "+str(myParamDict["conda_path"])+"/etc/profile.d/conda.sh && conda activate ete3 && python3 "+pathToSnake+"/plot_tree.py --file RAxML_bestTree."+str(myParamDict["name"])+" --output RAxML_bestTree_"+str(myParamDict["name"])+" --format svg && python3 "+pathToSnake+"/plot_tree.py --file RAxML_bipartitions."+str(myParamDict["name"])+" --output RAxML_bipartitions_"+str(myParamDict["name"])+" --format svg && conda deactivate"
			subprocess.run(cmd, shell=True, executable='/bin/bash')

			cmd="source "+str(myParamDict["conda_path"])+"/etc/profile.d/conda.sh && conda activate ete3 && python3 "+pathToSnake+"/plot_tree.py --file RAxML_bestTree."+str(myParamDict["name"])+" --output RAxML_bestTree_"+str(myParamDict["name"])+" --format pdf && python3 "+pathToSnake+"/plot_tree.py --file RAxML_bipartitions."+str(myParamDict["name"])+" --output RAxML_bipartitions_"+str(myParamDict["name"])+" --format pdf && conda deactivate"
			subprocess.run(cmd, shell=True, executable='/bin/bash')

		####
		print("Run multiqc")
		p = os.popen("sleep 20 && multiqc . --dirs-depth 2")
		print(p.read())



	###---------------------
	# Run covid.smk file
	if (stepf==steps_list[3]):
		print("Running all the pipeline : pycoQC --> pip_core --> mafft, raxmlHPC and pangolin ;)")
		if "input_sequence_summary" in myParamDict.keys():
			runPipeline("pycoQC",coref,mySampleDict,myParamDict)
		runPipeline("pip_core",coref,mySampleDict,myParamDict)
		runPipeline("m_r_p",coref,mySampleDict,myParamDict)




####################################################
####################################################


def main():
	# Get some infos
	tracemalloc.start()

	time_start = time.perf_counter()
	printInfo()

	##############################################################################

	parser = argparse.ArgumentParser(description='run_ONTdeCIPHER.py')
	parser.add_argument('--step', help='Select step: pycoQC, pip_core, m_r_p or all. pycoQC: runs pycoQC tool. pip_core: runs ontdeCIPHER.smk. m_r_p: runs mafft, raxmlHPC and pangolin', required=True)
	parser.add_argument('--params', help='Config file containing some parameters to run the pipeline.', required=True)
	parser.add_argument('--samples', help='Config file to associate barcodes with sample names.', required=True)
	parser.add_argument('--threads','-t', type=int , default= 4, help='The number of cores to run the pipeline')
	args = parser.parse_args()

	mySampleDict, myParamDict = readConfigFile(args.samples,args.params,args.threads)

	if args.step in steps_list:
		if checkfiles(myParamDict, mySampleDict):
			runPipeline(args.step,args.threads, mySampleDict, myParamDict)

	else:
		print("ERROR: ",args.step," is not supported. Please select one these:",steps_list)
		sys.exit(1)

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
