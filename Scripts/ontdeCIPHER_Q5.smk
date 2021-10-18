#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
# deCIPHER.smk
#
# Copyright (C) 2021 Fatou Seck Thiam , Mohammad Salma & Emira Cherif
#
###############################################################################
############################## Tools ##########################################
import glob, os, re, subprocess
###############################################################################
dataPath=os.getcwd()

############################# Variables #######################################
# Get config info

mySampleDict={}
with open("Step1_usedConfigs/config_samplename.tsv","r") as Samplename:
	for line in Samplename:
		if not line.startswith('#'):
			try:
				columns=line.split("\t")
				if columns[0] not in mySampleDict:
					mySampleDict.update({columns[0].strip(" "):columns[1].strip("\n")})
			except:
				pass
			

# Get parameters
myParamDict={}
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
			

sequencing_summary=""	
if "input_sequence_summary" in myParamDict.keys():
	sequencing_summary="--sequencing-summary "+str(myParamDict['input_sequence_summary'])


sniffles="sniffles"
if "sniffles" in myParamDict.keys():
	sniffles=str(myParamDict['sniffles'])


samples_dir = glob.glob(str(myParamDict['input_fastq'])+"/barcode*")


samples=[]
for smpl in samples_dir:
	smpl=smpl.split("/")[-1]
	smpl=mySampleDict[smpl]+"_"+smpl
	samples.append(smpl)


#Define the target files
FINAL_FILES = ['Summary/'+sample+'_snake_Q5.txt' for sample in samples]




foldersList= ["Summary",
				"Step2_artic_guppyplex_filter", "Step3_artic_medaka_result","Step4_artic_nanopolish_result",
				"Step5_snpEff_result", "Step5_snpEff_result/BedFiles", "Step5_snpEff_result/HtmlFiles",
				"Step6_plotCoverage_result", "Step7_bamCoverage_result","Step8_sniffles_result"
				]
#Creat output directories
for folder in foldersList:
	subprocess.call(["mkdir", "-p",folder])
###############################################################################
########################### Rules #############################################
###############################################################################

rule all:
	input: FINAL_FILES
###############################################################################

rule artic_guppyplex:
	input:
		fastq_Dir= lambda wildcards: expand(str(myParamDict['input_fastq'])+"/{name}/", name=wildcards.sample.split("_")[-1])

	output:
		out1="Step2_artic_guppyplex_filter/{sample}.fastq"
	

	shell:
		"artic guppyplex --skip-quality-check --min-length "+myParamDict['min']+" --max-length "+myParamDict['max']+" --directory {input.fastq_Dir} --output {output.out1} "


rule artic_minion_medaka:
	input:
		fast5_Dir = lambda wildcards: expand(str(myParamDict['input_fast5'])+"/{name}/", name=wildcards.sample.split("_")[-1]),
		artic_guppyplex_fastq="Step2_artic_guppyplex_filter/{sample}.fastq"

	output:
		out1 = "Step3_artic_medaka_result/{sample}.minion.log.txt",
		out2 = "Step3_artic_medaka_result/{sample}.pass.vcf.gz",
		out3 = "Step3_artic_medaka_result/{sample}.trimmed.rg.sorted.bam",
		out4 = "Step3_artic_medaka_result/{sample}.sorted.bam",
		out5 = "Step3_artic_medaka_result/{sample}.sorted.bam.bai"
	
	threads: round(float(myParamDict['cpu'])/2)
	shell:
		"cd Step3_artic_medaka_result/ && artic minion --medaka --medaka-model "+str(myParamDict['medaka_model'])+" "
		"--normalise 200 --threads {threads} --strict --scheme-directory "+str(myParamDict['scripts_path'])+"/artic-ncov2019_data/primer_schemes "+sequencing_summary+
		"--read-file ../{input.artic_guppyplex_fastq} --fast5-directory {input.fast5_Dir} "+str(myParamDict['primers'])+" {wildcards.sample}"


rule sniffels:
	input:
		bam="Step3_artic_medaka_result/{sample}.sorted.bam",
		bai="Step3_artic_medaka_result/{sample}.sorted.bam.bai"
	output:
		out="Step8_sniffles_result/{sample}.vcf"
	
	threads: round(float(myParamDict['cpu'])/2)
	run:
		shell(sniffles+" -m {input.bam} -t {threads} -v {output.out}")

rule artic_minion_nanopolish:
	input:
		fast5_Dir = lambda wildcards: expand(str(myParamDict['input_fast5'])+"/{name}/", name=wildcards.sample.split("_")[-1]),
		artic_guppyplex_fastq="Step2_artic_guppyplex_filter/{sample}.fastq"

	output:
		out1="Step4_artic_nanopolish_result/{sample}.minion.log.txt"
	
	threads: round(float(myParamDict['cpu'])/2)
	shell:
		"cd Step4_artic_nanopolish_result/ && artic minion --normalise 200 --threads {threads} --strict --scheme-directory "+str(myParamDict['scripts_path'])
		"/artic-ncov2019_data/primer_schemes --read-file ../{input.artic_guppyplex_fastq}  --fast5-directory {input.fast5_Dir} "+sequencing_summary+
		" "+str(myParamDict['primers'])+" {wildcards.sample}"


rule snpEff:
	input:
		vcf = "Step3_artic_medaka_result/{sample}.pass.vcf.gz"

	output:
		bed="Step5_snpEff_result/BedFiles/{sample}.ann.bed",
		html="Step5_snpEff_result/HtmlFiles/{sample}.html"

	shell:
		"snpEff ann -v "+str(myParamDict['reference_genome_snpEff'])+" {input.vcf} -s {output.html} -o bed > {output.bed}"


rule plotCoverage:
	input:
		"Step3_artic_medaka_result/{sample}.trimmed.rg.sorted.bam"

	output:
		pdf="Step6_plotCoverage_result/{sample}.pdf",
		txt="Step6_plotCoverage_result/{sample}.txt",
		pdf2="Step6_plotCoverage_result/{sample}_2.pdf",
		txt2="Step6_plotCoverage_result/{sample}_2.txt"

	threads: round(float(myParamDict['cpu'])/2)
	run:
		shell("plotCoverage -b {input} -o {output.pdf} --smartLabels -T {wildcards.sample} --outRawCounts "
		"{output.txt} --outCoverageMetrics {output.txt} --plotFileFormat pdf -p {threads} --plotWidth 25 --plotHeight 7"),  #$bar_sequencing_depth  not used
		shell("plotCoverage -b {input} -o {output.pdf2} --labels {wildcards.sample} -T {wildcards.sample} --outRawCounts "
		"{output.txt2} --outCoverageMetrics {output.txt2} --plotFileFormat pdf -p {threads} --plotWidth 25 --plotHeight 7")


rule bamCoverage:
	input:
		"Step3_artic_medaka_result/{sample}.trimmed.rg.sorted.bam"

	output:
		bedGraph="Step7_bamCoverage_result/{sample}.bedgraph",
		bigwig="Step7_bamCoverage_result/{sample}.bigwig"

	threads: round(float(myParamDict['cpu'])/2)
	run:
		shell("bamCoverage -b {input} -o {output.bedGraph} -of 'bedgraph' -p {threads} --effectiveGenomeSize 29903 --normalizeUsing RPGC "),
		shell("bamCoverage -b {input} -o {output.bigwig} -of 'bigwig' --binSize 5 -p {threads} --effectiveGenomeSize 29903 --normalizeUsing RPGC ")

################################################################################
################################################################################
rule finalizing:
	input:
		fastq_Dir= lambda wildcards: expand(str(myParamDict['input_fastq'])+"/{name}/", name=wildcards.sample.split("_")[-1]),
		input2= "Step2_artic_guppyplex_filter/{sample}.fastq",
		input3= "Step3_artic_medaka_result/{sample}.minion.log.txt",
		input4= "Step4_artic_nanopolish_result/{sample}.minion.log.txt",
		input5= "Step5_snpEff_result/BedFiles/{sample}.ann.bed",
		input6= "Step6_plotCoverage_result/{sample}.pdf",
		bedGraph="Step7_bamCoverage_result/{sample}.bedgraph",
		bigwig="Step7_bamCoverage_result/{sample}.bigwig",
		outsniffles="Step8_sniffles_result/{sample}.vcf"
	output:
		"Summary/{sample}_snake_Q5.txt"
	message: "Finalizing and removing tmp file"
	run:
		shell("touch {output}"),
		shell("echo '#The number of reads by barcode before filtring: ' >> {output}"),
		shell("cat {input.fastq_Dir}/*.fastq* | seqkit stats -T | csvtk pretty -t >> {output}"),
		shell("echo '#-----------------------------------------------' >> {output}"),
		shell("echo '#The number of reads by barcode after filtring: ' >> {output}"),
		shell("seqkit stats {input.input2} -T | csvtk pretty -t >> {output}"),
		shell("echo '#Processing {wildcards.sample} is done!' >> {output}")
###############################################################################
###############################################################################
