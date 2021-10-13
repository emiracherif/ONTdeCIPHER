#!/usr/bin/env python

#############################################################################################
# run_snakFile.py
#
# Copyright (C) 2021 Fatou Seck Thiam , Mohammad Salma & Emira Cherif
#
# This script is used to run the pipeline steps
#
# usage:
# python3 plot_tree.py [-h] --file file
#############################################################################################
#############################################################################################

import os, sys, pwd, argparse , subprocess, re
from ete3 import Tree, TreeStyle

def cleanTree(file):
	
	with open(file ,'r') as brutTree:
		cleanTreeFile="clean_tree_"+file
		cmd="rm -f "+cleanTreeFile
		subprocess.Popen(cmd, shell=True, executable='/bin/bash')
		with open(cleanTreeFile, 'a') as outcleanTreeFile:
			for line in brutTree:
				line = re.sub("_barcode..\/ARTIC\/medaka", "", line)
				outcleanTreeFile.write(line)


def plot_tree(file, output, format):
	cleanTree(file)
	cleanTreeFile="clean_tree_"+file
	t=Tree(cleanTreeFile)
	ts = TreeStyle()
	ts.show_leaf_name = True
	ts.show_branch_support = True
	ts.show_leaf_name = True
	ts.show_branch_length = True
	ts.show_branch_support = True
	outputfile="Summary/"+output+"."+format
	t.render(outputfile, w=1600, tree_style = ts)
	ts.mode = "c"
	ts.arc_start = -180 # 0 degrees = 3 o'clock
	ts.arc_span = 180
	outputfile="Summary/c_"+output+"."+format
	t.render(outputfile, w=1600, tree_style = ts)

def main():
	##############################################################################

	parser = argparse.ArgumentParser(description='plot_tree.py')
	parser.add_argument('--file', help='Select a newick file.', required=True)
	parser.add_argument('--output', type=str , default= "best_tree", help='Output file name.')
	parser.add_argument('--format', type=str , default= "pdf", help='pdf or png, default: pdf.')
	args = parser.parse_args()
	if args.format not in args:
		plot_tree(args.file,args.output,"pdf")
		plot_tree(args.file,args.output,"png")
	else:
		plot_tree(args.file,args.output,args.format)

	##############################################################################

if __name__ == '__main__':
	main()
sys.exit(0)


	