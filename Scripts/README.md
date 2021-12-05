# Companion script

ONTdeCIPHER is an Oxford Nanopore Technology (ONT) amplicon-based sequencing pipeline to perform key downstream analyses on raw sequencing data from quality testing to SNPs effect to phylogenetic analysis. 
For user who would like to aligne ONTdeCIPHER Csonsensus fasta file `Step9_consensus_fasta/all_fasta.fasta` to custom fasta reference file, we provide a `companion_script.py`.
 
## To use companion_script.py :

You have to run the master script `companion_script.py` from working directory by:
1) activating ONTdeCIPHER conda environment:

```sh
conda activate ontdecipher
```

2) running the master script

```sh
python3 absolute_path_to_script_directory/companion_script.py [-h] --databases user_refrence.fasta --fastaCsonsensus Step9_consensus_fasta/all_fasta.fasta --params config.txt --threads 4 --output user_custom
```

`--databases` : a custom refrence fasta file.

`--fastaCsonsensus` : it can be Step9_consensus_fasta/all_fasta.fasta or any fasta consensus file.

`--params` : the config.txt file containing the parameters to run the script.
#### Example:

	### SARS-CoV-2 reference in the snpEff data base.
	reference_genome_snpEff ="MN908947.3"

	### primer schemes 
	primers="nCoV-2019/V3"

	### mafft config 
	maxambiguous=0.4

	### raxmlHPC config 
	number_distinct_starting_trees=100

`--threads/-t` : Maximum number of threads to use. Default: 4

`--output` : a output name. Default: user_custom

### companion_script.py output results

After running ONTdeCIPHER steps you will hava in your working directory the following files and folders.


	├── RAxML_bestTree.testComScript2
	├── RAxML_bipartitions.testComScript2
	├── RAxML_bipartitionsBranchLabels.testComScript2
	├── RAxML_bootstrap.testComScript2
	├── RAxML_info.testComScript2
	├── Step1_usedConfigs
	├── Step8_consensus_fasta
	├	├── testComScript2_all_alignment.fasta
	├	├── testComScript2_all_alignment.fasta.reduced
	├	└── testComScript2_reference.fasta
	└── Summary
		├── RAxML_bestTree_testComScript2.pdf
		├── RAxML_bestTree_testComScript2.svg
		├── RAxML_bipartitions_testComScript2.pdf
		├── RAxML_bipartitions_testComScript2.svg
		├── c_RAxML_bestTree_testComScript2.pdf
		├── c_RAxML_bestTree_testComScript2.svg
		├── c_RAxML_bipartitions_testComScript2.pdf
		└── c_RAxML_bipartitions_testComScript2.svg

Data to test ONTdeCIPHER and the results you will obtain are available here:

https://osf.io/jd2vz/?view_only=6d333ddc5a3045d297d5e3cc59e7e461
## Tips

**MAFFT**

If you are amplifying and analyzing environmental or wastewater samples with fragmented genomes, you may need to modify some mafft and pangolin options.

So you can choose the percentage of ambiguous nucleotides to tolerate by modifying `maxambiguous =` (mafft config), `max-ambig =` (pangolin config) in the `config.txt` file. 

## References
### Pre-processing and quality control
Adrien Leger, Tommaso Leonardi, February 28, 2019, pycoQC, interactive quality control for Oxford Nanopore Sequencing
(https://tleonardi.github.io/pycoQC/)

Wei Shen, Shuai Le, Yan Li, Fuquan Hu, October 5, 2016, SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation
(https://bioinf.shenwei.me/seqkit/)  

Fidel Ramírez, Friederike Dündar, Sarah Diehl, Björn A. Grüning, Thomas Manke, May 05, 2014, deepTools: a flexible platform for exploring deep-sequencing data
(https://deeptools.readthedocs.io/en/develop/content/installation.html)

Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report
(https://multiqc.info/)

### Genome reconstruction and genomic analysis
Nick Loman, Andrew Rambaut, Jannuary 22, 2020, nCoV-2019 novel coronavirus bioinformatics environment setup 
(https://artic.network/ncov-2019/ncov2019-it-setup.html)

Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences
(https://github.com/lh3/minimap2)

Sedlazeck, F.J., Rescheneder, P., Smolka, M. et al. 2018.Accurate detection of complex structural variations using single-molecule sequencing
(https://github.com/fritzsedlazeck/Sniffles)

Pablo Cingolani et al, April 01, 2012, SnpEff: A program for annotating and predicting the effects of single nucleotide polymorphisms  
(http://pcingola.github.io/SnpEff/se_introduction/)

Katoh, Rozewicki, Yamada, 2019, MAFFT online service: multiple sequence alignment, interactive sequence choice and visualization    
(https://mafft.cbrc.jp/alignment/software/)

Alexandros Stamatakis, January 21, 2014, RAxML: a tool for phylogenetic analysis and post-analysis of large phylogenies 
(https://www.metagenomics.wiki/tools/phylogenetic-tree/construction/raxml) 

Huerta-Cepas, J., Serra, F., & Bork, P. 2016. ETE 3: reconstruction, analysis, and visualization of phylogenomic data
(http://etetoolkit.org/)

Andrew Rambaut et al, July 15 2020, A dynamic nomenclature proposal for SARS-CoV-2 lineages to assist genomic epidemiology 
(https://cov-lineages.org/pangolin.html)


## Project contributors

* Fatou Seck Thiam (BCD/ISEM/IRD), Marine Combe (ISEM/IRD),Georgina Rivera-Ingraham (ISEM/IRD), Fabienne Justy (ISEM/UM), Damien Breugnot (ISEM/IRD), Mohammad Salma (IGMM), Marie-ka Tilak (ISEM), Jean-Claude Doudou (IRD Cayenne), Rodolphe E. Gozlan (ISEM/IRD) and Emira Cherif (ISEM/IRD).
* **Written by Mohammad Salma, Fatou Seck Thiam and Emira Cherif.**


## Contact 
emira.cherif@ird.fr

## License
Licencied under CeCill-C (http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html) and GPLv3.                                                  

Intellectual property belongs to IRD and authors.

 

