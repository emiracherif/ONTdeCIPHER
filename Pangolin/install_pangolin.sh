#!/bin/bash

rm -rf pangolin

# Clone Pangolin along with it's submodules
git clone --recursive https://github.com/cov-lineages/pangolin.git


# Go into Pangolin fold --> create pangolin env --> install pangolin dependencies
cd pangolin
conda env create --name pangolin --file=environment.yml

pathconda=$(echo $CONDA_PREFIX | sed 's/envs\/ontdecipher//g')

echo $pathconda
source $pathconda'etc/profile.d/conda.sh'


# Activate Pangolin env --> install Pangolin 
conda activate pangolin

conda install -c bioconda pysam -y
pip install . 


# Activate Pangolin env --> install Pangolin
which snakemake
which pangolin
pangolin -v
pangolin -pv 
conda install -c anaconda jinja2 networkx pygments -y
conda install -c conda-forge pygraphviz -y
#####
conda deactivate
conda create -n ete3 python=3.5
conda activate ete3
conda install -c etetoolkit ete3 ete_toolchain -y

#####
conda deactivate
conda env create --name pycoqc --file=../Environments/pycoqc_env.yml



