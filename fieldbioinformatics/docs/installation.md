---
title: installation
summary: The installation guide.
authors:
  - Will Rowe
  - Nick Loman
date: 2020-03-30
---

# Installation

As of [release 1.1.0](https://github.com/artic-network/fieldbioinformatics/releases/tag/1.1.0), it is probably easiest to install the `artic pipeline` using conda. Alternatively, you can install the pipeline itself via source/pip but you will have to also satisfy the pipeline dependencies.

## Via conda

```sh
conda install -c bioconda artic
```

## Via source

### 1. installing the pipeline

```sh
git clone https://github.com/artic-network/fieldbioinformatics
cd fieldbioinformatics
python setup.py install
```

### 2. installing dependencies

The `artic pipeline` has several [software dependencies](https://github.com/artic-network/fieldbioinformatics/blob/master/environment.yml). You can solve these dependencies using the minimal conda environment we have provided:

```sh
conda env create -f environment.yml
conda activate artic
```

### 3. testing the pipeline

First check the pipeline can be called:

```
artic -v
```

To check that you have all the required dependencies, you can try the pipeline tests with both workflows:

```
./test-runner.sh nanopolish
./test-runner.sh medaka
```

For further tests, such as the variant validation tests, see [here](http://artic.readthedocs.io/en/latest/tests?badge=latest).