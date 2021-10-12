---
title: release-notes
summary: The release notes (copied from repo changelog).
authors:
  - Will Rowe
  - Nick Loman
date: 2020-03-30
---

**Note**: this document is copied from the repository changelog.

***

# Release notes

## v1.1.0-rc1:
* Support for read groups:
    * Support ‘pool’ read groups taken from BED file, e.g.:
        * nCoV2019_1
        * nCoV2019_2
    * These can be helpfully viewed in the Tablet viewer.
    * Do variant calling separately on each group, permits detection of incongruent mutations e.g. from contamination
* Update to nanopolish 0.12.5 fo more informative VCF output
* Change nanopore filtering to support new VCF output (e.g. strand-bias)
* Support indels with bcftools consensus
* Longshot added to Medaka to permit filter VCF on depth
* New workflow model - all commands run "one per sample", script nanopolish index
* Add 'artic gupplyplex' for support for demultiplexed Guppy output to replace 'artic gather'
* Fix for align_trim to remove erroneous CIGAR strains that can cause medaka to fail
* Support for amplicons v3 format file (Will Rowe)
* Test suite (Will Rowe)
* Documentation (Will Rowe)
* Remove old deprecated code (Will Rowe)

## v0.13:
* ARTIC fork

## v0.12: 
* Update to nanopolish v0.8.4 for Albacore 2.0+ support:
    * run nanopolish -d /path/to/fast5 input.fasta before zibra.py minion

## v0.11-dev:
* Add stringent pipeline

## v0.11 (25th July 2017):
* Revert back to Ryan Wick's porechop repo as --untrimmed now implemented
* New script: align_trim_fasta.py and reconstitute.py as basis for more principled demultiplexing (align first, then trim to known amplicon end points +/- 40 bases) to help prevent chimeric reads resulting in incorrect barcode identification

## v0.10 (15th July 2017):
* Better handling of confident deletions (previously ignored, now N-masked)
* Update to nanopolish HEAD
* Fixes to support latest nanopolish VCF format changes
* Ability to choose higher values of max-haplotypes to support more divergent references

## v0.9:
* New command line interface zibra.py to ease extraction and demultiplexing with Albacore and multiple groups

## v0.8:
* Track latest Porechop with better adaptor trimming

## v0.7:
* Add Illumina wrapper scripts

## v0.6:
* Update poretools

## v0.5:
* Update nanopolish to latest HEAD - to fix Albacore 1.1 support (7de633d01cc35a58e5537af5dd1024ae0040d15c)
* Add Lassa virus scheme

## v0.4:
* Update nanopolish to latest HEAD - aaaf90beb4efe37243c651e62fcfd11374fe2176 for Albacore 1 support
* Initial support for Yellow Fever (500 and 1000 bp schemes)
    * Move to poretools for extraction with named basecaller

## v0.3:
* add Porechop

## v0.2:
* update nanopolish to 0.6 - adds default support for various nanopore models
* update samtools to 1.4
