# HiCBin: Binning metagenomic contigs and recovering metagenome-assembled genomes using Hi-C contact maps

## Introduction
HiCBin is a new open-source metagenomic Hi-C-based binning pipeline to recover high-quality MAGs. HiCBin employs the HiCzin normalization method and the Leiden community detection algorithm, and includes the spurious contact detection into binning pipelines for the first time.

## Install and Setup
### conda
We recommend using conda to run HiCBin.

After installing Anaconda (or miniconda), Users can clone the repository with git
```
git clone --recursive https://github.com/dyxstat/HiCBin.git
```

Once complete and assuming you are in the repository folder, then create a HiCBin environment using conda:
```
# enter the HiCBin folder
cd HiCBin
# Construct environment
conda env create -f HiCBin_env.yaml
# Enter the environment
conda activate HiCBin_env
```

Normalization method in HiCBin depends on R package 'glmmTMB', which is installed in R:
```
# Enter the R
R
# download the R package 
install.packages("glmmTMB", type="source")
```
Finally, you can test the pipeline, and testing result are in test/out/hicbin.log:
```
python ./hicbin.py test test/out
```

## Initial data preparation
### 1.Preprocessing Raw reada
Adaptor sequences are removed by bbduk from the BBTools suite with parameter ‘ktrim=r k=23 mink=11 hdist=1 minlen=50 tpe tbo’ and reads are quality-trimmed using bbduk with parameters ‘trimq=10 qtrim=r ftm=5 minlen=50’. Then, the first 10 nucleotides of each read are trimmed by bbduk with parameter ‘ftl=10’.
### 2.Shotgun assembly
For the shotgun library, de novo metagenome assembly is produced by MEGAHIT with parameters ‘-min-contig-len 300 -k-min 21 -k-max 141 -k-step 12 -merge-level 20,0.95’ and contigs shorter than 1 kb are discarded.
### 3.Calculating the coverage of assembled contigs
Firstly, BBmap from the BBTools suite is applied to map the shotgun reads to the assembled contigs with parameters ‘bamscript=bs.sh; sh bs.sh’. The coverage of contigs is computed using MetaBAT2 v2.12.1 script: ‘jgi summarize bam contig depths’.
### 4.Alignment of Hi-C paired-end reads
Hi-C paired-end reads are mapped to assembled contigs using BWA-MEM with parameters ‘-5SP’. Then, samtools with parameters ‘view -F 0x904’ is applied on the resulting BAM files to remove unmapped reads (0x4) and supplementary (0x800) and secondary (0x100) alignments. 
### 5.Assign taxonomy to contigs by TAXAassign
The taxonomic assignment of contigs was resolved using NCBI’s Taxonomy and its nt database by TAXAassign(v0.4) with parameters ‘-p -c 20 -r 10 -m 98 -q 98 -t 95 -a “60,70,80,95,95,98” -f’. 

## HiCBin analysis
### Usage of the HiCBin pipeline
```
python ./hicbin.py pipeline [parameters] FASTA_file BAM_file TAX_file COV_file OUT_dir
```
### Parameters
```
-e: Case-sensitive enzyme name. Use multiple times for multiple enzymes
--min-len: Minimum acceptable contig length (default 1000)
--min-mapq: Minimum acceptable alignment quality (default 30)
--min-match: Accepted alignments must at least be N matches (default 30)
--min-signal: Minimum acceptable signal (default 2)
--thres: Maximum acceptable fraction of incorrectly identified valid contacts in spurious contact detection (default 10%)
--min-binsize: Minimum bin size used in output (default 150000)
```









