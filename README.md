# HiCBin: Binning metagenomic contigs and recovering metagenome-assembled genomes using Hi-C contact maps

## Introduction
HiCBin is a new open-source metagenomic Hi-C-based binning pipeline to recover high-quality MAGs. HiCBin employs the HiCzin normalization method and the Leiden community detection algorithm, and includes the spurious contact detection into binning pipelines for the first time.

## Install and Setup
### conda
We recommend using conda to run HiCBin.

After installing Anaconda (or miniconda), Users can clone the repository with git:
```
git clone https://github.com/dyxstat/HiCBin.git
```

Once complete, you can enter the repository folder and then create a HiCBin environment using conda:
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
### 1.Preprocess Raw reads
Adaptor sequences are removed by bbduk from the BBTools suite with parameter ‘ktrim=r k=23 mink=11 hdist=1 minlen=50 tpe tbo’ and reads are quality-trimmed using bbduk with parameters ‘trimq=10 qtrim=r ftm=5 minlen=50’. Then, the first 10 nucleotides of each read are trimmed by bbduk with parameter ‘ftl=10’.
### 2.Shotgun assembly
For the shotgun library, de novo metagenome assembly is produced by MEGAHIT (v1.2.9).
```
megahit -1 SG1.fastq.gz -2 SG2.fastq.gz -o ASSEMBLY --min-contig-len 1000 --k-min 21 --k-max 141 --k-step 12 --merge-level 20,0.95
```
### 3.Align Hi-C paired-end reads to assembled contigs
Hi-C paired-end reads are mapped to assembled contigs using BWA-MEM with parameters ‘-5SP’. Then, samtools with parameters ‘view -F 0x904’ is applied on the resulting BAM files to remove unmapped reads (0x4) and supplementary (0x800) and secondary (0x100) alignments. 
```
bwa mem -5SP final.contigs.fa hic_read1.fastq.gz hic_read2.fastq.gz| \
    samtools view -F 0x904 -bS - | \
    samtools sort -o hic_map.bam -
```
### 4.Assign taxonomy to contigs by TAXAassign
The taxonomic assignment of contigs was resolved using NCBI’s Taxonomy and its nt database by TAXAassign(v0.4) with parameters ‘-p -c 20 -r 10 -m 98 -q 98 -t 95 -a “60,70,80,95,95,98” -f’. 
### 5.Calculate the coverage of assembled contigs
Firstly, BBmap from the BBTools suite is applied to map the shotgun reads to the assembled contigs with parameters ‘bamscript=bs.sh; sh bs.sh’. The coverage of contigs is computed using script: ‘jgi summarize bam contig depths’ from MetaBAT2 v2.12.1.
```
bbmap.sh in1=SG1.fastq.gz in2=SG2.fastq.gz ref=final.contigs.fa out=SG_map.sam bamscript=bs.sh; sh bs.sh
jgi_summarize_bam_contig_depths --outputDepth coverage.txt --pairedContigs pair.txt SG_map_sorted.bam
```

## HiCBin analysis
### Usage of the HiCBin pipeline
```
python ./hicbin.py pipeline [parameters] FASTA_file BAM_file TAX_file COV_file OUTPUT_directory
```
### Parameters
```
-e: Case-sensitive enzyme name. Use multiple times for multiple enzymes
--min-len: Minimum acceptable contig length (default 1000)
--min-mapq: Minimum acceptable alignment quality (default 30)
--min-match: Accepted alignments must be at least N matches (default 30)
--min-signal: Minimum acceptable signal (default 2)
--thres: Maximum acceptable fraction of incorrectly identified valid contacts in spurious contact detection (default 10%)
--min-binsize: Minimum bin size used in output (default 150000)
```









