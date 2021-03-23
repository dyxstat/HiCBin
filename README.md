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
bwa index final.contigs.fa
bwa mem -5SP final.contigs.fa hic_read1.fastq.gz hic_read2.fastq.gz > MAP.sam
samtools view -F 0x904 -bS MAP.sam > MAP_UNSORTED.bam
samtools sort MAP_UNSORTED.bam -o . > MAP_SORTED.bam
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
### Pipeline of the HiCBin 
```
python ./hicbin.py pipeline [parameters] FASTA_file BAM_file TAX_file COV_file OUTPUT_directory
```
#### Parameters
```
-e: Case-sensitive enzyme name. Use multiple times for multiple enzymes 
(Optional; If no enzyme is input, HiCzin_LC mode will be employed to do normalization)
--min-len: Minimum acceptable contig length (default 1000)
--min-mapq: Minimum acceptable alignment quality (default 30)
--min-match: Accepted alignments must be at least N matches (default 30)
--min-signal: Minimum acceptable signal (default 2)
--thres: Maximum acceptable fraction of incorrectly identified valid contacts in spurious contact detection (default 10%)
--min-binsize: Minimum bin size used in output (default 150000)
-v: Verbose output
```
#### Input File
```
FASTA_file: a fasta file of the assembled contig (e.g. final.contigs.fa)
BAM_file: a bam file of the Hi-C alignment (e.g. MAP_SORTED.bam)
TAX_file: a csv file of contigs' taxonomy assignment by TAXAassign (e.g. contig_tax.csv)
COV_file: a txt file of contigs' coverage information computed by script: ‘jgi summarize bam contig depths’ from MetaBAT2 (e.g. depth.txt)
```

#### Example
```
python ./hicbin.py pipeline -e HindIII -e NcoI -v final.contigs.fa MAP_SORTED.bam contig_tax.csv depth.txt out
```
If the restriction enzymes employed in the experiment are unspecified, use
```
python ./hicbin.py pipeline -v final.contigs.fa MAP_SORTED.bam contig_tax.csv depth.txt out
```
The results of the pipeline action are all in the 'out' directory. The draft genomic bins are in 'out/BIN' and 'hicbin.log' file contains the specific implementation information of HiCBin.

### Post-processing step of HiCBin
This is used to process the partially contaminated bins with completeness larger than 50% and contamination larger than 10%.
```
python ./hicbin.py recluster --cover -v FASTA_file Contaminated_Bins_file OUTPUT_directory
```
#### Input File
```
FASTA_file: a fasta file of the assembled contig (e.g. final.contigs.fa).
Contaminated_Bins_file: a csv file of the names of the partially contaminated bins; Bin names are arranged 
                        in columns and don't include the file formats .fa at the end of each name
```
Example of a Contaminated_Bins_file:
```
BIN0000
BIN0001
BIN0005
...
```
##### Please make sure that the OUTPUT_directory is the same as your directory in the pipeline action.

#### Example
```
python ./hicbin.py recluster --cover -v final.contigs.fa contaminated_bins.csv out
```
The generated sub bins are in 'out/SUB_BIN' and the specific implementation details of the post-processing step will be added to the 'hicbin.log' file.


## Contacts and bug reports
If you have any questions or suggestions, welcome to contact Yuxuan Du (yuxuandu@usc.edu).

## Copyright and License Information
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.






