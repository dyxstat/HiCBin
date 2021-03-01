# HiCBin: Binning metagenomic contigs and recovering metagenome-assembled genomes using Hi-C contact maps

## Introduction
HiCBin is a new open-source metagenomic Hi-C-based binning pipeline to recover high-quality MAGs. HiCBin employs the HiCzin normalization method and the Leiden community detection algorithm, and includes the spurious contact detection into binning pipelines for the first time.

## Setup
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

