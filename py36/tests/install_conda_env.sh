#!/bin/bash

wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh
export PATH="~/miniconda2/bin:$PATH"
source ~/.bashrc
conda --version
echo "channels:" > ~/.condarc
echo "  - bioconda" >> ~/.condarc
echo "  - r" >> ~/.condarc
echo "  - defaults" >> ~/.condarc
conda create --name RepEnrich_py36_0 python=3.6 bowtie=1.2.0 samtools=0.1.19 bedtools=2.20.1 biopython=1.69
