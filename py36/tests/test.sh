#!/bin/bash
# usage: bash tool_wrapper.sh

echo "RUN bash install_conda_env.sh FIRST"
export PATH=~/miniconda2/bin:$PATH
. activate RepEnrich_py36_0
input_base=Samp
baseReference=chrM

bowtie-build ${baseReference}.fa ${baseReference}

python ../RepEnrich_setup.py ${baseReference}_repeatmasker.txt ${baseReference}.fa \
       setup_folder_${baseReference}

bowtie $baseReference -p 16 -t -m 1 -S --max ${input_base}_multimap.fastq \
       ${input_base}.fastq ${input_base}_unique.sam

samtools view -bS ${input_base}_unique.sam > ${input_base}_unique.bam
samtools sort ${input_base}_unique.bam ${input_base}_unique_sorted
mv ${input_base}_unique_sorted.bam ${input_base}_unique.bam
samtools index ${input_base}_unique.bam
rm ${input_base}_unique.sam

python ../RepEnrich.py ${baseReference}_repeatmasker.txt ${input_base} ${input_base} \
        setup_folder_${baseReference} ${input_base}_multimap.fastq ${input_base}_unique.bam --cpus 16
