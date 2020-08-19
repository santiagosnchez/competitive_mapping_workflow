#!/bin/bash

# ONLY internal
# merge and rename fastq files
# ls --color=auto | cut -d'.' -f5 | sed 's/_R1//' | sort -u > samples.txt
# cat samples.txt | parallel -j3 'cat *{}* > {}.fastq.gz' &
# rm HI.1217*
if [ -f samples.txt ]; then

    # trim reads
    cat samples.txt | parallel --plus -j3 trimmomatic SE -threads 3 {}.fastq.gz {}.trimmed.fastq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 & 
    cat samples.txt | parallel --plus mv {}.fastq.gz {}.raw.fastq.gz
    mkdir raw trimmed
    mv *trimmed.fastq.gz trimmed
    mv *raw.fastq.gz raw
else
    echo "a file named \"samples.txt\" with base sample names is required"
    echo "for example, your files look like this:"
    echo -e "9F1.fastq.gz\n9F2.fastq.gz\n9F3.fastq.gz\nBF1.fastq.gz\n...\nHM3.fastq.gz"
    echo "then, your sample.txt should look like:"
    echo -e "9F11\n9F2\n9F3\n...\nHM3"
    echo "C. briggsae files: B"
    echo "C. nigoni files: 9"
    echo "F1 hybrids: H"
    echo "males: M"
    echo "females: F"
    echo "replicate: [1-3]"
fi

