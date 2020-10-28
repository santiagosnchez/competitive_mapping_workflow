#!/bin/bash
args="--readFilesCommand zcat --runThreadN 4 --runMode genomeGenerate --sjdbGTFfeatureExon CDS --genomeSAindexNbases 12"
export args
echo -e "briggsae\nnigoni" | parallel star $args --genomeDir references/genomes/{} --genomeFastaFiles references/genomes/{}.chr.fa.gz --sjdbGTFfile references/annotations/{}.chr.gff

