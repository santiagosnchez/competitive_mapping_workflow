mkdir STAR 
mkdir STAR/no_pcr_duplicates STAR/marked_alignments STAR/star_logs
args="--runMode alignReads --runThreadN 4 --readFilesCommand zcat --alignIntronMin 40 --alignIntronMax 15600 --outSAMtype BAM SortedByCoordinate"
export args
cat samples.txt | grep -E "^B|H" | parallel -j4 "star $args --genomeDir references/genomes/briggsae --readFilesIn trimmed/{}.trimmed.fastq.gz --outFileNamePrefix STAR/{}_Cbr"
ls STAR | grep ".bam" | parallel --plus mv STAR/{} STAR/{/Aligned.sortedByCoord.out/_Cbr}
cat samples.txt | grep -E "^9|H" | parallel -j4 "star $args --genomeDir references/genomes/nigoni --readFilesIn trimmed/{}.trimmed.fastq.gz --outFileNamePrefix STAR/{}_Cni"
ls STAR | grep ".bam" | parallel --plus mv STAR/{} STAR/{/Aligned.sortedByCoord.out/_Cni}
mv *.bam BAM
# mark PCR duplicates
#ls STAR | grep ".bam" | parallel -j6 --plus  "picard MarkDuplicates -I STAR/{} -O STAR/marked_alignments/{..}.markdup.bam -M STAR/star_logs/{..}.duplicates.txt"
#rm STAR/*.bam
# remove duplicate reads
#ls STAR/marked_alignments | parallel -j6 --plus "samtools view -F0x400 -@3 -b STAR/marked_alignments/{} > STAR/no_pcr_duplicates/{...}.dedup.bam"
ls STAR/BAM/* | parallel -j6 "samtools index {}" 

