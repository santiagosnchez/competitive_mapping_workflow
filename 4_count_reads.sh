mkdir counts
# run featureCounts
featureCounts -t transcript -a references/annotations/briggsae.chr.transcript.ortho.gff -o counts/Cbr_fc.counts.txt STAR/BAM/B*_Cbr.bam &> /dev/null &
featureCounts -t transcript -a references/annotations/nigoni.chr.transcript.ortho.gff -o counts/Cni_fc.counts.txt STAR/BAM/9*_Cni.bam &> /dev/null &
# merge files and fix header
paste Cbr_fc.counts.txt Cni_fc.counts.txt | tail -n +2 | head -1 | cut -f1,7-12,19-24 | sed -E 's/STAR\/BAM\/|\.bam//g' > head.txt
paste Cbr_fc.counts.txt Cni_fc.counts.txt | tail -n +3 | cut -f7-12,19-24 | paste ../references/annotations/gene_names.txt - > species_counts.txt
cat head.txt species_counts.txt > tmp && mv tmp species_counts.txt

paste counts/Cbr_fc.counts.txt counts/Cni_fc.counts.txt | \
       tail -n +2 | paste ../reference	cut -f1,7-12,19-24 | tail -n +2 | sed -E 's/alignments\/|hf1_split\/|\.sorted\.bam//g' > real_counts.txt
cat real_counts.txt | cut -f1,5-7,11-13 > real_counts_HF1.txt
cat real_counts.txt | cut -f1,2-4,8-10 > real_counts_species.txt

# run CompMap
cd counts
ls ../STAR/BAM | grep "^H.*.bam$" | head -12 | sed -E 's/_Cbr.*|_Cni.*//' | paste - - | cut -f 1 | parallel \
"CompMap -1 ../STAR/BAM/{}_Cbr.bam -2 ../STAR/BAM/{}_Cni.bam -b1 ../references/annotations/briggsae.chr.transcript.ortho.bed -b2 ../references/annotations/nigoni.chr.transcript.ortho.bed --base {} -s1 Cbr -s2 Cni --star"
paste H*_counts.txt | cut -f1-3,5-6,8-9,11-12,14-15,17-18 > ase_counts.txt
cd ..

