#!/bin/bash

# manage data
mkdir sim_read_data
ls species_briggsae_simulated_reads/ | grep fasta | head -3 | parallel --keep cp species_briggsae_simulated_reads/{} sim_read_data/Cbr_{#}.fasta
ls species_nigoni_simulated_reads/ | grep fasta | tail -3 | parallel --keep cp species_nigoni_simulated_reads/{} sim_read_data/Cni_{#}.fasta
ls allele_briggsae_simulated_reads | grep fasta | pr -2 -t -s | parallel --colsep=" " --keep "cat allele_briggsae_simulated_reads/{1} allele_nigoni_simulated_reads/{2} > sim_read_data/HF1_{#}.fasta"

# make bwa indeces
bwa index briggsae.cds.small.fa
bwa index nigoni.cds.small.fa

# align
mkdir alignments
ls sim_read_data | parallel --keep --plus \
"bwa mem briggsae.cds.small.fa sim_read_data/{} | samtools view -b > alignments/{..}_Cbr.bam"
ls sim_read_data | parallel --keep --plus \
"bwa mem nigoni.cds.small.fa sim_read_data/{} | samtools view -b > alignments/{..}_Cni.bam"
# Hybrids real allele-specific alignments
ls sim_read_data/hf1_split/ | grep "Cbr" | parallel --keep --plus \
"bwa mem briggsae.cds.small.fa sim_read_data/hf1_split/{} | samtools view -b > alignments/hf1_split/{..}.bam"
ls sim_read_data/hf1_split/ | grep "Cni" | parallel --keep --plus \
"bwa mem nigoni.cds.small.fa sim_read_data/hf1_split/{} | samtools view -b > alignments/hf1_split/{..}.bam"

ls sim_read_data/hf1_split/ | parallel --plus samtools sort -o sim_read_data/hf1_split/{..}.sorted.bam sim_read_data/hf1_split/{}


# get read names
ls sim_read_data/ | parallel --plus -j1 'cat sim_read_data/{} | grep "^>" | awk "{ x = substr(\$0,2); print x}" > readnames/{..}.txt'

# sort and index bam files
ls alignments | parallel --plus samtools sort -o alignments/{..}.sorted.bam alignments/{}
ls alignments | grep "sorted" | parallel --plus samtools index alignments/{}
ls alignments/hf1_split/ | parallel --plus samtools sort -o alignments/hf1_split/{..}.sorted.bam alignments/hf1_split/{}
ls alignments/hf1_split/ | grep "sorted" | parallel --plus samtools index alignments/hf1_split/{}

# generate GFF file
cat briggsae.cds.small.fa | awk ' \
BEGIN{ h=""; s="" } \
{ if ($0 ~ /^>/){ if (length(s) > 0){ print h"\t.\ttranscript\t1\t"(length(s))"\t.\t+\t.\t"h }; \
h = substr($0,2); s="" } else { s = s$0 }} \
END{ print h"\t.\ttranscript\t1\t"(length(s))"\t.\t+\t.\t"h }' > briggsae.cds.small.gff
cat nigoni.cds.small.fa | awk '
BEGIN{ h=""; s="" } \
{ if ($0 ~ /^>/){ if (length(s) > 0){ print h"\t.\ttranscript\t1\t"(length(s))"\t.\t+\t.\t"h }; \
h = substr($0,2); s="" } else { s = s$0 }} \
END{ print h"\t.\ttranscript\t1\t"(length(s))"\t.\t+\t.\t"h }' > nigoni.cds.small.gff
paste briggsae.cds.small.gff nigoni.cds.small.gff | cut -f9,18 | tr '\t' ',' | sed 's/,/__/' > gene_names
cat gene_names | sed 's/^/gene_id=/' | paste briggsae.cds.small.gff - | cut -f1-8,10 > tmp && mv tmp briggsae.cds.small.gff
cat gene_names | sed 's/^/gene_id=/' | paste nigoni.cds.small.gff - | cut -f1-8,10 > tmp && mv tmp nigoni.cds.small.gff

# generate bed files
cat briggsae.cds.small.gff | awk '{ print $1"\t"($4-1)"\t"$5"\t"(substr($9,9))}' > briggsae.cds.small.bed
cat nigoni.cds.small.gff | awk '{ print $1"\t"($4-1)"\t"$5"\t"(substr($9,9))}' > nigoni.cds.small.bed

# competitive mapping
# mkdir split_compmap
# cd split_compmap
# ls ../sim_read_data/ | parallel --plus python "~/Google\\ Drive/GitHub/CompMap/CompMap.py" -1 ../alignments/{..}_Cbr.sorted.bam -2 ../alignments/{..}_Cni.sorted.bam -r ../readnames/{..}.txt --base {..} --AS_tag AS --NM_tag NM
# cd ..

# competitive mapping counts
mkdir counts
cd counts
ls ../sim_read_data/hf1_split/ | grep Cbr | head -3 | sed 's/_Cbr.*//' | parallel \
"CompMap -1 ../alignments/{}_Cbr.sorted.bam -2 ../alignments/{}_Cni.sorted.bam -b1 ../briggsae.cds.small.bed -b2 ../nigoni.cds.small.bed --base {} --NM_tag NM -s1 Cbr -s2 Cni"
paste HF1_*_counts.txt | cut -f1-3,5,6,8,9 > compmap_counts_HF1.txt
cd ../

# featureCounts
featureCounts -t transcript -a briggsae.cds.small.gff -o counts/Cbr_fc.counts.txt alignments/Cbr_*_Cbr.sorted.bam alignments/hf1_split/HF1_*_Cbr.sorted.bam
featureCounts -t transcript -a nigoni.cds.small.gff -o counts/Cni_fc.counts.txt alignments/Cni_*_Cni.sorted.bam alignments/hf1_split/HF1_*_Cni.sorted.bam
cd counts
paste Cbr_fc.counts.txt Cni_fc.counts.txt | cut -f1,7-12,19-24 | tail -n +2 | sed -E 's/alignments\/|hf1_split\/|\.sorted\.bam//g' > real_counts.txt
cat real_counts.txt | cut -f1,5-7,11-13 > real_counts_HF1.txt
cat real_counts.txt | cut -f1,2-4,8-10 > real_counts_species.txt

# htseq counts
# htseq-count --secondary-alignments ignore --nonunique none --supplementary-alignments ignore -n 4 -t transcript --idattr gene -c counts/countsCbr.txt split_compmap/*_Cbr.bam briggsae.cds.small.gff
# htseq-count --secondary-alignments ignore --nonunique none --supplementary-alignments ignore -n 4 -t transcript --idattr gene -c counts/countsCni.txt split_compmap/*_Cni.bam nigoni.cds.small.gff
# htseq-count --secondary-alignments ignore --nonunique none --supplementary-alignments ignore -n 4 -t transcript --idattr gene -c counts/amb_countsCbr.txt split_compmap/ambiguous/*_Cbr.bam briggsae.cds.small.gff
# htseq-count --secondary-alignments ignore --nonunique none --supplementary-alignments ignore -n 4 -t transcript --idattr gene -c counts/amb_countsCni.txt split_compmap/ambiguous/*_Cni.bam nigoni.cds.small.gff
# htseq-count --secondary-alignments ignore --nonunique none --supplementary-alignments ignore -n 4 -t transcript --idattr gene -c counts/op_countsCbr.txt split_compmap/other/*_Cbr.bam briggsae.cds.small.gff
# htseq-count --secondary-alignments ignore --nonunique none --supplementary-alignments ignore -n 4 -t transcript --idattr gene -c counts/op_countsCni.txt split_compmap/other/*_Cni.bam nigoni.cds.small.gff
# echo -e "\tCbr_1\tCbr_2\tCbr_3\tHF1_Cbr_1\tHF1_Cbr_2\tHF1_Cbr_3" > head
# cat head counts/countsCbr.txt > tmp && mv tmp counts/countsCbr.txt
# echo -e "\tCni_1\tCni_2\tCni_3\tHF1_Cni_1\tHF1_Cni_2\tHF1_Cni_3" > head
# cat head counts/countsCni.txt > tmp && mv tmp counts/countsCni.txt
# echo -e "\tamb_Cbr_1\tamb_Cbr_2\tamb_Cbr_3\tamb_HF1_Cbr_1\tamb_HF1_Cbr_2\tamb_HF1_Cbr_3" > head
# cat head counts/amb_countsCbr.txt > tmp && mv tmp counts/amb_countsCbr.txt
# echo -e "\tamb_Cni_1\tamb_Cni_2\tamb_Cni_3\tamb_HF1_Cni_1\tamb_HF1_Cni_2\tamb_HF1_Cni_3" > head
# cat head counts/amb_countsCni.txt > tmp && mv tmp counts/amb_countsCni.txt
# echo -e "\tCbr_1_Cni\tCbr_2_Cni\tCbr_3_Cni" > head
# cat head counts/op_countsCbr.txt > tmp && mv tmp counts/op_countsCbr.txt
# echo -e "\tCni_1_Cbr\tCni_2_Cbr\tCni_3_Cbr" > head
# cat head counts/op_countsCni.txt > tmp && mv tmp counts/op_amb_countsCni.txt
# rm head
# # real counts
# htseq-count --secondary-alignments ignore --nonunique none --supplementary-alignments ignore -n 4 -t transcript --idattr gene -c counts/real_countsCbr.txt alignments/Cbr_*_Cbr.sorted.bam alignments/hf1_split/*_Cbr.sorted.bam briggsae.cds.small.gff
# htseq-count --secondary-alignments ignore --nonunique none --supplementary-alignments ignore -n 4 -t transcript --idattr gene -c counts/real_countsCni.txt alignments/Cni_*_Cni.sorted.bam alignments/hf1_split/*_Cni.sorted.bam nigoni.cds.small.gff
# echo -e "\tCbr_1\tCbr_2\tCbr_3\tHF1_Cbr_1\tHF1_Cbr_2\tHF1_Cbr_3" > head
# cat head counts/real_countsCbr.txt > tmp && mv tmp counts/real_countsCbr.txt
# echo -e "\tCni_1\tCni_2\tCni_3\tHF1_Cni_1\tHF1_Cni_2\tHF1_Cni_3" > head
# cat head counts/real_countsCni.txt > tmp && mv tmp counts/real_countsCni.txt
# rm head

# ortholog pairs
cat briggsae.cds.small.fa nigoni.cds.small.fa | grep ">" | sed 's/>//' | pr -2 -t -s, -l2000 > orthologs.small.txt
