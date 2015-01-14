# From inside /mnt/ingolialab/mcglincy/NINM001/Sample_NINM01_index1

fastx-split -o split-20150113 -p NN -x NNNNNIIIII --min-insert=14 -s samples.csv sample_clipped.fastq

cd split-20150113

bowtie2 -p 8 --very-sensitive --quiet --un 112A_norrna.fq -x /mnt/ingolialab/mcglincy/yeast_annotation/sc-rrna -U 112A.fastq | rrna-stats -o 112A_rrna --tam --maxread 51 --lenrange 24,36 - &
bowtie2 -p 8 --very-sensitive --quiet --un 112B_norrna.fq -x /mnt/ingolialab/mcglincy/yeast_annotation/sc-rrna -U 112B.fastq | rrna-stats -o 112B_rrna --tam --maxread 51 --lenrange 24,36 - &
bowtie2 -p 8 --very-sensitive --quiet --un D19A_norrna.fq -x /mnt/ingolialab/mcglincy/yeast_annotation/sc-rrna -U D19A.fastq | rrna-stats -o D19A_rrna --tam --maxread 51 --lenrange 24,36 - &
bowtie2 -p 8 --very-sensitive --quiet --un D19B_norrna.fq -x /mnt/ingolialab/mcglincy/yeast_annotation/sc-rrna -U D19B.fastq | rrna-stats -o D19B_rrna --tam --maxread 51 --lenrange 24,36 - & 

nohup tophat -G /mnt/ingolialab/mcglincy/yeast_annotation/sac_cer_yassour.gtf --transcriptome-index=transcriptome_data/known /mnt/ingolialab/mcglincy/yeast_annotation/saccharomyces_cerevisiae
tophat -o 112A_norrna_v_genome -p4 --no-novel-juncs -G /mnt/ingolialab/mcglincy/yeast_annotation/sac_cer_yassour.gtf --transcriptome-index=transcriptome_data/known -M /mnt/ingolialab/mcglincy/yeast_annotation/saccharomyces_cerevisiae 112A_norrna.fq &
tophat -o 112B_norrna_v_genome -p4 --no-novel-juncs -G /mnt/ingolialab/mcglincy/yeast_annotation/sac_cer_yassour.gtf --transcriptome-index=transcriptome_data/known -M /mnt/ingolialab/mcglincy/yeast_annotation/saccharomyces_cerevisiae 112B_norrna.fq &
tophat -o D19A_norrna_v_genome -p4 --no-novel-juncs -G /mnt/ingolialab/mcglincy/yeast_annotation/sac_cer_yassour.gtf --transcriptome-index=transcriptome_data/known -M /mnt/ingolialab/mcglincy/yeast_annotation/saccharomyces_cerevisiae D19A_norrna.fq &
tophat -o D19B_norrna_v_genome -p4 --no-novel-juncs -G /mnt/ingolialab/mcglincy/yeast_annotation/sac_cer_yassour.gtf --transcriptome-index=transcriptome_data/known -M /mnt/ingolialab/mcglincy/yeast_annotation/saccharomyces_cerevisiae D19B_norrna.fq &

samtools view -h 112A_norrna_v_genome/accepted_hits.bam | grep -E '(NM:i:0)|(^@)' | samtools view -S -b - > 112A_norrna_v_genome.bam &
samtools view -h 112B_norrna_v_genome/accepted_hits.bam | grep -E '(NM:i:0)|(^@)' | samtools view -S -b - > 112B_norrna_v_genome.bam &
samtools view -h D19A_norrna_v_genome/accepted_hits.bam | grep -E '(NM:i:0)|(^@)' | samtools view -S -b - > D19A_norrna_v_genome.bam &
samtools view -h D19B_norrna_v_genome/accepted_hits.bam | grep -E '(NM:i:0)|(^@)' | samtools view -S -b - > D19B_norrna_v_genome.bam &

# OK, SO THIS PICS PERFECT ALIGNMENTS, BUT DOES NOT EXCLUDE MULTI-MAPPED READS

samtools index 112A_norrna_v_genome.bam &
samtools index 112B_norrna_v_genome.bam &
samtools index D19A_norrna_v_genome.bam &
samtools index D19B_norrna_v_genome.bam &

mkdir Statistics
fp-framing -o Statistics/112A_norrna_v_genome -b /mnt/ingolialab/mcglincy/yeast_annotation/sac_cer_yassour.bed 112A_norrna_v_genome.bam &
fp-framing -o Statistics/112B_norrna_v_genome -b /mnt/ingolialab/mcglincy/yeast_annotation/sac_cer_yassour.bed 112B_norrna_v_genome.bam &
fp-framing -o Statistics/D19A_norrna_v_genome -b /mnt/ingolialab/mcglincy/yeast_annotation/sac_cer_yassour.bed D19A_norrna_v_genome.bam &
fp-framing -o Statistics/D19B_norrna_v_genome -b /mnt/ingolialab/mcglincy/yeast_annotation/sac_cer_yassour.bed D19B_norrna_v_genome.bam &

fp-count -o Statistics/112A_norrna_v_genome_qexpr.txt -b /mnt/ingolialab/mcglincy/yeast_annotation/sac_cer_yassour.bed -a Statistics/asites.txt 112A_norrna_v_genome.bam &
fp-count -o Statistics/112B_norrna_v_genome_qexpr.txt -b /mnt/ingolialab/mcglincy/yeast_annotation/sac_cer_yassour.bed -a Statistics/asites.txt 112B_norrna_v_genome.bam &
fp-count -o Statistics/D19A_norrna_v_genome_qexpr.txt -b /mnt/ingolialab/mcglincy/yeast_annotation/sac_cer_yassour.bed -a Statistics/asites.txt D19A_norrna_v_genome.bam &
fp-count -o Statistics/D19B_norrna_v_genome_qexpr.txt -b /mnt/ingolialab/mcglincy/yeast_annotation/sac_cer_yassour.bed -a Statistics/asites.txt D19B_norrna_v_genome.bam &

# 
# Move to R script