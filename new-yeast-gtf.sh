# CREATING A NEW YEAST GTF FILE BASED ON A FUSION OF NICK'S MOST RECENT GTF AND INFORMATION FROM THE STEINMETZ TIF-SEQ PAPER, PELECHANO ET AL. (2013) NATURE
# 
# CREATE TWO VERSIONS:
# 1. "VANILLA"
# mTIFs COVERING ONE INTACT ORF
# mTIFs INTERGENIC TRANSCRIPTS
# 2. "CHOCOLATE"
# mTIFs COVERING ONE INTACT ORF
# mTIFs INTERGENIC TRANSCRIPTS
# mTIFs OVERLAP 5P OF ONE ORF
# mTIFs OVERLAP 3P OF ONE ORF
# 
# COPY NICK'S OLD GTF TO MY GITHUB REPO
# /Users/nmcglincy/Documents/computing/github/nmcglincy/profiling-analysis
# 
scp mcglincy@compute1.ingolia-lab.org:/mnt/ingolialab/ingolia/Genomes/Saccharomyces_cerevisiae/YeastGenome/sac_cer_yassour.gtf .
# 
# FOR THE mTIFs WE WANT PELECHANO ET AL. (2013) SUPP. DATA 2
# S2_tcd_mTIFAnno.txt
# 
# MANUALLY REMOVED COLUMN HEADERS TO ANOTHER FILE mTIF-anno-headers
# 
# GRAB mTIF TYPES INTO SEPARATE FILES USING GREP
grep -w "Covering one intact ORF" S2_tcd_mTIFAnno.txt | awk -f reformater.awk > covering-one-intact-orf
grep -w "Overlap 3' of one ORF" S2_tcd_mTIFAnno.txt | awk -f reformater.awk > overlap-3p-oneOrf
grep -w "Overlap 5' of one ORF" S2_tcd_mTIFAnno.txt | awk -f reformater.awk > overlap-5p-oneOrf
# 
# THIS IS A SPECIAL CASE, SINCE WE NEED TO CREATE IDS FOR THE INTERGENIC TRANSCRIPTS
grep -w "Intergenic transcripts" S2_tcd_mTIFAnno.txt | awk -f reformater-igt-tifs.awk > intergenic-trcpts
# 
# FUSE INTO FILES BY FLAVOUR
cat covering-one-intact-orf intergenic-trcpts > mTIF-anno-vanilla
cat covering-one-intact-orf intergenic-trcpts overlap-5p-oneOrf overlap-3p-oneOrf > mTIF-anno-choco
# 
# ONLY TAKE mTIFs WITH SOME EXPRESSION IN THE YPD CONDITION; REMOVE FIELD WITH EXPRESSION IN YGAL CONDITION
# CAN'T BRING MYSELF TO REMOVE EXPRESSION IN YPD CONDITION YET, MIGHT BE USEFUL LATER ON, BUT I'LL HAVE TO 
# SEE WHETHER THERE'S A PLACE FOR IT IN THE FINAL GTF.
awk '{ if($5>=1) {print $0} }' mTIF-anno-vanilla | cut -f 6 > mTIF-anno-vanilla-ypd
awk '{ if($5>=1) {print $0} }' mTIF-anno-choco | cut -f 6 > mTIF-anno-choco-ypd
# 
# TODO - I NOTICED THAT SOME OF THE mTIFs FROM ALL CLASSES HAVE NA AS A GENE IDENTIFIER, I'LL HAVE TO DO SOME THING ABOUT THAT.
# TODO - LOCATION DATA IS IN THE WRONG FORMAT.
# TODO - HAVE 5PUTR & 3PUTR AS SEPARATE ANNOTATIONS IN GTF.
# 