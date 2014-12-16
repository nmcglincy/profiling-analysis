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
grep -w "Covering one intact ORF" S2_tcd_mTIFAnno.txt > covering-one-intact-orf

awk '{ for (i=1; i<=6; i++) printf("%s%s", $(i), i<6 ? OFS : "\n"); }'



grep -w "Intergenic transcripts" S2_tcd_mTIFAnno.txt > intergenic-trcpts
grep -w "Overlap 3' of one ORF" S2_tcd_mTIFAnno.txt > overlap-3p-oneOrf
grep -w "Overlap 5' of one ORF" S2_tcd_mTIFAnno.txt > overlap-5p-oneOrf
# 
# THE INTERGENIC TRANSCRIPTS DON'T SEEM TO HAVE ANY ID, IT WOULD BE USEFUL IF THEY DID.


# FUSE INTO FILES BY FLAVOUR
cat covering-one-intact-orf intergenic-trcpts > mTIF-anno-vanilla
cat covering-one-intact-orf intergenic-trcpts overlap-5p-oneOrf overlap-3p-oneOrf > mTIF-anno-choco
# 
# TODO - LOCATION DATA IS IN THE WRONG FORMAT
# TODO - HAVE 5PUTR & 3PUTR AS SEPARATE ANNOTATIONS IN GTF
# 