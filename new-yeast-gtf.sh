# CREATING A NEW YEAST GTF FILE BASED ON A FUSION OF NICK'S MOST RECENT GTF AND INFORMATION FROM THE STEINMETZ TIF-SEQ PAPER, PELECHANO ET AL. (2013) NATURE
# 
# CREATE TWO VERSIONS:
# 1. "VANILLA"
# mTIFs COVERING ONE INTACT ORF
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