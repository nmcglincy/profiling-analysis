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
awk '{ if($5>=1) {print $0} }' mTIF-anno-vanilla | cut -f 1-5,7- > mTIF-anno-vanilla-ypd
awk '{ if($5>=1) {print $0} }' mTIF-anno-choco | cut -f 1-5,7- > mTIF-anno-choco-ypd
# 
# TODO - LOCATION DATA IS IN THE WRONG FORMAT.
awk -f reorder-loc.awk mTIF-anno-vanilla-ypd > mTIF-anno-vanilla-ypd1
# 
# THIS WILL BE 0 IS EVERYTHING WENT FINE
awk -f loc-double-check.awk mTIF-anno-vanilla-ypd1 | grep -cw "oh shit"
# 
# FOR THE CHOCOLATE ANNOTATION:
awk -f reorder-loc.awk mTIF-anno-choco-ypd > mTIF-anno-choco-ypd1
awk -f loc-double-check.awk mTIF-anno-choco-ypd1 | grep -cw "oh shit"
# 
# DOUBLE DOUBLE CHECK
wc -l mTIF-anno-choco*
# 
# TODO - I NOTICED THAT SOME OF THE mTIFs FROM ALL CLASSES HAVE NA AS A GENE IDENTIFIER, I'LL HAVE TO DO SOME THING ABOUT THAT.
grep -wc "NA" mTIF-anno-vanilla-ypd1
# 6279 ENTRIES, SEEM TO BE EXCLUSIVELY FROM COVERING-ONE-INTACT-ORF CLASS
# PUT IT INTO BED FORMAT
grep -w "NA" mTIF-anno-vanilla-ypd1 > orphans
awk '{print "chr" $1 "\t" $3 "\t" $4}' orphans > orphans.bed
# 
# GROAN, THE CHR NUMBERS ARE SPECIFIED IN LATIN NUMERAL RATHER THAN NUMBERS, SO I'M GETTING NO OVERLAP
sed -e 's/[[:<:]]chr1[[:>:]]/chrI/g' -e 's/[[:<:]]chr2[[:>:]]/chrII/g' -e 's/[[:<:]]chr3[[:>:]]/chrIII/g' -e 's/[[:<:]]chr4[[:>:]]/chrIV/g' -e 's/[[:<:]]chr5[[:>:]]/chrV/g' -e 's/[[:<:]]chr6[[:>:]]/chrVI/g' -e 's/[[:<:]]chr7[[:>:]]/chrVII/g' -e 's/[[:<:]]chr8[[:>:]]/chrVIII/g' -e 's/[[:<:]]chr9[[:>:]]/chrIX/g' -e 's/[[:<:]]chr10[[:>:]]/chrX/g' -e 's/[[:<:]]chr11[[:>:]]/chrXI/g' -e 's/[[:<:]]chr12[[:>:]]/chrXII/g' -e 's/[[:<:]]chr13[[:>:]]/chrXIII/g' -e 's/[[:<:]]chr14[[:>:]]/chrXIV/g' -e 's/[[:<:]]chr15[[:>:]]/chrXV/g' -e 's/[[:<:]]chr16[[:>:]]/chrXVI/g' orphans.bed > orphansRoman.bed
# REPEAT GALAXY ANALYSIS, WOULD BE COOL IF I COULD LEARN HOW TO DO THIS AT THE COMMAND LINE OR IN R
# 
# ACTUALLY PUT INTO GTF FORMAT
awk -f anno-to-gtf.awk mTIF-anno-vanilla-ypd1 > mTIF-vanilla-ypd.gtf
awk -f anno-to-gtf.awk mTIF-anno-choco-ypd1 > mTIF-choco-ypd.gtf
# 
# OK, AFTER ALL THE CONFUSION, I THINK THE BEST APPROACH IS TO MAKE AN EXTRA-VANILLA CONSISTING OF JUST THOSE COVERING-ONE-INTACT-ORF WITH A VALID GENE_ID 
# 
# MAKING GTF WITHOUT ANY NAS WITHOUT INTEGENIC TRANSCRIPTS
grep -wv "NA" mTIF-anno-vanilla-ypd1 | head
grep -wv "NA" mTIF-anno-vanilla-ypd1 | cut -f 7| uniq
# 
# THE ACTUAL COMMAND
grep -wv "NA" mTIF-anno-vanilla-ypd1 | grep -v "int_trcpt_*" > mTIF-exVanilla-ypd1
# 
# DOUBLE CHECK BY GREP C COUNTING
grep -cw "NA" mTIF-exVanilla-ypd1
grep -c "int_trcpt_*" mTIF-exVanilla-ypd1
# 
# IT ALL LOOKS GOOD
awk -f anno-to-gtf.awk mTIF-exVanilla-ypd1 > mTIF-exVanilla.gtf
# 
# TRY A DIFFERENT TAC FROM NICK'S NOTES, IMPORTING A BED FILE
awk -f anno-to-bed.gtf mTIF-exVanilla-ypd1 > mTIF-exVanilla.bed
