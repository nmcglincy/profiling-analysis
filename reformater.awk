#!/bin/awk -f
BEGIN {OFS="\t"}
{
	for (i=1; i<=6; i++) printf("%s%s", $(i), i<6 ? OFS : "\n");
	for (i=7; i<=NF-1; i++) printf("%s%s", $(i), i<NF-1 ? "_" : "\n")
#    print NF-1
    print $NF
}