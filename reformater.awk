#!/bin/awk -f
BEGIN {OFS="\t"}
{
	for (i=1; i<=6; i++) printf("%s%s", $(i), OFS);
	for (i=7; i<=NF-1; i++) printf("%s%s", $(i), i<NF-1 ? "_" : OFS)
    print $NF
}