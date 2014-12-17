#!/bin/awk -f
# AN AWK SCRIPT TO CHANGE LOCATION DETAILS TO MAKE IT CONSISTENT WITH NICK'S GTF.
# SMALLER COORDINATE NEEDS TO BE FIRST EVEN IF IT'S ON THE NEGATIVE STRAND
BEGIN {OFS = "\t"}
{
    if ( $2 == "-") {
        print $1 $2 $4 $3 $5 $6 $7
    } else {
        print $0
    }
}