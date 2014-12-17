#!/bin/awk -f
# AN AWK SCRIPT TO CHANGE LOCATION DETAILS TO MAKE IT CONSISTENT WITH NICK'S GTF.
# SMALLER COORDINATE NEEDS TO BE FIRST EVEN IF IT'S ON THE NEGATIVE STRAND
{
    if ( $2 == "-") {
        print "hello"
    } else {
        print "goodbye"
    }
}