#!/bin/bash
cat $1 | awk -v OFS="\t" '{ split($9,a,";"); print $1,($4-1),$5,a[1] }' 
