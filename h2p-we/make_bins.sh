#!/bin/bash
# make linear bins
# 1 = lower bound
# 2 = interval
# 3 = upper bound

printf "\n['-inf'"
for bin in $(seq $1 $2 $3) ; do
    printf " ,$bin"
done
printf ", 'inf']\n\n"
