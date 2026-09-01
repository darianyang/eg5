#!/bin/bash

filein=$1
fileout=$2

paste <(cat $filein | tail -n +2 | awk {'print $10'}) <(cat $filein | tail -n +2 | awk {'print $19'}) <(cat $filein | tail -n +2 | awk {'print $28'}) > $fileout
