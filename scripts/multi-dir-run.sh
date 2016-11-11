#!/bin/bash

#
# This script executes a given python script (as input) within different spin directions (i.e. spin-*/)
#

script="$1"

if [ -z "$1" ]
  then
    echo "Please supply the python script you want to run, as an argument."
    exit 1
fi

#echo $script

OWD=$( pwd )
echo $OWD              # Original working directory
for i in $( ls -d spin*/ ); do
   cd $i
   echo i am here: $( pwd )
   python ../$script
   cd $OWD
done
