#!/bin/bash

source ~/.bash_profile
source ~/tmp/cantera/usr/local/bin/setup_cantera

files=$(find . -type f -name "*.xml")
echo $files
for f in $files ; do
   echo "$f"
   cmd="python /Users/james/tmp/cantera/usr/local/lib/python3.9/site-packages/cantera/ctml2yaml.py $f"
   #echo $cmd
   $cmd
done

