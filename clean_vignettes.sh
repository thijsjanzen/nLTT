#!/bin/bash
cd inst/doc

for filename in `ls *.R` 
do
  echo $filename
  sed -i -e 's/## ---/# ---/g' $filename
done
