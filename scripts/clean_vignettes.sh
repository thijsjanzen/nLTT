#!/bin/bash
# This script satisfies @lintr-bot about the knitr intermediate R code
#
# When knitr creates vignettes, it creates R code as an intermediate.
# lintr-bot does not like those intermediates by default.
# For example, knitr adds lines like `## ---` which are
# interpreted as commented code. 

cd inst/doc

for filename in `ls *.R` 
do
  echo $filename

  # Replace 'commented code'
  sed -i -e 's/## ---/# ---/g' $filename

  # Remove trailing newlines, from http://stackoverflow.com/a/7359879/3364162
  sed -i -e :a -e '/^\n*$/{$d;N;};/\n$/ba' $filename

  # Remove trailing whitespace, from http://stackoverflow.com/a/4438318/3364162
  sed -i 's/[ \t]*$//' $filename

done
