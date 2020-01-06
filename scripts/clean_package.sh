#!/bin/bash
# This script satisfies @lintr-bot about this package its R code

# Replace all tabs by two spaces, from http://stackoverflow.com/a/11094620/3364162
find . -name '*.R' ! -type d -exec bash -c 'expand -t 2 "$0" > /tmp/e && mv /tmp/e "$0"' {} \;

cd R

for filename in `ls *.R` 
do
  echo $filename

  # Remove trailing newlines, from http://stackoverflow.com/a/7359879/3364162
  sed -i -e :a -e '/^\n*$/{$d;N;};/\n$/ba' $filename

  # Remove trailing semicolons
  sed -i 's/;$//' $filename

  # Remove trailing whitespace, from http://stackoverflow.com/a/4438318/3364162
  sed -i 's/[ \t]*$//' $filename

done

cd ../tests/testthat

for filename in `ls *.R` 
do
  echo $filename

  # Remove trailing newlines, from http://stackoverflow.com/a/7359879/3364162
  sed -i -e :a -e '/^\n*$/{$d;N;};/\n$/ba' $filename

  # Remove trailing semicolons
  sed -i 's/;$//' $filename

  # Remove trailing whitespace, from http://stackoverflow.com/a/4438318/3364162
  sed -i 's/[ \t]*$//' $filename

done
