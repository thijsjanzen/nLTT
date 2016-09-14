#!/bin/bash
# This script satisfies @lintr-bot about this package its R code

# Replace all tabs by two spaces, from http://stackoverflow.com/a/11094620/3364162
find . -name '*.R' ! -type d -exec bash -c 'expand -t 2 "$0" > /tmp/e && mv /tmp/e "$0"' {} \;
