#!/bin/bash
# Convert a Config.sh file to an Arepoflags.mk file
#
# Takes each line starting with an alphanumeric character and prepends it with
# 'AREPOFLAGS += -D'.
#
# Typical usage is:
#     bash configToArepoflags.sh Config.sh > Arepoflags.mk

if [ "$#" -ne 1 ]
then
  echo "Usage: bash $0 input-file > output-file"
  exit 1
fi

sed 's/\(^[[:alpha:]].*\)/AREPOFLAGS += -D\1/' $1
