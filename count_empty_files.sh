#!/bin/bash
#
# check recursively a directory:
#
function checkDirr() {
  for x in $1/*
  do

    #
    # if it's a sub-directory:
    #
    if [ -d $x ]; then
      #echo "$x is a directory."
      ndirs=$((ndirs + 1))
      #
      # glob all entries in an array:
      #
      entries_array=($x/*)
      #
      # check entries array count:
      #
      if [ ${#entries_array[*]} -eq 0 ]; then
        #echo "empty dir"
        nemptydirs=$((nemptydirs + 1))
      else
        #
        # call recursively myself to check non-empty sub-directory:
        #
        checkDirr $x

      fi

    #
    # otherwise it's a file:
    #
    else
      #echo "$x is a file."
      nfiles=$((nfiles + 1))
      #
      # check empty file:
      #
      if [ ! -s $x ]; then
        nemptyfiles=$((nemptyfiles + 1))
      fi
    fi

  done
}

#
# main():
#
#
# initialize globals:
#
ndirs=0
nfiles=0

nemptydirs=0
nemptyfiles=0

#
# use given directory path or the current directory .:
#
dir0=${1:-.}

#
# set shell nullglob option to avoid dir/* string when dir is empty:
#
shopt -s nullglob

#
# now check recursively the directory:
#
checkDirr $dir0

#
# unset shell nullglob option:
#
shopt -u nullglob

#
# send statistics:
#
echo $dir0 has:
echo
echo $nemptyfiles empty files
echo $((nfiles - nemptyfiles)) files with data
echo $nemptydirs empty directories
echo $((ndirs - nemptydirs)) non-empty directories
