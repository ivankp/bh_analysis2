#!/bin/bash

if [[ $# -ne 1 ]]; then
  echo "usage: $0 script.sh"
  exit 1
fi

base="${1%.*}"
ext="${1##*.}"

if [ "$ext" != "sh" ]; then
  echo "Argument must have .sh extension"
  exit 1
fi

if [ ! -x "${base}.sh" ]; then
  echo "File \"${base}.sh\" does not exists or is not executable"
  exit 1
fi

echo "
Universe   = vanilla
Executable = ${base}.sh
Output     = ${base}.out
Error      = ${base}.err
Log        = ${base}.log
getenv = True
Queue
" > ${base}.condor
#| condor_submit - > /dev/null

