#!/bin/bash

if [[ $# -ne 1 ]]; then
  echo "usage: $0 script"
  exit 1
fi

base="${1%.*}"

if [ ! -x "$1" ]; then
  echo "File \"$1\" does not exists or is not executable"
  exit 1
fi

echo "
Universe   = vanilla
Executable = $1
Output     = ${base}.out
Error      = ${base}.err
Log        = ${base}.log
getenv = True
Queue
" | condor_submit - > /dev/null
# " > ${base}.condor

