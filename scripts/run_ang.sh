#!/bin/bash

if [ -z "$1" ]
then
  echo "usage: $0 number_of_jets"
  exit
fi

# ./bin/hist_Hjets_ang -j$1 -o ang_H${1}j.root \
./bin/var_Hjets_angles -j$1 -o H${1}j_angles.root \
  $(sqlite3 sql/ntuples.db -separator "/" \
    "select distinct dir,file from ntuples where instr(info,\"mtop\") and njets=$1 and part=\"B\"")
