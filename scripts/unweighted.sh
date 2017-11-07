#!/bin/bash

# max_weight[1]=1876.79
# max_weight[2]=17363.2
# max_weight[3]=23848.6

max_weight[1]=1000.00
max_weight[2]=3162.28
max_weight[3]=3162.28

./bin/unweighted -j$1 \
  ${max_weight[$1]+-w${max_weight[$1]}} \
  ${max_weight[$1]+-oH${1}j_mtop_unweighted.root} \
  $(sqlite3 sql/ntuples.db -separator "/" \
    "select distinct dir,file from ntuples where instr(info,\"mtop\") and njets=$1 and part=\"B\"")
