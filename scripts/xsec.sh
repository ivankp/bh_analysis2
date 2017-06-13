#!/bin/bash

nj=$1

for P in B RS I V; do
echo $P
for R in 1 2 3 4 5 6 7 8 9 10; do
echo $R
  ./bin/xsec_Hjets ${nj}j antikt${R} \
    /msu/data/t3work4/luisonig/H${nj}jets_ggf/NTuplesFiles/*_${P}_6500*10{0,1,2,3,4}.root \
    w ~/disk2/wt/${nj}_H${nj}j/H${nj}j_13TeV_${P}_10{0,1,2,3,4}_HT2-unc_CT10nlo.root \
    > H${nj}j_AntiKt${R}_${P}.out
done
done
echo "DONE"
