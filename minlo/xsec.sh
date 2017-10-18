#!/bin/bash

dir_w=/msu/data/t3work8/ivanp/scale_dep/HT2
for j in 3; do
  dir_bh=/msu/data/t3work4/luisonig/H${j}jets_ggf/EDNTuplesFiles
  for p in B RS I V; do
    fs=$(ls ${dir_w}/H${j}.0j_GGFHT_${p}_6500_pt25.0_eta4.5_r100_*.root \
        | sed 's/.*\/\(.*\)/\1/')
    echo $fs
    /home/ivanp/work/bh_analysis2/bin/xsec_Hjets ${j}j \
      $(echo $fs | tr ' ' '\n' | sed "s|^|${dir_bh}/|") \
      w $(echo $fs | tr ' ' '\n' | sed "s|^|${dir_w}/|") \
      > minlo_H${j}j_${p}.txt
  done
done
