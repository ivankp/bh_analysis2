#!/bin/bash

pwd -P

./bin/hist_Hjets_mtop ${1}j out:inf_${1}j_raw.root \
  /msu/data/t3work4/luisonig/H${1}jets_ggf/NTuplesFiles/H${1}.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_10*.root

~ivanp/work/bh_analysis/bin/merge_parts inf_${1}j.root inf_${1}j_raw.root
