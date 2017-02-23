#!/bin/bash

./bin/bh_analysis2/bin/hist_Hjets_mtop 2j out:inf_raw.root \
  /msu/data/t3work4/luisonig/H2jets_ggf/NTuplesFiles/H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_10*.root

~ivanp/work/bh_analysis/bin/merge_parts inf.root inf_raw.root
