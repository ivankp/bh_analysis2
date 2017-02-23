#!/bin/bash

./bin/bh_analysis2/bin/hist_Hjets_mtop 2j out:finite_raw.root \
  `sed 's|.*|/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/&|' ~ivanp/finite_H2j_ok.txt`

~ivanp/work/bh_analysis/bin/merge_parts finite.root finite_raw.root
