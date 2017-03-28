#!/bin/bash

./bin/hist_Hjets_mtop ${1}j out:finite_raw_${1}j.root \
  $(sed "s|.*|/msu/data/t3work4/luisonig/H${1}jets_ggf_mtop/NTuplesFiles/&|" ~ivanp/finite_H${1}j_ok.txt)

~ivanp/work/bh_analysis/bin/merge_parts finite_${1}j.root finite_raw_${1}j.root
