#!/bin/bash

for f in `ls | grep 'njets_test_.*\.root'`; do

rxplot -i $f -o `echo $f | sed 's/\.root$/_w.pdf/'` \
  -r 'sng/w_(particles|jets_[^_]+)(_(.*))?/\3/draw=histtext0' \
     'nl/w_(particles|jets_[^_]+).*/\1' \
     'x/.*/njets' 'y/.*/weight' \
     'gt' 't/_/, ' 't/nocut/no cuts' 't/pt([0-9]+)/pT > \1 GeV' 't/y([0-9]+)([0-9])/y < \1.\2' \
  -ltr --ytitle-offset=1.3

rxplot -i $f -o `echo $f | sed 's/\.root$/_n.pdf/'` \
  -r 'sng/n_(particles|jets_[^_]+)(_(.*))?/\3/draw=histtext0' \
     'nl/n_(particles|jets_[^_]+).*/\1' \
     'x/.*/njets' 'y/.*/num entries' \
     'gt' 't/_/, ' 't/nocut/no cuts' 't/pt([0-9]+)/pT > \1 GeV' 't/y([0-9]+)([0-9])/y < \1.\2' \
  -ltr --ytitle-offset=1.3

done
