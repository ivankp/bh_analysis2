#!/bin/bash

merger=~/work/bh_analysis/bin/merge_parts

if [ "$1" == "--all" ]; then
  for raw in *_raw.root; do
    cooked=$(sed 's/_raw//' <<< $raw)
    if [ ! -s "$cooked" ]; then
      $merger $cooked $raw
    fi
  done
fi

for name in $(ls | sed -n 's/\(.*\)_\(B\|RS\|I\|V\)_raw\.root$/\1/p' | sort -u)
do
  if [ ! -s "${name}_NLO.root" ]; then
    $merger ${name}_NLO.root ${name}_*_raw.root
  fi
done
