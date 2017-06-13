#!/bin/bash

for f in ATLASdata/NNLO*; do
  sed 's/.*\/.\+\.\(.\+\)\.txt/\1 {/' <<< $f | tr '\n' ' '
  <$f cut -d' ' -f1 | xargs printf '%f\n' | sed 's/\.\?0\+$//' | tr '\n' ' '
  <$f tail -1 | cut -d' ' -f3 | xargs printf '%f\n' | sed 's/\.\?0\+$/ }/'
done
