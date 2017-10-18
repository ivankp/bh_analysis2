#!/bin/bash

ntuples=/msu/data/t3work4/luisonig/H3jets_ggf/EDNTuplesFiles
weights=/msu/data/t3work8/ivanp/weights/H3j/minlo-HT2-CT10nlo
out=$(sed 's:/H[0-9]j/minlo-:/minlo/:' <<< $weights)

mkdir -p $out

i=0
diff <(ls $ntuples) <(ls $weights) | sed -n 's/< //p' \
| while read f; do
  echo $f
  base="${f%.*}"
  scale=$(sed 's:.*/minlo-\(HT[12]\)-.*:\\$SCALE_DEF\:=\\$SCALE_\1:' <<< $weights)
  njets=$(sed 's/H\([1-9]\)\.0j_.*/\\$NJETS:=\1/' <<< $f)
  nlopart=$(sed 's/.*_\(B\|RS\|I\|V\)_.*/\\$NLOPART:=\1/' <<< $f)
  sid=$(sed 's/.*_\([0-9]\+\)\.root/\\$SID:=\1/' <<< $f)

  script=${out}/${base}.sh
  echo "#!/bin/bash
cd $out
/msu/data/t3work5/ivanp/4/sherpa/bin/Sherpa \\
  -f /home/ivanp/work/bh_analysis2/minlo/Run.dat \\
  \\\$OUT_PATH:=$out \\
  $njets $nlopart $sid $scale
" > $script
chmod +x $script

  echo "
Universe   = vanilla
Executable = $script
Output     = ${out}/${base}.out
Error      = ${out}/${base}.err
Log        = ${out}/${base}.log
getenv = True
Queue
" | condor_submit - > /dev/null

  if ((++i >= 100)); then break; fi
  sleep 0.1

done
