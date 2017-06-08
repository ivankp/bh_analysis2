#!/bin/bash

DIR=$(pwd -P)
# OUT='/msu/data/t3work5/ivanp/4/out'
OUT='/msu/data/t3work8/ivanp/minlo4'
DB=$DIR/../sql/ntuples.db

max_jobs=100

for ntuple in $(sqlite3 $DB "SELECT dir,file,njets,part,sid FROM ntuples WHERE 
particle='H' and energy=13 and instr(info,'ED') and njets=3 and sid=100")
do

  # wait if too many jobs
  while [ "$(condor_q ivanp | tail -1 |\
           sed 's/\([0-9]*\) jobs.*/\1/')" -ge "$max_jobs" ]
  do sleep 1800; done

  IFS='|'; read dir file njets part sid <<< "$ntuple"

  base=$(basename $file .root)
  name=${OUT}/$base
  script=${name}.sh
  # pipe=${name}.wts.gz

  # skip if file exists
  if [ -a "${name}.wts.gz" ]; then continue; fi
  # if [ -a "${OUT}/../minlo/${base}.wts.gz" ]; then continue; fi

  echo $file

  # write script
  echo "#!/bin/bash" > $script
  chmod +x $script
  echo "cd $OUT" >> $script
  # echo "mkfifo $pipe" >> $script

  echo "/msu/data/t3work5/ivanp/4/sherpa/bin/Sherpa -f $DIR/Run2.dat" \
       " \\\$NJETS:=$njets \\\$NLOPART:=$part \\\$SID:=$sid" >> $script
  # echo "sleep 5" >> $script

  # echo "zcat $pipe | /home/ivanp/work/sherpa_tools/bin/wts2root ${name}.root" >> $script
  # echo "rm $pipe" >> $script

echo "
Universe   = vanilla
Executable = $script
Error      = ${name}.err
Output     = ${name}.out
Log        = ${name}.log
getenv = True
Queue
" | condor_submit - > /dev/null

  sleep 0.1

done
