#!/bin/bash

exe=/home/ivanp/work/bh_analysis2/bin/hist_Hjets
dir=/msu/data/t3work2/ivanp
out=$dir/condor
histdir=$dir/hist/raw

group_size=50

mkdir -p $out
mkdir -p $histdir

db="$dir/ntuples.db"

#####################################################################

for jalg in AntiKt4
do

for set in `sqlite3 $db "
  SELECT distinct particle, njets, energy, part, wt.scales, wt.pdf, dset
  FROM bh
  JOIN wt ON bh.id = wt.bh_id
  WHERE particle = 'H' and energy = 13 and wt.pdf = 'CT10nlo' and wt.scales = 'HT2-unc' and dset = 3
"`; do

arr=(`echo $set | tr '|' ' '`)
arr[2]=`printf "%.0f" ${arr[2]}` # round real-valued energy
base="${arr[0]}${arr[1]}j_${arr[2]}TeV_${arr[3]}_${arr[4]}_${jalg}_${arr[5]}"

sql="
  FROM bh
  JOIN wt ON bh.id = wt.bh_id
  WHERE particle='${arr[0]}' and njets=${arr[1]} and
        energy=${arr[2]} and part='${arr[3]}' and
        wt.scales='${arr[4]}' and wt.pdf='${arr[5]}' and
        dset=${arr[6]}
"

min=`sqlite3 $db "SELECT min(sid) $sql"`
max=`sqlite3 $db "SELECT max(sid) $sql"`

(( begin = min     ))
(( end   = min - 1 ))

while true; do # loop over sid subgroup

(( begin = end + 1 ))
(( end = begin + group_size - 1 ))

if [ $((max - end)) -le $((group_size / 2)) ]; then
  (( end = max ))
fi

name="${base}_${begin}-${end}"

files="`sqlite3 $db "
  SELECT bh.dir, bh.file, wt.dir, wt.file $sql 
  and sid >= $begin and sid <= $end
" | sed 's/\(.*\)|\(.*\)|\(.*\)|\(.*\)/bh \1\/\2 w \3\/\4/'`"

if [ -z "$files" ]; then
  continue
fi

echo $name

# Form arguments

args="${arr[1]}j $jalg /home/ivanp/work/bh_analysis2/Hjets.bins out:$histdir/$name.root $files"

# Form temporary wrapper script
GCC_PATH=/msu/data/t3work9/ivanp/gcc-6.2.0

echo "#!/bin/bash" > $out/$name.sh
echo "LD_LIBRARY_PATH=$GCC_PATH/lib64:$GCC_PATH/lib:\$LD_LIBRARY_PATH" >> $out/$name.sh
echo 'echo $LD_LIBRARY_PATH' >> $out/$name.sh
echo "ldd $exe" >> $out/$name.sh
echo $exe $args >> $out/$name.sh
chmod +x $out/$name.sh

# Form condor script

echo "
Universe   = vanilla
Executable = $out/$name.sh
Output     = $out/$name.out
Error      = $out/$name.err
Log        = $out/$name.log
getenv = True
Queue
" | condor_submit - > /dev/null

if [ $end -ge $max ]; then break; fi

done

done
done

echo 'DONE!'

