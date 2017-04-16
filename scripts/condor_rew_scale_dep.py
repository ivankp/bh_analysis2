#!/usr/bin/env python

import subprocess, glob

path = '/home/ivanp/work/bh_analysis2'
out = path+'/out/scale_dep'

def condor(basename, ntuples):
  return '''\
Universe   = vanilla
Executable = {2}/bin/scale_dep
Arguments  = {3}/{0}.dat {1}
Output     = {3}/{0}.out
Error      = {3}/{0}.err
Log        = {3}/{0}.log
getenv = True
Queue
'''.format(basename, ' '.join(ntuples),path,out)

# print condor('scale_dep_H2j_V',['~/disk2/links/H2j/H2.0j_GGFHT_V_6500_pt25.0_eta4.5_r100_1{1..5}*.root'])

for P in ['B','RS','I']:
  basename = 'scale_dep_H2j_%s' % P
  print basename
  cf = out+basename+'.condor'
  f = open(cf,'w')
  f.write(condor(basename,
      glob.glob('/msu/data/t3work4/luisonig/H2jets_ggf/NTuplesFiles/H2.0j_GGFHT_%s_6500_pt25.0_eta4.5_r100_%s.root' % (
        P, '100' if P!='V' else '1[0-4][0-9]' ))))
  f.close()
  # print cmd
  subprocess.Popen('condor_submit '+cf, shell=True)

# for P in ['B', 'RS', 'I', 'V']:

# echo "
# " | condor_submit - > /dev/null

