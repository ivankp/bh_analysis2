#!/usr/bin/env python

import subprocess, glob

path = '/home/ivanp/work/bh_analysis2'
out = path+'/out/scale_dep'

def condor(basename, args):
  return '''\
Universe   = vanilla
Executable = {2}/bin/scale_dep
Arguments  = {3}/{0}.dat {1}
Output     = {3}/{0}.out
Error      = {3}/{0}.err
Log        = {3}/{0}.log
getenv = True
Queue
'''.format(basename,args,path,out)

for pt in [30,50,100]:
  for Nj in [1,2,3]:
    for P in ['B','RS','I','V']:
      basename = 'H{0}j_13TeV_CT10nlo_antikt4_jetpt{2}_HThpp_{1}'.format(Nj,P,pt)
      print basename
      cf = out+'/'+basename+'.condor'
      f = open(cf,'w')
      f.write(condor(basename,
          "{0}j jetpt{1} antikt4 ".format(Nj,pt) + ' '.join(
          glob.glob('/msu/data/t3work4/luisonig/H{0}jets_ggf/NTuplesFiles/H{0}.0j_GGFHT_{1}_6500_pt25.0_eta4.5_r100_{2}.root'.format(
            Nj, P, '100' if (P!='V' or Nj==1) else ( '1[0-4][0-9]' if Nj==2 else '1[0-9][0-9]' ) )))))
      f.close()
      subprocess.Popen('condor_submit '+cf, shell=True)

