#!/usr/bin/env python

import sys, os, sqlite3
from subprocess import Popen, PIPE
from collections import defaultdict

out = '/home/ivanp/work/bh_analysis2/hgam'
exe = 'hist_hgam'
chunk_size = 20

path = '/home/ivanp/work/bh_analysis2'
db  = sqlite3.connect(path+'/sql/ntuples.db')
cur = db.cursor()

def mkdir(d):
    if not os.path.isdir(d):
        os.makedirs(d)
for d in ['/condor','/raw']:
    mkdir(out+d)

sets = defaultdict(list)
cur.execute('''
SELECT particle,njets,part,
       ntuples.dir,ntuples.file,
       weights.dir,weights.file
FROM weights JOIN ntuples ON weights.ntuple_id = ntuples.id
WHERE weights.info=\'hgam 2017\' and njets=1
''')
for f in cur.fetchall():
    sets[f[0:3]].append((f[3]+'/'+f[4],f[5]+'/'+f[6]))

def chunks(l, n):
    for i in xrange(0, len(l), n):
        yield l[i:i + n]

def condor(s,i,g):
    name = '{}{}j{}'.format(*s[:3]) + '_' + str(i)
    print name
    base = out + '/condor/' + name
    with open(base+'.sh','w') as f:
        f.write('#!/bin/bash\n')
        f.write('\n')
        f.write('export LD_LIBRARY_PATH=/msu/data/t3work3/ivanp/gcc-7.2.0/hep/root-6.10.02/lib:/msu/data/t3work3/ivanp/gcc-7.2.0/hep/lib:/msu/data/t3work3/ivanp/gcc-7.2.0/gcc/lib64:/msu/data/t3work3/ivanp/gcc-7.2.0/gcc/lib:/msu/data/t3work3/ivanp/gcc-7.2.0/lib64:/msu/data/t3work3/ivanp/gcc-7.2.0/lib')
        f.write('\n')
        f.write(path+'/bin/' + exe + ' {}j '.format(s[1]) +\
                # ' --no-photon-cuts' +\
                ' \\\nout:' + out + '/raw/' + name + '.root \\\n' +\
                ' \\\n'.join([ 'bh {} w {}'.format(*ff) for ff in g ]))
    os.chmod(base+'.sh',0o775)
    return '''\
Universe   = vanilla
Executable = {0}.sh
Output     = {0}.out
Error      = {0}.err
Log        = {0}.log
+IsLongJob = true
getenv = True
Queue
'''.format(base)

for s in sorted(
    [ k+(chunks(f,chunk_size),) for k, f in sets.iteritems() ],
    key = lambda tup: tup[:3]
):
    i = 1
    for g in s[-1]:
        # condor(s,i,g)
        p = Popen(('condor_submit','-'), stdin=PIPE, stdout=PIPE)
        p.stdin.write( condor(s,i,g) )
        p.communicate()
        p.stdin.close()

        i += 1


