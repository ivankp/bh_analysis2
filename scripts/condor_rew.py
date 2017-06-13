#!/usr/bin/env python

import sys, os, sqlite3
from subprocess import Popen, PIPE
from collections import defaultdict

path = '/home/ivanp/work/bh_analysis2'
db  = sqlite3.connect(path+'/sql/ntuples.db')
cur = db.cursor()

sets = defaultdict(list)
cur.execute('''
SELECT particle,njets,part,scales,pdf,weights.dir,ntuples.dir,ntuples.file
FROM weights JOIN ntuples ON weights.ntuple_id = ntuples.id
WHERE weights.info=\'hgam 2017\'
''')
for f in cur.fetchall():
    if not os.path.isfile(f[5]+'/'+f[7]):
        sets[f[0:6]].append(f[6]+'/'+f[7])

def chunks(l, n):
    for i in xrange(0, len(l), n):
        yield l[i:i + n]

def condor(s,i,ff):
    base = '{}{}j{}'.format(*s[:3]) + '_' + str(i)
    print base
    out = s[5]
    if not os.path.isdir(out + '/condor/'):
        os.makedirs(out + '/condor/')
    base = out + '/condor/' + base
    with open(base+'.sh','w') as f:
        f.write('#!/bin/bash\n')
        f.write('cd {}\n'.format(out))
        f.write(path+'/bin/reweigh --pdf-variations' +\
                ' -o '+out + ' --scale='+s[3] + ' --pdf='+s[4] +\
                ' -i ' + ' '.join(ff) + '\n')
    os.chmod(base+'.sh',0o775)
    return '''\
Universe   = vanilla
Executable = {0}.sh
Output     = {0}.out
Error      = {0}.err
Log        = {0}.log
getenv = True
Queue
'''.format(base)

for s in sorted(
    [ k+(chunks(f,20),) for k, f in sets.iteritems() ],
    key = lambda tup: tup[0:3]
):
    i = 1
    for ff in s[-1]:
        condor(s,i,ff)
        p = Popen(('condor_submit','-'), stdin=PIPE, stdout=PIPE)
        p.stdin.write( condor(s,i,ff) )
        p.communicate()
        p.stdin.close()

        i += 1

