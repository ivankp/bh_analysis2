#!/usr/bin/env python

import sys, os, sqlite3
from subprocess import Popen, PIPE
from common import *

conf('cards/minlo.yml')

group_size = 100000
db = sqlite3.connect(conf['database'])

if not os.path.isdir(conf['out']): os.makedirs(conf['out'])

# sql query =========================================================
db.row_factory = sqlite3.Row
cur = db.cursor()
# FROM ntuples JOIN weights ON weights.ntuple_id = ntuples.id
cur.execute('''
SELECT particle,njets,part,nevents,dir,file
FROM ntuples
WHERE '''+conf['select'])
for f in cur.fetchall():
    dset = '{}{}j_{}'.format(*f)
    ntuple = ( f['dir'], f['file'] )
    datasets[dset].add( f['nevents'], ntuple )

# form groups of files for condor jobs ==============================
groups = []
for k, d in datasets.iteritems():
    d.group(group_size)
    # print k
    i = 0
    for g in d.gs:
        i += 1
        groups.append((i,k,g))
        # print '  {}: {}'.format(i,len(g))

# submit jobs =======================================================
njobs = len(groups)
i = 0
for g in sorted(groups, key=lambda g: g[0]):
    i += 1
    base = '{}_{}'.format(g[1], g[0])

    print '{:3}/{} : {}'.format(i,njobs,base)

    script = conf['out']+'/'+base+'.sh'
    args = ' '.join([ f[0]+'/'+f[1] for f in g[2] ])
    make_script(script, args)
        
    p = Popen((condor_sh,script))

