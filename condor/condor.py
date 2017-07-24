#!/usr/bin/env python

import sys, os, sqlite3
from subprocess import Popen, PIPE
from common import *

try:
    if len(sys.argv)>1: getN(sys.argv[1])
except ValueError:
    print "Invalid argument for number of events:", sys.argv[1]
    sys.exit(1)

exe = 'hist_'
group_size = 100000
path = '/home/ivanp/work/bh_analysis2'
# out = '/msu/data/t3work8/ivanp/scale_dep/HT2'
out = 'test'
db  = sqlite3.connect(path+'/sql/ntuples.db')

if not os.path.isdir(out): os.makedirs(out)

# sql query =========================================================
db.row_factory = sqlite3.Row
cur = db.cursor()
cur.execute('''
SELECT nevents,particle,njets,part,dir,file,id
FROM ntuples
WHERE energy=13 and particle=\'H\' and instr(info,\'ED\')
and part=\'V\'
''')
for f in cur.fetchall():
    dset = f['particle']+f['njets']+'j_'+f['part']
    ntuple = ( f['dir'], f['file'], f['id'] )
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

    script = out+'/'+base+'.sh'
    args = ' '.join([ f[0]+'/'+f[1] for f in g[2] ])
    make_script(script, exe, args)
        
    p = Popen((condor_sh(),script))

# for f in sorted(
#     [ k+f for k,fs in datasets.iteritems() for f in fs.files ],
#     key = lambda tup: (tup[4],tup[0],tup[1])
# ):
#     print '{}j {:<2}: {}'.format(f[0],f[1],f[3])
#     print f
#
#     base = f[3][:-5]
#     if os.path.isfile(out+'/'+f[3]):
#         print 'already exists'
#         continue

    # make_script(out+'/'+base+'.sh',exe,'')

    # p = Popen(('condor_submit','-'), stdin=PIPE, stdout=PIPE)
    # p.stdin.write( condor(base,
    #         '\\$OUT_PATH:={} '.format(out) +\
    #         '\\$NJETS:={} \\$NLOPART:={} \\$SID:={}'.format(f[0],f[1],f[4])))
    # p.communicate()
    # p.stdin.close()

