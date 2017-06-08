#!/usr/bin/env python

import sys, sqlite3, locale, os, subprocess
from subprocess import Popen, PIPE
from collections import defaultdict

# N = 10000000
N = (0 if len(sys.argv)==1
    else int(sys.argv[1]) if not sys.argv[1][-1].isalpha()
    else int(sys.argv[1][:-1])
         * (1000000 if sys.argv[1][-1]=='M' else 1000))

locale.setlocale(locale.LC_ALL, 'en_US')
print 'N = ' + locale.format('%d', N, grouping=True) + '\n'

exe = 'hist_Hjets'
path = '/home/ivanp/work/bh_analysis2'
card = path+'/minlo/Run2.dat'
out = '/msu/data/t3work8/ivanp/minlo4'
db  = sqlite3.connect(path+'/sql/ntuples.db')
cur = db.cursor()

class collector:
    def __init__(self):
        self.files = []
        self.nevents = 0
    def add(self,n,ntuple):
        if N==0 or self.nevents < N:
            self.nevents += n
            self.files.append(ntuple)
    def __str__(self):
        return '{}:{}'.format(self.files,self.nevents)

sets = defaultdict(collector)
cur.execute('''
SELECT nevents,njets,part,dir,file,sid
FROM ntuples
WHERE energy=13 and particle=\'H\' and instr(info,\'ED\')
''')
for f in cur.fetchall():
    sets[f[1:3]].add( f[0], (f[3],f[4],f[5]) )

def condor(basename, args):
    sh = '{}/{}.sh'.format(out,basename)
    with open(sh,'w') as f:
        f.write('#!/bin/bash\n')
        f.write('cd {}\n'.format(out))
        f.write('/msu/data/t3work5/ivanp/4/sherpa/bin/Sherpa '\
                '-f {} {}\n'.format(card,args))
    os.chmod(sh,0o775)
    return '''\
Universe   = vanilla
Executable = {2}
Output     = {1}/{0}.out
Error      = {1}/{0}.err
Log        = {1}/{0}.log
getenv = True
Queue
'''.format(basename,out,sh)

for f in sorted(
    [ k+f for k in sets for f in sets[k].files ],
    key = lambda tup: (tup[4],tup[0],tup[1])
):
    print '{}j {:<2}: {}'.format(f[0],f[1],f[3])

    base = f[3][:-5]
    if os.path.isfile(out+'/'+base+'.wts.gz'):
        print 'already exists'
        continue

    p = Popen(('condor_submit','-'), stdin=PIPE, stdout=PIPE)
    p.stdin.write( condor(base,
            '\\$OUT_PATH:={} '.format(out) +\
            '\\$NJETS:={} \\$NLOPART:={} \\$SID:={}'.format(f[0],f[1],f[4])))
    p.communicate()
    p.stdin.close()

