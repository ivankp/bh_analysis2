#!/usr/bin/env python

import sys, sqlite3
from subprocess import Popen
from collections import defaultdict

exe = 'hist_Hjets'
# N = 10000000
N = (0 if len(sys.argv)==1
    else int(sys.argv[1]) if not sys.argv[1][-1].isalpha()
    else int(sys.argv[1][:-1])
         * (1000000 if sys.argv[1][-1]=='M' else 1000))

path = '/home/ivanp/work/bh_analysis2'
db  = sqlite3.connect(path+'/sql/ntuples.db')
cur = db.cursor()

class collector:
    def __init__(self):
        self.files = []
        self.nevents = 0
    def add(self,n,ntuple,weights):
        if N==0 or self.nevents < N:
            self.nevents += n
            self.files.append((ntuple,weights))
    def __str__(self):
        return '{}:{}'.format(self.files,self.nevents)

sets = defaultdict(collector)
cur.execute('''
SELECT nevents,njets,part,ntuples.dir,ntuples.file,weights.dir,weights.file
FROM weights JOIN ntuples ON weights.ntuple_id = ntuples.id
WHERE energy=13 and weights.scales=\'minlo HT-hat\'\'\' and pdf=\'CT10nlo\'
and njets=2
''')
for f in cur.fetchall():
    sets[f[1:3]].add( f[0], f[3]+'/'+f[4], f[5]+'/'+f[6] )

for k in sets:
    # print '{}: {}'.format(k,sets[k])
    name = path+'/out/default_{}j_{}'.format(*k)
    with open(name+'.out','w') as out, open(name+'.err','w') as err:
        p = Popen(
            [ path+'/bin/'+exe,
              '{}j'.format(k[0]), 'out:'+name+'_raw.root' ] +\
            # [ arg for f in sets[k].files for arg in ['bh', f[0], 'w', f[1]] ],
            [ f[0] for f in sets[k].files ],
            stdout=out, stderr=err)
    print name

