#!/usr/bin/env python

import sys, os, sqlite3, signal
from subprocess import Popen, PIPE
from common import *

def signal_handler(signal, frame): # exit on interupt
    print '\033[31;1m'+sys.argv[0]+' interupted\033[0m'
    db.commit() # save database
    sys.exit(1)
signal.signal(signal.SIGINT, signal_handler)

if len(sys.argv) == 1:
    print 'usage: '+sys.argv[0]+' card.yml'
    sys.exit(0)

conf(sys.argv[1])

group_size = 10000000
db = sqlite3.connect(conf['database'])

if not os.path.isdir(conf['out']): os.makedirs(conf['out'])

# sql query =========================================================
# db.row_factory = sqlite3.Row
cur = db.cursor()
cur.execute(conf['query'])
dset_form = conf['dset_form']
n_dset_form_fields = dset_form.count('{}')
for f in cur.fetchall():
    dset = dset_form.format(*f)
    ntuple = f[n_dset_form_fields+1:]
    f1 = os.path.join(*ntuple[:2])
    if not os.path.isfile(f1):
        print '\033[33mWarning\033[0m: no file {}'.format(f1)
        continue
    f2 = os.path.join(*ntuple[2:])
    if not os.path.isfile(f2):
        print '\033[33mWarning\033[0m: no file {}'.format(f2)
        continue
    datasets[dset].add( f[n_dset_form_fields], ntuple )

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
args_form = conf['args_form']
for g in sorted(groups, key=lambda g: g[0]):
    i += 1
    base = '{}_{}'.format(g[1], g[0])

    print '{:3}/{} : {}'.format(i,njobs,base)

    script = conf['out']+'/'+base+'.sh'

    args_dict = {'base': base, 'njets': base[1]} # FIXME
    args = conf['args'].format(**dict(conf.yaml.items()+args_dict.items()))
    args += ' '+ ' '.join([ args_form.format(*f) for f in g[2] ])
    make_script(script, args)

    p = Popen((condor_sh,script))

