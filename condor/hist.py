#!/usr/bin/env python

import sys, os, sqlite3, signal
from subprocess import Popen, PIPE
from common import *

def signal_handler(signal, frame): # exit on interupt
    print '\033[31;1m'+sys.argv[0]+' interrupted\033[0m'
    db.commit() # save database
    sys.exit(1)
signal.signal(signal.SIGINT, signal_handler)

if len(sys.argv) == 1:
    print 'usage: '+sys.argv[0]+' card.yml'
    sys.exit(0)

conf(sys.argv[1])

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
num_jobs = conf['num_jobs']
group_size = sum([ d.nevents for k, d in datasets.iteritems() ])/num_jobs

print "mean entries per job:", group_size

# groups = defaultdict(list)
# for k, d in datasets.iteritems():
#     d.group(group_size)
#     # print k
#     # i = 0
#     for g in d.gs:
#         # i += 1
#         # groups[k].append((i,k,g))
#         groups[k].append(g)
#         num_jobs_count += 1
#         # print '  {}: {}'.format(i,len(g))

# groups = sorted( [ (k,i,g) for k, gs in groups.iteritems()
#                            for i, g in enumerate(gs,1) ],
#                  key=lambda g: len(g[2]), reverse=True )

groups = sorted( [ (k,(i,),g) for k, d in datasets.iteritems()
                   for i, g in enumerate((d.group(group_size), d.gs)[1],1) ],
                 key=lambda g: len(g[2]), reverse=True )

# for g in groups:
#     print g[0], g[1], len(g[2])
# print ""

ng = len(groups)
while ng<num_jobs:
    i = 0
    while i<ng:
        g = groups[i]
        fs = g[2]
        nf = len(fs)
        # print g[:-1], nf
        if nf>1:
            groups = groups[:i] \
                   + [ (g[0],g[1]+(0,),fs[:nf/2]), (g[0],g[1]+(1,),fs[nf/2:]) ] \
                   + groups[i+1:]
            ng += 1
            if ng==num_jobs: break
            i += 1
        i += 1
# print ng
# print ""

groups_d = defaultdict(list)
for g in groups:
    l = groups_d[g[0]]
    l.append((len(l)+1,)+g[1:])

groups = sorted(
        [ (k,g[0],g[2]) for k, gs in groups_d.iteritems() for g in gs ],
        key=lambda g: g[1] )

# submit jobs =======================================================
i = 0
args_form = conf['args_form']
for g in groups:
    i += 1
    base = '{}_{}'.format(*g[:2])

    print '{:3}/{} : {} : {}'.format(i,ng,base,len(g[2]))

    script = conf['out']+'/'+base+'.sh'

    args_dict = {'base': base, 'njets': base[1]} # FIXME
    args = conf['args'].format(**dict(conf.yaml.items()+args_dict.items()))
    args += ' '+ ' '.join([ args_form.format(*f) for f in g[2] ])
    make_script(script, args)

    p = Popen((condor_sh,script))

