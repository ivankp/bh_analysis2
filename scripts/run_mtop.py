#!/usr/bin/env python

import sys, os, sqlite3, signal
from subprocess import Popen

def until(l,total,pred):
    s = 0
    o = []
    for x in l:
        o.append(x)
        s += pred(x)
        if s >= total: break;
    return (s,o)

path = '/home/ivanp/work/bh_analysis2'
db  = sqlite3.connect(path+'/sql/ntuples.db')
cur = db.cursor()

select = 'SELECT dir,file,nevents FROM ntuples WHERE'

for njets in [1,2,3]:
    cur.execute(select +\
                ' instr(info,\'mtop\') and' +\
                ' energy=13 and njets={}'.format(njets))
    fin = cur.fetchall()

    cur.execute(select +\
                ' not instr(info,\'mtop\') and not instr(info,\'ED\') and' +\
                ' energy=13 and njets={}'.format(njets))
    inf = cur.fetchall()

    nfin = sum([ r[2] for r in fin ])
    ninf, inf = until(inf,max(nfin,50000000),lambda x: x[2])

    print 'Njets', njets
    print 'Nfin', nfin
    print 'Ninf', ninf

    name = path+'/mtop/inf_{}j'.format(njets)
    with open(name+'.out','w') as out, open(name+'.err','w') as err:
        p = Popen(
            [ path+'/bin/hist_Hjets_mtop',
              '{}j'.format(njets), 'out:'+name+'_raw.root' ] +\
            [ x[0]+'/'+x[1] for x in inf ],
            stdout=out, stderr=err)
    print name

    name = path+'/mtop/finite_{}j'.format(njets)
    with open(name+'.out','w') as out, open(name+'.err','w') as err:
        p = Popen(
            [ path+'/bin/hist_Hjets_mtop',
              '{}j'.format(njets), 'out:'+name+'_raw.root' ] +\
            [ x[0]+'/'+x[1] for x in fin ],
            stdout=out, stderr=err)
    print name

