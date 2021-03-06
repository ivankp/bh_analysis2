#!/usr/bin/env python

import sys, sqlite3, locale
from subprocess import Popen

def until(l,total,pred):
    s = 0
    o = []
    for x in l:
        o.append(x)
        s += pred(x)
        if s >= total: break;
    return (s,o)

exe = 'hist_hgam'
N = 30000000

locale.setlocale(locale.LC_ALL,'')
print 'N = ' + locale.format('%d', N, grouping=True)

path = '/home/ivanp/work/bh_analysis2'
db  = sqlite3.connect(path+'/sql/ntuples.db')
cur = db.cursor()

for njets in [1,2]:
    for part in ['B','RS','I','V']:
    # for part in ['RS']:
        cur.execute('SELECT dir,file,nevents FROM ntuples WHERE'\
                    ' not instr(info,\'mtop\') and'\
                    ' not instr(info,\'ED\') and'\
                    ' not instr(info,\'GGFHTLS\') and'\
                    ' particle=\'H\' and energy=13' +\
                    ' and njets={} and part=\'{}\''.format(njets,part))
        ntuples = cur.fetchall()

        ntuples = until(ntuples, N, lambda x: x[2])[1]

        name = path+'/out/hgam_{}j_{}'.format(njets,part)
        with open(name+'.out','w') as out, open(name+'.err','w') as err:
            p = Popen(
                [ path+'/bin/'+exe,
                  '{}j'.format(njets), 'out:'+name+'_raw.root' ] +\
                sys.argv[1:] +\
                [ x[0]+'/'+x[1] for x in ntuples ],
                stdout=out, stderr=err)
        print name

