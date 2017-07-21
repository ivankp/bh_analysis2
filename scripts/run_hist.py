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

N = 30000000

locale.setlocale(locale.LC_ALL, '')
print 'N = ' + locale.format('%d', N, grouping=True)

exe = 'hist_Hjets_yy'
path = '/home/ivanp/work/bh_analysis2'
db  = sqlite3.connect(path+'/sql/ntuples.db')
cur = db.cursor()

for njets in [2]:
    for part in ['B','RS','I','V']:
        # cur.execute('SELECT dir,file,nevents FROM ntuples WHERE'\
        #             ' not instr(info,\'mtop\') and'\
        #             ' not instr(info,\'ED\') and'\
        #             ' not instr(info,\'GGFHTLS\') and'\
        #             ' particle=\'H\' and energy=13' +\
        #             ' and njets={} and part=\'{}\''.format(njets,part))
        cur.execute('SELECT ntuples.dir,ntuples.file,nevents,'\
                    'weights.dir,weights.file FROM ntuples'\
                    ' JOIN weights ON weights.ntuple_id = ntuples.id'\
                    ' WHERE'\
                    ' weights.scales=\'HT1\''\
                    ' and particle=\'H\' and energy=13' +\
                    ' and njets={} and part=\'{}\''.format(njets,part))
        ntuples = cur.fetchall()

        ntuples = until(ntuples, N, lambda x: x[2])[1]

        name = path+'/out/gionata_{}j_{}'.format(njets,part)
        with open(name+'.out','w') as out, open(name+'.err','w') as err:
            cmd = \
                [ path+'/bin/'+exe, 'config/gionata.bins',
                  '{}j'.format(njets), 'out:'+name+'_raw.root' ] +\
                [ arg for x in ntuples
                      for arg in ['bh',x[0]+'/'+x[1],'w',x[3]+'/'+x[4]] ]
            p = Popen(cmd, stdout=out, stderr=err)
        print name

