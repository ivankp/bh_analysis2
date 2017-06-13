#!/usr/bin/env python

import sys, sqlite3, signal

db = sqlite3.connect('ntuples.db')
cur = db.cursor()

def signal_handler(signal, frame): # exit on interupt
    print '\033[31;1m'+sys.argv[0]+' interupted\033[0m'
    sys.exit(1)
signal.signal(signal.SIGINT, signal_handler)

cur.execute('SELECT id,file,particle,njets FROM ntuples WHERE' +\
            ' not instr(info,\'ED\') and' +\
            ' not instr(info,\'mtop\') and' +\
            ' not instr(info,\'GGFHTLS\') and' +\
            ' particle=\'H\' and energy=13')
ntuples = cur.fetchall() # list of tuples

out = '/msu/data/t3work8/ivanp/weights'
scale = 'HT1'
pdf = 'PDF4LHC15_nlo_30'

for x in ntuples:
    print x[1]
    cur.execute( "insert into weights " + \
        "(ntuple_id,file,dir,scales,pdf,info) values " + \
        "({},'{}','{}','{}','{}','{}')".format(
            x[0], x[1], out+'/'+x[2]+str(x[3])+'j/'+scale+'-'+pdf,
            scale, pdf, 'hgam 2017'
        )
    )

db.commit()

