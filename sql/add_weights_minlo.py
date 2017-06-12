#!/usr/bin/env python

import sys, os, sqlite3, signal

if len(sys.argv)!=2:
    print 'usage:', sys.argv[0], 'ntuples.db'
    sys.exit(1)

path = '/msu/data/t3work8/ivanp/HT2'

db = sqlite3.connect(sys.argv[1])
cur = db.cursor()

def signal_handler(signal, frame): # exit on interupt
    print '\033[31;1m'+sys.argv[0]+' interupted\033[0m'
    db.commit() # save database
    sys.exit(1)
signal.signal(signal.SIGINT, signal_handler)

for f in os.listdir(path):
    if not f.endswith('.root'): continue
    print f

    cur.execute('SELECT id FROM ntuples WHERE' +\
                ' instr(info,\'ED\') and' +\
                ' file=\'{}\''.format(f))
    ids = cur.fetchall()

    if len(ids)==0: raise Exception('No matching ntuples')
    if len(ids)> 1: raise Exception('More than 1 matching ntuples')
    id1 = ids[0][0] # ids is a list of tuples
    print id1

    cur.execute( "insert into weights " + \
        "(ntuple_id,file,dir,scales,pdf,info) values " + \
        "({},'{}','{}','{}','{}','{}')".format(
            id1, f, path, 'minlo HT-hat-pp', 'CT10nlo', ''
        )
    )

db.commit()

