#!/usr/bin/env python

import sys, os, sqlite3, subprocess, signal

if len(sys.argv)!=2:
    print 'usage:', sys.argv[0], 'ntuples.db'
    sys.exit(1)

path = '/msu/data/t3work4/gosam_diphoton_jet/'

db = sqlite3.connect(sys.argv[1])
cur = db.cursor()

# exit on interupt
def signal_handler(signal, frame):
    print '\033[31;1m'+sys.argv[0]+' interupted\033[0m'
    db.commit() # save database
    sys.exit(1)
signal.signal(signal.SIGINT, signal_handler)

for p in [("born","B"),("real","RS"),("int","I"),("virt","V")]:
    d = path + p[0]
    print '\033[34;1m'+d+'\033[0m'

    info = 'AAj'

    i = 0

    for f in os.listdir(d):
        if not f.endswith('.root'): continue

        print f

        try:
            tree = subprocess.check_output(
                '/home/ivanp/work/bh_analysis2/bin/check_tree'+' '+d+'/'+f,
                shell=True)
        except Exception as e:
            print e
            continue
        tree = tree.split()

        cur.execute( "insert into ntuples " + \
            "(file,dir,particle,njets,energy,part,sid,info," + \
            "tree,nentries,nevents) values " + \
            "('{}','{}','{}',{},{},'{}','{}','{}','{}',{},{})".format(
                f, d, "AA", 1, 13., p[1], f[len(p[0]):-5], "Nico",
                tree[0], tree[1], tree[2]
            )
        )

    db.commit()

