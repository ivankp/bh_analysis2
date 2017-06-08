#!/usr/bin/env python

import sys, os, sqlite3, subprocess, signal

if len(sys.argv)!=2:
    print 'usage:', sys.argv[0], 'ntuples.db'
    sys.exit(1)

path = '/msu/data/t3work4/luisonig'

db = sqlite3.connect(sys.argv[1])
cur = db.cursor()

# exit on interupt
def signal_handler(signal, frame):
    print '\033[31;1m'+sys.argv[0]+' interupted\033[0m'
    db.commit() # save database
    sys.exit(1)
signal.signal(signal.SIGINT, signal_handler)

dirs = []

for d1 in os.listdir(path):
    path_d1 = os.path.join(path,d1)
    if os.path.isfile(path_d1): continue
    for d2 in os.listdir(path_d1):
        if not 'NTuplesFiles' in d2: continue
        if '100' in d2: continue # skipp 100 TeV ntuples for now
        dirs.append(os.path.join(path_d1,d2))

for d in sorted(dirs):
    print '\033[34;1m'+d+'\033[0m'

    info = ''
    if '/ED' in d: info += 'ED '
    if '_mtop' in d: info += 'mtop '

    i = 0

    for f in os.listdir(d):
        if not f.endswith('.root'): continue
        t = f.split('_')

        print f

        rwgt = None
        if t[0].startswith('RWGT'):
            rwgt = t[0]
            t = t[1:]

        try:
            tree = subprocess.check_output(
                '/home/ivanp/work/bh_analysis2/bin/test_tree'+' '+d+'/'+f,
                shell=True)
        except Exception as e:
            print e
            continue
        tree = tree.split()

        cur.execute( "insert into ntuples " + \
            "(file,dir,particle,njets,energy,part,sid,info," + \
            "tree,nentries,nevents) values " + \
            "('{}','{}','{}',{},{},'{}','{}','{}','{}',{},{})".format(
                f, d, t[0][0], t[0][1], float(t[3])/500, t[2],
                (rwgt[4:]+'_' if rwgt else '') + t[-1].split('.')[0],
                info + ' '.join([t[1],t[4],t[5]]),
                tree[0], tree[1], tree[2]
            )
        )

    db.commit()

