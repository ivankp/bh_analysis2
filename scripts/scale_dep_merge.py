#!/usr/bin/env python

import os, sys, re
from collections import defaultdict

def read(fin):
    fac  = []
    ren  = []
    xsec = []
    # unc  = []
    with open(fin,'r') as f:
        for line in [ line.split() for line in f ]:
            if len(line)==2:
                n_events = int(line[0])
            else:
                nums = [float(x) for x in line]
                fac.append(nums[0])
                ren.append(nums[1])
                xsec.append(nums[2]/n_events)
                # unc.append(math.sqrt(nums[3])/n_events)
    return (fac,ren,xsec)

try:
    os.chdir(sys.argv[1])
except Exception as e:
    sys.stderr.write(str(e)+'\n')
    sys.stderr.write('cannot cd into '+sys.argv[1]+'\n')
    sys.exit(1)

groups = defaultdict(list)
for g in [
    f.rsplit('_',1) for f in os.listdir('.') if re.search(r'(B|RS|I|V)\.dat$',f)
]:
    groups[g[0]].append(g[1])

for g in groups:
    parts = [ read(g+'_'+suf) for suf in groups[g] ]

    fac  = parts[0][0]
    ren  = parts[0][1]
    xsec = [ sum(y) for y in zip(*[ x[2] for x in parts ]) ]

    with open(g+'_NLO.dat','w') as f:
        print f.name
        for point in zip(fac,ren,xsec):
            f.write('%14.8e %14.8e %15.8e\n' % (point[0],point[1],point[2]))
            # f.write(' '.join(['%14.8e'%p for p in point])+'\n')

