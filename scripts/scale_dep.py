#!/usr/bin/env python

import sys, math

def read(fin):
  fac  = []
  ren  = []
  xsec = []
  unc  = []
  with open(fin,'r') as f:
    for line in [ line.split() for line in f ]:
      if len(line)==2:
        n_events = int(line[0])
      else:
        nums = [float(x) for x in line]
        fac.append(nums[0])
        ren.append(nums[1])
        xsec.append(nums[2]/n_events)
        unc.append(math.sqrt(nums[3])/n_events)
  return (fac,ren,xsec,unc)

# parts = [ read('scale_dep_H2j_%s.txt'%p) for p in ['B','R','I','V'] ]
parts = [ read(f) for f in sys.argv[1:] ]

# if not all(len(x) == len(parts[0]) for x in parts[1:]):
#   sys.stderr.write('Files don\'t have the same length')
#   sys.exit()

fac  = parts[0][0]
ren  = parts[0][1]
xsec = [ sum(y) for y in zip(*[ x[2] for x in parts ]) ]

for point in zip(fac,ren,xsec):
  # if point[2] > 0:
  # if point[1] < 0.15: continue
  # if point[1] > 0.35: continue
  # if point[0] > 0.3: continue
  # if point[1] < 0.3: continue
  print point[0], point[1], point[2]

