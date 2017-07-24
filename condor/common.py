import os
from collections import defaultdict

# read enough files to precess at least N events
# 0 means process all
N = 0

def getN_impl(n_str):
    if len(n_str) > 1:
        s = n_str[-1]
        if s=='M': return int(n_str[:-1])*1000000;
        if s=='k': return int(n_str[:-1])*1000;
    return int(n_str)

def getN(n_str):
    global N
    N = getN_impl(n_str)
    print 'N = {}'.format(N)

class collector:
    def __init__(self):
        self.files = []
        self.nevents = 0

    def add(self, ne, ntuple): # ne is number of events
        if N==0 or self.nevents < N:
            self.nevents += ne
            self.files.append((ntuple,ne))

    def group(self, ng, key=None): # ng is group size
        self.gs = [] # groups
        if key!=None:
            self.files = sorted(self.files, key=lambda f: key(f[0]))
        g = []
        n = 0
        for f in self.files:
            n += f[1]
            g.append(f[0])
            if n >= ng:
                self.gs.append(g)
                g = []
                n = 0
        if len(g)>0: self.gs.append(g)

    def __str__(self):
        return '{}:{}'.format(self.files,self.nevents)

datasets = defaultdict(collector)

condor_sh_str = os.path.dirname(os.path.realpath(__file__))+'/condor.sh'
def condor_sh():
    return condor_sh_str

def make_script(name, exe, args, recursive=True):
    with open(name,'w') as f:
        f.write('#!/bin/bash\n')
        if recursive:
            f.write('''
if [ ! -x \"{}\" ]; then
  {} {}
  exit
fi
'''.format(exe,condor_sh_str,name))
        f.write('ldd '+exe+'\n\n')
        f.write(exe+' '+args+'\n')
    os.chmod(name,0o775)


