import os, yaml
from collections import defaultdict

def parse_num_events(n_str):
    if len(n_str) > 1:
        s = n_str[-1]
        if s=='M': return int(n_str[:-1])*1000000;
        if s=='k': return int(n_str[:-1])*1000;
    return int(n_str)

class Conf:
    def __init__(self): pass
    def __call__(self, name):
        with open(name) as f:
            self.yaml = yaml.load(f)
        self['exe'] = os.path.join(self['path'],'bin',self['exe'])
        self['database'] = os.path.join(self['path'],self['database'])
        # self['out'] = os.path.join(self['path'],self['out'])
        if type(self['max_num_events']) is str:
            self['max_num_events'] = parse_num_events(self['max_num_events'])
    def __getitem__(self, key):
        return self.yaml[key]
    def __setitem__(self, key, value):
        self.yaml[key] = value

conf = Conf()

class collector:
    def __init__(self):
        self.files = []
        self.nevents = 0

    def add(self, ne, ntuple): # ne is number of events
        if conf['max_num_events']==0 or self.nevents < conf['max_num_events']:
            self.nevents += ne
            self.files.append((ntuple,ne))

    def group(self, ng, pred=None): # ng is group size
        self.gs = [] # groups
        if pred!=None:
            self.files = sorted(self.files, key=lambda f: pred(f[0]))
        g = []
        n = 0
        for f in self.files:
            n += f[1]
            g.append(f[0])
            # g.append(f)
            if n >= ng:
                self.gs.append(g)
                g = []
                n = 0
        if len(g)>0: self.gs.append(g)

    def __str__(self):
        return '{}:{}'.format(self.files,self.nevents)

datasets = defaultdict(collector)

condor_sh = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),'condor.sh')

def make_script(name, args):
    with open(name,'w') as f:
        f.write('''#!/bin/bash

exe={0}

if [ ! -x \"$exe\" ]; then
  >&2 echo "$(uname -n) does not see $exe"
  exit
fi

ldd $exe

$exe {1}
'''.format(conf['exe'],args))
    os.chmod(name,0o775)

