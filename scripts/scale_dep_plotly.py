#!/usr/bin/env python

import os, sys

try:
    os.chdir(sys.argv[1])
except Exception as e:
    sys.stderr.write(str(e)+'\n')
    sys.stderr.write('cannot cd into '+sys.argv[1]+'\n')
    sys.exit(1)

groups = [ f.rsplit('_',1)[0] for f in os.listdir('.') if f.endswith('_NLO.dat') ]

done_name = 'done.lst'
done = [ l.split() for l in open(done_name,'r') ] \
       if os.path.isfile(done_name) else []

need = set(groups).difference([ d[0] for d in done ])
# print need

if len(need)==0:
    print "Nothing to be done"
    sys.exit(0)

import math
import plotly.plotly as py
import plotly.graph_objs as go

for name in need:
    print name
    data = [ [ float(x) for x in line.split() ] for line in open(name+'_NLO.dat','r') ]
    x = [ math.log(p[0],2) for p in data if p[2]>0 ]
    y = [ math.log(p[1],2) for p in data if p[2]>0 ]
    # x = [ p[0] for p in data ]
    # y = [ p[1] for p in data ]
    z = [ p[2] for p in data if p[2]>0 ]

    trace1 = go.Scatter3d(
        x=x, y=y, z=z,
        mode='markers',
        marker=dict(
            size=8,
            color=z,             # set color to an array/list of desired values
            colorscale='Viridis' # choose a colorscale
            # opacity=0.8
        )
    )

    data = [trace1]
    layout = go.Layout(
        title=name,
        scene=go.Scene(
            xaxis=go.XAxis(title='log2 fac [HT\'\']'),#type='log'),
            yaxis=go.YAxis(title='log2 ren [HT\'\']'),#type='log'),
            zaxis=go.ZAxis(title='xsec [pb]')
        )
    )
    fig = go.Figure(data=data, layout=layout)
    url = py.plot(fig, filename='scale_dep_'+name, auto_open=False)
    print url

    done.append([name,url])

with open(done_name,'w') as f:
    for line in done:
        f.write(' '.join(line)+'\n')

