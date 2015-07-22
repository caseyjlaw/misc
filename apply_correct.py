#!/opt/local/Library/Frameworks/Python.framework/Versions/2.6/Resources/Python.app/Contents/MacOS/Python
"""
Apply calibration parameters to a data set.
"""

import aipy as a
import optparse, sys, os

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, src=True)
o.add_option('--no_gain', action='store_true', help='Do not change the gain (i.e. no flux cal).')
o.add_option('--no_phs', action='store_true', help='Do not change the phase (i.e. no phase cal).')
opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
if not opts.src is None:
    if not opts.src.startswith('zen'):
        srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
        src = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs).values()[0]
    else: src = 'z'
else: src = None
del(uv)

import numpy as n

def mfunc(uv, p, d, f):
    uvw,t,(i,j) = p
    pol = a.miriad.pol2str[uv['pol']]
    aa.set_active_pol(pol)
    aa.set_jultime(t)
    if not src is None and not type(src) == str: src.compute(aa)
    if not opts.no_gain: d /= aa.passband(i, j)
    if not opts.no_phs: d = aa.phs2src(d, src, i, j)
    uvw = aa.get_baseline(i, j, src)
    p = uvw, t, (i, j)
    return p, d, f

for infile in args:
    outfile = infile+'C'
    print infile,'-->',outfile
    if os.path.exists(outfile):
        print 'File exists, skipping....'
        continue 

    uvi = a.miriad.UV(infile)
    uvo = a.miriad.UV(outfile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc, raw=True, append2hist='APPLY_CAL: ' + ' '.join(sys.argv) + '\n')
