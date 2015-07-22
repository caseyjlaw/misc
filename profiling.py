#!/usr/bin/env python

import pstats, cProfile
import numpy as n
import pyximport
pyximport.install()
import leanpipe_cython as lib
#import leanpipe_external as lib
import leanpipe
#import qimg_cython as qimg

nint = 200
nch = 256
nbl = 27
npol = 2
nants = 27
#arr = n.ma.ones(nint*nbl*nch*npol, dtype=n.complex128).reshape(nint,nbl,nch,npol)
arr=0
l1 = 0.001
m1 = 0.001
u = n.arange(nint*nbl, dtype=n.float).reshape(nint,nbl)
v = n.arange(nint*nbl, dtype=n.float).reshape(nint,nbl)
d = {'freq' : n.arange(1400,1400+2*nch,2)/1000.,
'freq_orig' : n.arange(1400,1400+2*nch,2)/1000.,
'l0' : 0., 
'm0' : 0.,
'chans': range(nch),
'nchan': nch,
'inttime': 0.01,
#'ants' : n.array([ 1,  2,  3,  4,  5,  6,  7,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28]),
'nants' : nants,
'ants' : n.arange(1,nants+1),
'blarr' : n.array([[ 1,  2], [ 1,  3], [ 1,  4], [ 1,  5], [ 1,  6], [ 1,  7], [ 1,  9], [ 1, 10], [ 1, 11], [ 1, 12], [ 1, 13], [ 1, 14], [ 1, 15], [ 1, 16], [ 1, 17], [ 1, 18], [ 1, 19], [ 1, 20], [ 1, 21], [ 1, 22], [ 1, 23], [ 1, 24], [ 1, 25], [ 1, 26], [ 1, 27], [ 1, 28], [ 2,  3], [ 2,  4], [ 2,  5], [ 2,  6], [ 2,  7], [ 2,  9], [ 2, 10], [ 2, 11], [ 2, 12], [ 2, 13], [ 2, 14], [ 2, 15], [ 2, 16], [ 2, 17], [ 2, 18], [ 2, 19], [ 2, 20], [ 2, 21], [ 2, 22], [ 2, 23], [ 2, 24], [ 2, 25], [ 2, 26], [ 2, 27], [ 2, 28], [ 3,  4], [ 3,  5], [ 3,  6], [ 3,  7], [ 3,  9], [ 3, 10], [ 3, 11], [ 3, 12], [ 3, 13], [ 3, 14], [ 3, 15], [ 3, 16], [ 3, 17], [ 3, 18], [ 3, 19], [ 3, 20], [ 3, 21], [ 3, 22], [ 3, 23], [ 3, 24], [ 3, 25], [ 3, 26], [ 3, 27], [ 3, 28], [ 4,  5], [ 4,  6], [ 4,  7], [ 4,  9], [ 4, 10], [ 4, 11], [ 4, 12], [ 4, 13], [ 4, 14], [ 4, 15], [ 4, 16], [ 4, 17], [ 4, 18], [ 4, 19], [ 4, 20], [ 4, 21], [ 4, 22], [ 4, 23], [ 4, 24], [ 4, 25], [ 4, 26], [ 4, 27], [ 4, 28], [ 5,  6], [ 5,  7], [ 5,  9], [ 5, 10], [ 5, 11], [ 5, 12], [ 5, 13], [ 5, 14], [ 5, 15], [ 5, 16], [ 5, 17], [ 5, 18], [ 5, 19], [ 5, 20], [ 5, 21], [ 5, 22], [ 5, 23], [ 5, 24], [ 5, 25], [ 5, 26], [ 5, 27], [ 5, 28], [ 6,  7], [ 6,  9], [ 6, 10], [ 6, 11], [ 6, 12], [ 6, 13], [ 6, 14], [ 6, 15], [ 6, 16], [ 6, 17], [ 6, 18], [ 6, 19], [ 6, 20], [ 6, 21], [ 6, 22], [ 6, 23], [ 6, 24], [ 6, 25], [ 6, 26], [ 6, 27], [ 6, 28], [ 7,  9], [ 7, 10], [ 7, 11], [ 7, 12], [ 7, 13], [ 7, 14], [ 7, 15], [ 7, 16], [ 7, 17], [ 7, 18], [ 7, 19], [ 7, 20], [ 7, 21], [ 7, 22], [ 7, 23], [ 7, 24], [ 7, 25], [ 7, 26], [ 7, 27], [ 7, 28], [ 9, 10], [ 9, 11], [ 9, 12], [ 9, 13], [ 9, 14], [ 9, 15], [ 9, 16], [ 9, 17], [ 9, 18], [ 9, 19], [ 9, 20], [ 9, 21], [ 9, 22], [ 9, 23], [ 9, 24], [ 9, 25], [ 9, 26], [ 9, 27], [ 9, 28], [10, 11], [10, 12], [10, 13], [10, 14], [10, 15], [10, 16], [10, 17], [10, 18], [10, 19], [10, 20], [10, 21], [10, 22], [10, 23], [10, 24], [10, 25], [10, 26], [10, 27], [10, 28], [11, 12], [11, 13], [11, 14], [11, 15], [11, 16], [11, 17], [11, 18], [11, 19], [11, 20], [11, 21], [11, 22], [11, 23], [11, 24], [11, 25], [11, 26], [11, 27], [11, 28], [12, 13], [12, 14], [12, 15], [12, 16], [12, 17], [12, 18], [12, 19], [12, 20], [12, 21], [12, 22], [12, 23], [12, 24], [12, 25], [12, 26], [12, 27], [12, 28], [13, 14], [13, 15], [13, 16], [13, 17], [13, 18], [13, 19], [13, 20], [13, 21], [13, 22], [13, 23], [13, 24], [13, 25], [13, 26], [13, 27], [13, 28], [14, 15], [14, 16], [14, 17], [14, 18], [14, 19], [14, 20], [14, 21], [14, 22], [14, 23], [14, 24], [14, 25], [14, 26], [14, 27], [14, 28], [15, 16], [15, 17], [15, 18], [15, 19], [15, 20], [15, 21], [15, 22], [15, 23], [15, 24], [15, 25], [15, 26], [15, 27], [15, 28], [16, 17], [16, 18], [16, 19], [16, 20], [16, 21], [16, 22], [16, 23], [16, 24], [16, 25], [16, 26], [16, 27], [16, 28], [17, 18], [17, 19], [17, 20], [17, 21], [17, 22], [17, 23], [17, 24], [17, 25], [17, 26], [17, 27], [17, 28], [18, 19], [18, 20], [18, 21], [18, 22], [18, 23], [18, 24], [18, 25], [18, 26], [18, 27], [18, 28], [19, 20], [19, 21], [19, 22], [19, 23], [19, 24], [19, 25], [19, 26], [19, 27], [19, 28], [20, 21], [20, 22], [20, 23], [20, 24], [20, 25], [20, 26], [20, 27], [20, 28], [21, 22], [21, 23], [21, 24], [21, 25], [21, 26], [21, 27], [21, 28], [22, 23], [22, 24], [22, 25], [22, 26], [22, 27], [22, 28], [23, 24], [23, 25], [23, 26], [23, 27], [23, 28], [24, 25], [24, 26], [24, 27], [24, 28], [25, 26], [25, 27], [25, 28], [26, 27], [26, 28], [27, 28]])
#'blarr' : n.array([[ 1,  2], [ 1,  3], [ 1,  4], [ 1,  5], [ 1,  6], [ 1,  7], [ 1,  8], [ 1, 9], [ 2,  3], [ 2,  4], [ 2,  5], [ 2,  6], [ 2,  7], [ 2,  8], [ 2, 9], [3,  4], [ 3,  5], [ 3,  6], [ 3,  7], [ 3,  8], [ 3, 9], [ 4,  5], [ 4,  6], [ 4,  7], [ 4,  8], [ 4, 9], [ 5,  6], [ 5,  7], [ 5,  8], [ 5, 9], [6,  7], [ 6,  8], [ 6, 9], [7,  8], [ 7, 9], [ 8, 9]])
#'blarr' : n.array([[ 1,  2], [ 1,  3], [ 2,  3]])
}
d['datadelay'] = lib.calc_delay(d, 0)

#cProfile.runctx("[leanpipe.phaseshift(arr, d, l1, m1, u, v) for (l1,m1) in [(0.01,0.01),(-0.01,0.01),(0.1,0.1),(-0.01,0.01),(0.1,0.1),(-0.01,0.01),(0.1,0.1),(-0.01,0.01),(0.1,0.1),(-0.01,0.01)]]", globals(), locals(), "Profile.prof")
#cProfile.runctx("leanpipe.pipe_thread(filename='sysstartL.56526.02837960648.ms', telcalfile='sysstartL.56526.02837960648.GN', nints=10, iterint=10, nskip=10, scan=0, selectpol=['RR','LL'], datacol='data', size=51200, res=100, dmarr=range(0,1), fwhmfield=0.02, fwhmsurvey=0.02)", globals(), locals(), "Profile.prof")
#cProfile.runctx("d = leanpipe.pipe_thread(telcalfile='TRSR0001_sb24704013_1_000.56530.82972260417.GN', filename='TRSR0001_sb24704013_1_000.56530.82972260417.ms', datacol='data', size=25600, res=50, chans=range(128), dmarr=range(0,100,20), searchtype='imageallstat', nints=1000, iterint=200, nskip=100, scan=2, spw=[0,1], selectpol=['RR','LL'])", globals(), locals(), "Profile.prof")
cProfile.runctx("d = leanpipe.pipe_thread(telcalfile='TRSR0001_sb24704013_1_128ch.56531.76190895833.GN', filename='TRSR0001_few_scans.ms', datacol='data', size=25600, res=50, chans=range(5,123)+range(133,251), dmarr=range(0,3000,25), searchtype='imageallstat', nints=16000, iterint=200, nskip=0, scan=2, spw=[0,1], selectpol=['RR','LL'], sigma_image=10)", globals(), locals(), "Profile.prof")
#cProfile.runctx("d = leanpipe.pipe_thread(filename='12A-336_1b_j0628.ms', nints=3000, iterint=256, nskip=11000, scan=2, chans=range(5,59), spw=[0], datacol='corrected_data', selectpol=['I'], dmarr=range(0,330,33), sigma=7., fwhmfield=0.03, fwhmsurvey=0.03)", globals(), locals(), "Profile.prof")
#cProfile.runctx("[leanpipe.calc_delay(d, dm) for dm in range(100)]", globals(), locals(), "Profile.prof")
#cProfile.runctx("[lib.dedisperse(arr, d, dm, verbose=1) for dm in range(0,330,33)]", globals(), locals(), "Profile.prof")
#cProfile.runctx("triples = lib.make_triples(d); d['ntr'] = len(triples); bispectra = lib.make_bispectra(arr, triples)", globals(), locals(), "Profile.prof")
cProfile.runctx("triples = lib.make_triples(d); d['ntr'] = len(triples); lib.phaseshift(arr, d, l1, m1, u, v); lib.dedisperse(arr, d, 33); bispectra = lib.make_bispectra(arr, triples); bispectra[-d['datadelay'].max():] = n.ma.masked; leanpipe.detect_bispectra(bispectra, d)", globals(), locals(), "Profile.prof")
#cProfile.runctx("lib.dedisperse(arr, d, 33); qimg.imgallstat(u[0], v[0], arr.mean(axis=3).reshape(nint, nbl*nch), 25600,100)", globals(), locals(), "Profile.prof")
#cProfile.runctx("qimg.imgall(u[0], v[0], arr.mean(axis=3).reshape(nint, nbl*nch), 8*25600,100)", globals(), locals(), "Profile.prof")

def pipe(arr=arr, d=d):
    triples = lib.make_triples(d)
    d['ntr'] = len(triples)
    for dm in range(0,1000,100):
        lib.dedisperse(arr, d, dm)
#        arr = n.ma.masked_equal(arr, 0j, copy=False)
        bispectra = lib.make_bispectra(arr, triples)
        bispectra = n.ma.masked_array(bispectra, bispectra==0j)
        bispectra[-d['datadelay'].max():] = n.ma.masked
        leanpipe.detect_bispectra(bispectra, d)
    return d, triples, arr, bispectra

s = pstats.Stats("Profile.prof")
s.strip_dirs().sort_stats("time").print_stats()
