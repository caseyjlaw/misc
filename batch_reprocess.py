#!/usr/bin/env python

import leanpipet, qrtt
import sys, os, string
import numpy as n

arg0 = sys.argv.index('batch_reprocess.py')
msfile = sys.argv[arg0+1]
if len(sys.argv[arg0+2:]):
    ndm = int(sys.argv[arg0+2])
else:
    ndm = 1
if len(sys.argv[arg0+3:]):
    npix = int(sys.argv[arg0+3])
else:
    npix = 0
if len(sys.argv[arg0+4:]):
    res = int(sys.argv[arg0+4])
else:
    res = 0

threshold = 6.5         # 6.5sigma should produce 7.5 false+ for 512x512 imgs, 30 dms, and 23900 ints (i.e., a typical run on a node)
filtershape = 'b'
chans = range(6,122)+range(134,250)       # 128 ch/spw version cutting 5% of channels from each edge plus bad chans at bottom
flagmode = 'ring2medcht1.5badbp2blstd3'
dmarrall = [0,19.2033,38.4033,57.6025,76.8036,96.0093,115.222,134.445,153.68,172.93,192.198,211.486,230.797,250.133,269.498,288.894,308.323,327.788,347.292,366.837,386.426,406.062,425.747,445.484,465.276,485.125,505.033,525.005,545.042,565.147,585.322,605.571,625.896,646.3,666.786,687.355,708.012,728.759,749.598,770.532,791.565,812.699,833.936,855.28,876.733,898.299,919.979,941.778,963.697,985.741,1007.91,1030.21,1052.64,1075.21,1097.92,1120.76,1143.76,1166.9,1190.19,1213.63,1237.23,1260.99,1284.92,1309.01,1333.27,1357.7,1382.31,1407.09,1432.06,1457.22,1482.56,1508.1,1533.83,1559.76,1585.89,1612.23,1638.77,1665.53,1692.51,1719.7,1747.12,1774.77,1802.64,1830.75,1859.1,1887.69,1916.53,1945.61,1974.95,2004.54,2034.39,2064.51,2094.9,2125.56,2156.49,2187.71,2219.21,2250.99,2283.07,2315.45,2348.13,2381.12,2414.41,2448.02,2481.94,2516.19,2550.77,2585.68,2620.92,2656.51,2692.44,2728.72,2765.36,2802.36,2839.72,2877.45,2915.55,2954.04,2992.91]
fileroot = string.join(msfile.split('_')[:2], '_')
gainfile = fileroot + '.g2'
bpfile = fileroot + '.b1'
telcalfile = '13B-409_sb25785146_2.56560.079294421295.GN'

# size of data in ints
datadelay = [0,   3,   6,   9,  12,  15,  18,  21,  23,  26,  29,  32,  35, 38,  41,  44,  47,  50,  53,  56,  59,  62,  65,  68,  71,  74, 77,  80,  83,  86,  89,  92,  95,  99, 102, 105, 108, 111, 114, 118, 121, 124, 127, 130, 134, 137, 140, 144, 147, 150, 154, 157, 161, 164, 167, 171, 174, 178, 182, 185, 189, 192, 196, 200, 203, 207, 211, 215, 218, 222, 226, 230, 234, 238, 242, 246, 250, 254, 258, 262, 266, 271, 275, 279, 284, 288, 292, 296, 300, 306, 310, 314, 320, 324, 328, 334, 338, 342, 348, 352, 358, 362, 368, 372, 378, 384, 388, 394, 400, 404, 410, 416, 422, 426, 432, 438, 444, 450, 456]
canddm = int(msfile.split('.')[-2].split('-')[-1])
nints = datadelay[canddm] + 200   # +- 100 added to cand dm track to make it work robustly
maxiter = 200
iterint = nints/(((nints-n.mod(nints,maxiter))/maxiter)+2)    # keeps iterint under maxiter
#nints = 550    # max snippet size
#blocks = 2
#iterint = nints/blocks

dmc = dmarrall[canddm]
if canddm > 1:
    ddm = dmarrall[canddm] - dmarrall[canddm-1]
    dmarr = dmc + n.array([n.sum([2*(ddm/ndm)*i*(-1)**i for i in range(ndm)][:j]) for j in range(1,ndm+1)])   # yeah, but it works!
elif canddm == 1:
    ddm = dmarrall[canddm+1] - dmarrall[canddm]
    dmarr = range(dmarrall[canddm]-ddm, dmarrall[canddm]+ddm, 2*ddm/ndm)
else:
    ddm = dmarrall[canddm+1] - dmarrall[canddm]
    dmarr = range(0, 2*ddm, 2*ddm/ndm)

print 'Searching %s for dmrange %s with %d ints by %d ints.' % (msfile, dmarr, nints, iterint)
d = leanpipet.pipe_thread(filename=msfile, nints=nints, nskip=0, iterint=iterint, spw=[0,1], chans=chans, dmarr=dmarr, fwhmsurvey=0.5, fwhmfield=0.5, selectpol=['RR','LL'], scan=0, datacol='data', size=npix*res, res=res, sigma_image=threshold, searchtype='imageall', gainfile=gainfile, bpfile=bpfile, filtershape=filtershape, savecands=False, flagmode=flagmode, nthreads=min(14,ndm))
#d = leanpipet.pipe_thread(filename=msfile, nints=nints, nskip=0, iterint=iterint, spw=[0,1], chans=chans, dmarr=dmarr, fwhmsurvey=0.5, fwhmfield=0.5, selectpol=['RR','LL'], scan=0, datacol='corrected_data', size=npix*res, res=res, sigma_image=threshold, searchtype='imageall', telcalfile=telcalfile, filtershape=filtershape, savecands=False, flagmode=flagmode, nthreads=min(14,ndm))
