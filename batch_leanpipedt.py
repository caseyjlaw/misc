#!/usr/bin/env python
#
# Python script for batch conversion and transient searching

import leanpipedt, parseasdm
import sys, os, argparse
import numpy as n

parser = argparse.ArgumentParser()
parser.add_argument("asdmfile", help="input asdm file name")
parser.add_argument("msfile", help="root ms name for output and cal file name")
parser.add_argument("scan", help="scan to select from asdmfile")
parser.add_argument("--dms", help="dmranges to search over (comma-delimited; 0,5,10 => [0,5], [5,10])")
parser.add_argument("--nthreads", help="number of threads to use during search")
parser.add_argument("--nints", help="number of integrations to search (good for interactive testing)")
args = parser.parse_args(); asdmfile = args.asdmfile; msfile = args.msfile; scan = args.scan
if args.dms:
    dms = [int(i) for i in args.dms.split(',')]
else:
    dms = [-1]
if args.nthreads:
    nthreads = int(args.nthreads)
else:
    nthreads = 14

# associated processing details
threshold = 6.0    # 6.5sigma should produce 7.5 false+ for 512x512 imgs, 30 dms, and 23900 ints
nskip = 300      # data reading parameters
iterint = 200    # helps if multiple of 2,4,8. must not allow much fringe rotation.
spw = [0,1]
chans = range(6,122)+range(134,250)
filtershape = 'z'      # time filter
res = 58   # imaging parameters. set res=size=0 to define from uv coords
npix = 512
size = npix*res
searchtype = 'lowmem'
secondaryfilter = 'fullim'
dtarr = [1,2,4,8]   # integer to integrate in time for independent searches
dmarrall = [0,19.2033,38.4033,57.6025,76.8036,96.0093,115.222,134.445,153.68,172.93,192.198,211.486,230.797,250.133,269.498,288.894,308.323,327.788,347.292,366.837,386.426,406.062,425.747,445.484,465.276,485.125,505.033,525.005,545.042,565.147,585.322,605.571,625.896,646.3,666.786,687.355,708.012,728.759,749.598,770.532,791.565,812.699,833.936,855.28,876.733,898.299,919.979,941.778,963.697,985.741,1007.91,1030.21,1052.64,1075.21,1097.92,1120.76,1143.76,1166.9,1190.19,1213.63,1237.23,1260.99,1284.92,1309.01,1333.27,1357.7,1382.31,1407.09,1432.06,1457.22,1482.56,1508.1,1533.83,1559.76,1585.89,1612.23,1638.77,1665.53,1692.51,1719.7,1747.12,1774.77,1802.64,1830.75,1859.1,1887.69,1916.53,1945.61,1974.95,2004.54,2034.39,2064.51,2094.9,2125.56,2156.49,2187.71,2219.21,2250.99,2283.07,2315.45,2348.13,2381.12,2414.41,2448.02,2481.94,2516.19,2550.77,2585.68,2620.92,2656.51,2692.44,2728.72,2765.36,2802.36,2839.72,2877.45,2915.55,2954.04,2992.91]
flagmode = 'standard'  # flagging
# set up cal files
gainfile = msfile[:-3] + '.g2'
bpfile = msfile[:-3] + '.b1'
telcalfile = ''

# set up ms and define nints
msfile2 = parseasdm.asdm2ms(asdmfile, msfile, scan)
goodscans = parseasdm.getscans(asdmfile)    # only consider scans that have bdfs
intlist = [sc[2] for sc in goodscans if sc[0] == int(scan)]
if len(intlist) > 0:
    if args.nints:
        nints = int(args.nints)
    else:
        nints = intlist[0]
else:
    print 'That scan is not in filtered scan list. Not a target field?'
    exit()

for i in range(len(dms)):
    dmbin0 = dms[i]
    try:
        if len(dms) == 1:
            dmbin1 = dms[i]+1
        else:
            dmbin1 = dms[i+1]
    except IndexError:
        break

    if dms[i] == -1:      # option to search all dm defined in script
        dmbin0 = 0
        dmbin1 = len(dmarrall)

    dmarr = dmarrall[dmbin0:dmbin1]

    candsfile = 'cands_' + msfile2[:-3] + '_dm' + str(dmbin0) + '-' + str(dmbin1) + '.pkl'     # cands filename
    print 'Searching %s for dmbin from %d to %d' % (msfile2, dmbin0, dmbin1)
    if os.path.exists(candsfile):
        print '%s candidate file already exists. Stopping.' % candsfile
    else:
        try:
            d = leanpipedt.pipe_thread(filename=msfile2, nints=nints, nskip=nskip, iterint=iterint, spw=spw, chans=chans, dmarr=dmarr, dtarr=dtarr, fwhmsurvey=0.5, fwhmfield=0.5, selectpol=['RR','LL'], scan=0, datacol='data', size=size, res=res, sigma_image=threshold, searchtype=searchtype, secondaryfilter=secondaryfilter, gainfile=gainfile, bpfile=bpfile, filtershape=filtershape, savecands=True, candsfile=candsfile, flagmode=flagmode, nthreads=nthreads)
        except:
            print 'Processing of %s failed with %s exception.' % (msfile2, sys.exc_info()[0])

# tell node manager that we're done here...
finishedfile = 'tracking_dir/' + os.uname()[1] + '.finished'
open(finishedfile, 'a').close()
