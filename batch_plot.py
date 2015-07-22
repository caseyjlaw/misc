#!/usr/bin/env python
#
# Script to be run within casapy for batch plotting of candidates
#
# Usage:
# batch_plot.py msfileroot
#
# asdmfile is input asdm, msroot is converted ms filename root (includes '.ms', but stuff gets added before that).
# next is scan numbers (0 based, comma-delimited), MS file name will use MS scan number instead. iterates over scans.

import candspipe
import os, glob, argparse

parser = argparse.ArgumentParser()
parser.add_argument("msroot", help="root ms name for files")
parser.add_argument("--remove", help="dict defines times to remove from summary products as scan:[int0,int1]", type=eval, default={})
parser.add_argument("--threshold", help="SNR threshold to filter single cand inspection", type=float, default=0)
parser.add_argument("--candnums", help="comma-delimited list of candnums for individual inspection")
parser.add_argument("--reproducecand", help="flag to reproduce candidate", action='store_true')
parser.add_argument("--plotcand", help="flag to plot candidate", action='store_true')
args = parser.parse_args(); 
if args.candnums:
    candnums = [int(i) for i in args.candnums.split(',')]
else:
    candnums = []
threshold = args.threshold

msroot = args.msroot.split('.')[0]  # assure that suffix is excluded
mergec = 'cands_%s_mergec.pkl' % (msroot)

# merge cands (will not overwrite)
candspipe.mergecands(glob.glob('cands_%s_s*.pkl' % (msroot)))

# make new clean merge file and corresponding plots. can iteratively remove bad times. will not overwrite, unless remove specified
if ( (not os.path.exists(mergec)) or (len(args.remove.keys()) != 0) ):
    candspipe.mergeremove('cands_%s_merge.pkl' % (msroot), mergec, remove=args.remove)
    candspipe.mergeplot(mergec)
else:
    print '%s exists. Not regenerating mergc or plots.' % (mergec)

# reproduce candidates and make individual candidate plots
if threshold:
    cands = candspipe.filter_candidates(mergec, threshold)
    if len(cands):
        print '**Candidates**'
        for i in range(len(cands[0])):
            print cands[0][i], cands[1][i], cands[2][i]

# extract measurement sets for candidates. candspipe.extractms does not overwrite.
for candnum in candnums:
    outms = candspipe.extractms(mergec, threshold, candnum)

    if args.reproducecand:
        d = candspipe.reproducecand(mergec, threshold, candnum, filename=outms, searchtype='image', sigma_image=5.5)

        if args.plotcand:
            ims, snrs, ints = candspipe.reproducecanddata(d)
            if len(ints) == 1:
                print 'Yay! Got one candidate as expected (SNR=%.1f).' % snrs[0]
                im, data = candspipe.reproducecanddata(d, candint=ints[0])
                candspipe.singleplot2(mergec, threshold, candnum, im, data, 'cands_'+outms[:-3]+'.png')
            else:
                print 'Got snrs %s and ints %s. WTF?' % (snrs, ints)
