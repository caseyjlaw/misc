#!/usr/bin/env python
#
# Script to be run within casapy for refined candidate processing and searching
#
# Usage:
# batch_leanpipet.py asdmfile msfileroot scanlist dmlist
#
# asdmfile is input asdm, msroot is converted ms filename root (includes '.ms', but stuff gets added before that).
# next is scan numbers (0 based, comma-delimited), MS file name will use MS scan number instead. iterates over scans.
# similarly, the dmbins are given as comma-delimited list. a single run will use pairs of these numbers. "0,30,60" will run with 0-30, then 30-60.
# iteration over dms happens within iteration over scans.
# next two arguments define range of dmbins to search (0-32 here).
# working directory should have asdm file and telcal file named asdmfil+".GN"

import leanpipet_dev
import sys, os, qrtt
import cPickle as pickle

#!@# arg0 = sys.argv.index('batch_refine.py')
#!@# asdmfile = sys.argv[arg0+1]
inpicklefile='cands_13B-409_13sep19v1_s35_merge.pkl' #!@#
#!@# inpicklefile = sys.argv[arg0+1]
#!@# operation_mode = sys.argv[arg0+2]
#scans = sys.argv[arg0+3]
#scans = [int(i) for i in scans.split(',')]
#dms = sys.argv[arg0+4]
#dms = [int(i) for i in dms.split(',')]
nthreads = 14
threshold = 6.5         # 6.5sigma should produce 7.5 false+ for 512x512 imgs, 30 dms, and 23900 ints (i.e., a typical run on a node)
workdir = os.getcwd()
nskip = 300

# Read header and candidate list from the pickle file
inpickle = open(inpicklefile,'rb')
head = pickle.load(inpickle)
(loc,prop) = pickle.load(inpickle)
# cand[0]: location (unused/boxcar width,unused/phasing reference value,integration num, dm bin)
# cand[1]: "properties"


# Extract required header information
old_dmarr = head['dmarr']
t_samp = head['inttime']
freqarr = head['freq']
b_chan = freqarr[1]-freqarr[0]
lo_freq = freqarr[0]
hi_freq = freqarr[-1]
bw = hi_freq - lo_freq
ctr_freq = (hi_freq + lo_freq)/2

print 'Found these parameters:\ntsamp %f\nchansize %f\nfreqrange %f - %f\nCtr freq %f\nTotal BW %f' % (t_samp,b_chan,lo_freq,hi_freq,ctr_freq,bw)


#!H NOW HERE, NEED TO:
# X Read the pickle
# X Parse the pickle; read dmarr, header freqs/tsamps, candidate SNRs and times. Also can read selected chans etc.
# X Calculate refined DM trial steps between iDM-4 and iDM+4
# - Set up to read the ASDM file...
# - Run leanpipet as usual with refined DM trials, and only for a fragment of the data.

# Read candidate DM and iteration data 
usecand = loc[0]
cand_iter = usecand[2] # Iteration number
cand_i_dm = usecand[3] # DM index
cand_dm = old_dmarr[cand_i_dm]


# Determine new DM search range for candidate
# Here we've semi-arbitrarily set it to search +/- 4 steps from the detected DM
i_dm_range = 4 
if cand_i_dm < i_dm_range:
    dm_lo = old_dmarr[0]
else:
    dm_lo = old_dmarr[cand_i_dm - i_dm_range]

if cand_i_dm > len(old_dmarr):
    dm_hi = old_dmarr[-1]
else:
    dm_hi = old_dmarr[cand_i_dm + i_dm_range]


print 'Building DM array between DM %f -> %f' % (dm_lo,dm_hi)

dmarr = calc_dmlist(dm_lo,dm_hi,t_samp,0.0,b_chan, ctr_freq,(hi_freq-lo_freq)/b_chan,1.1)


print 'New DM array is:'
print dmarr

exit

# !!! Need to calculate candidate delay here to pick range of time to use.




telcalfile=asdmfile + '.GN'
gainfile = msfile[:-3] + '.g1'
bpfile = msfile[:-3] + '.b1'
flagmode = 'medcht1.5badbp2'
if asdmfile[-1] == '/':
    asdmfile = asdmfile[:-1]
nskip = 300
iterint = 100

# run prep and search
goodscans = qrtt.prep(asdmfile, workdir=workdir, intentfilter=intentfilter)    # set up asdm and telcal files, get scans of interest (tuple with scannum, name, nints)

if os.path.exists(telcalfile) or os.path.exists(gainfile):    # need to have telcal ready to proceed
    for scan in scans:
        for i in range(len(dms)-1):
            dmbin0 = dms[i]
            dmbin1 = dms[i+1]
            dmarr = dmarrall[dmbin0:dmbin1]
            msfile2 = qrtt.asdm2ms(asdmfile, msfile, str(goodscans[scan][0]))            # convert asdm to ms for a scan. returns new ms file name
            candsfile = 'cands_' + msfile2[:-3] + '_dm' + str(dmbin0) + '-' + str(dmbin1) + '.pkl'     # cands filename
            print 'Searching %s for dmbins from %d to %d' % (msfile2, dmbin0, dmbin1)
            if os.path.exists(candsfile):
                print '%s candidate file already exists. Stopping.' % candsfile
            else:
                d = leanpipet.pipe_thread(filename=msfile2, nints=goodscans[scan][2]-(nskip+iterint), nskip=nskip, iterint=iterint, spw=[0,1], chans=chans, dmarr=dmarr, fwhmsurvey=0.5, fwhmfield=0.5, selectpol=['RR','LL'], scan=0, datacol='data', size=0, res=0, sigma_image=threshold, searchtype='imageall', gainfile=gainfile, bpfile=bpfile, filtershape=None, savecands=True, candsfile=candsfile, flagmode=flagmode, nthreads=nthreads)    # skip 200 ints to avoid mystery data at start

# tell node manager that we're done here...
if os.uname()[1] == 'jvla.local':
    pass
else:
    finishedfile = 'tracking_dir/' + os.uname()[1] + '.finished'
    open(finishedfile, 'a').close()
