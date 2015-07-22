#!/usr/bin/env python
#
# Script to be run within casapy for batch data conversion
#
# Usage:
# batch_filler.py asdmfile msfileroot scanlist
#
# asdmfile is input asdm, msroot is converted ms filename root (includes '.ms', but stuff gets added before that).
# next is scan numbers (0 based, comma-delimited), MS file name will use MS scan number instead. iterates over scans.
# similarly, the dmbins are given as comma-delimited list. a single run will use pairs of these numbers. "0,30,60" will run with 0-30, then 30-60.
# iteration over dms happens within iteration over scans.
# next two arguments define range of dmbins to search (0-32 here).
# working directory should have asdm file and telcal file named asdmfil+".GN"

import leanpipedt, parseasdm
import sys, os, argparse
import numpy as n

parser = argparse.ArgumentParser()
parser.add_argument("asdmfile", help="input asdm file name")
parser.add_argument("msfile", help="root of output ms file name")
parser.add_argument("scan", help="scan to select from asdmfile")
args = parser.parse_args(); asdmfile = args.asdmfile.rstrip('/'); msfile = args.msfile; scan = args.scan

# run prep and search
goodscans = parseasdm.getscans(asdmfile, namefilter='J1911')
print goodscans
#msfile2 = parseasdm.asdm2ms(asdmfile, msfile, str(goodscans[int(scan)][0]))   # if scans indexed zero-based
msfile2 = parseasdm.asdm2ms(asdmfile, msfile, scan)    # if scans indexed by scan number
d = leanpipedt.pipe_thread(filename=msfile2, nints=100, nskip=0, iterint=100, spw=[0,1], chans=range(64), dmarr=[0], fwhmsurvey=0.5, fwhmfield=0.5, selectpol=['RR','LL'], scan=0, datacol='data', size=25600, res=50, sigma_image=10, searchtype='', filtershape=None, savecands=False, candsfile='', flagmode='')    # make pkl file only

# tell node manager that we're done here...
finishedfile = 'tracking_dir/' + os.uname()[1] + '.ready_' + str(n.random.randint(100))
open(finishedfile, 'a').close()
