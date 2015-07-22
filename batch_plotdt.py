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
import sys, os

arg0 = sys.argv.index('batch_plot.py')
msfile = sys.argv[arg0+1]

candspipe.mergecands_plotall('cands_' + msfile[:-3] + '*.pkl', saveroot='cands_' + msfile[:-3])
