#!/usr/bin/env python

import tasklib as tl
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("msfile", help="Input ms file name")
parser.add_argument("msfile2", help="Output ms file name")
parser.add_argument("scan", help="CASA string to select scans (e.g., 1~5; inclusive)")
args = parser.parse_args()
if args.msfile:
    msfile = args.msfile.rstrip('/')
if args.msfile2:
    msfile2 = args.msfile2.rstrip('/')

cfg = tl.SplitConfig()  # configure split
cfg.vis = msfile
cfg.out = msfile2
cfg.col = 'data'
cfg.scan=args.scan  # discard autos
tl.split(cfg)  # run task
