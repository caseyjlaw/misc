import candspipe, sys
import pylab as p
import numpy as n

# parse arguments
arg0 = sys.argv.index('compresscands.py')
print 'Args:', sys.argv[arg0:]
pkllist = sys.argv[arg0+1:-1]
outname = sys.argv[-1]

candspipe.compresscands(pkllist, outname)
