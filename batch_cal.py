#!/usr/bin/env python
#
# Script to be run within casapy for batch preparation of calibration solutions
#
# Usage:
# batch_cal.py asdmfile msfileroot
#
# asdmfile is input asdm, msroot is converted ms filename root (includes '.ms', but stuff gets added before that).
# imagetest is 0/1 to say whether imaging test should be run on cal scans afterwards
# gainscan and bpscan are optional override of scans if appropriate intents are not found. comma-delimited string list of scan numbers

import leanpipet, applycals, qrtt
import pylab as p
import numpy as n
import sys, os

# parse arguments
arg0 = sys.argv.index('batch_cal.py')
print 'Args:', sys.argv[arg0:]
asdmfile = sys.argv[arg0+1]
msfile = sys.argv[arg0+2]
imagetest = int(sys.argv[arg0+3])
plotcals = int(sys.argv[arg0+4])
if len(sys.argv[arg0+5:]):
    gainstr = sys.argv[arg0+5]
else:
    gainstr=''
if len(sys.argv[arg0+6:]):
    bpstr = sys.argv[arg0+6]
else:
    bpstr=''

# other parameters
workdir = os.getcwd()
flagmode = 'ring2medcht1.5badbp2blstd3'
threshold = 10.
refant = 'ea10'    # arbitrary for now
if asdmfile[-1] == '/':
    asdmfile = asdmfile[:-1]
telcalfile = asdmfile + '.GN'
fluxmodel = ''

# associated calibration details
chans = range(6,122)+range(134,250)       # 128 ch/spw version cutting 5% of channels from each edge plus bad chans at bottom
#spw0 = '0~1:25~35'
#spw1 = '0~1:5~60'    # avoids some RFI
spw0 = '0~1:60~70'
spw1 = '0~1:6~122'    # avoids some RFI
fps = 200  # frames per second. needed to find middle int of 1s data based on int number of high fps data.

allscans = qrtt.prep(asdmfile, workdir=workdir, intentfilter='')    # set up asdm and telcal files, get scans of interest (tuple with scannum, name, nints)

# find scans needed for gain and (optionally) bandpass
if not gainstr:
    gainints = []
    for scan in allscans:
        if 'CALIBRATE_PHASE' in scan[3] and 'BANDPASS' not in scan[3]:
            gainstr = gainstr+str(scan[0])+','
            gainints.append(scan[2])
    gainstr = gainstr[:-1]
    print 'Gain calibration with scans %s' % gainstr
else:
    gainints = []
    for scan in allscans:
        if 'CALIBRATE_PHASE' in scan[3] and 'BANDPASS' not in scan[3]:
            gainints.append(scan[2])

if not bpstr:
    for scan in allscans:
        if 'BANDPASS' in scan[3]:
            bpstr = bpstr+str(scan[0])+','
    bpstr = bpstr[:-1]
    if bpstr == '':
        bpstr = gainstr
    print 'BP calibration with scans %s' % bpstr

# assume data in place, look for cal tables 
g0name = msfile[:-3]+'.g0'
b1name = msfile[:-3]+'.b1'
g1name = msfile[:-3]+'.g1'
if fluxmodel:        # if doing flux scale, define name
    fluxfile = g0name
else:
    fluxfile = ''

if bpstr:   # doing full bp cal first
    bpms = msfile[:-3] + '_bp.ms'
    bpms2 = qrtt.asdm2ms(asdmfile, bpms, bpstr, inttime='1s')   # integrate down to 1s during split
    if not os.path.exists(b1name):

        # flag data with quack, rflag, and clip of zeros
        print 'Starting flagging of bp file...'
        flagdata(vis=bpms2, mode='unflag')
        flagdata(vis=bpms2, mode='manual', antenna='ea06')
        flagdata(vis=bpms2, mode='shadow')
        flagdata(vis=bpms2, mode='clip', clipzeros=True)
        flagdata(vis=bpms2, mode='rflag')
        flagdata(vis=bpms2, mode='extend', growaround=True, extendpols=True)
        flagdata(vis=bpms2, mode='quack', quackinterval=15)
        flags = flagdata(vis=bpms2, mode='summary')
        print 'Gain flag summary:'
        print flags

        # cal narrow for each spw, then apply for bandpass
        print 'Starting bp cal...'
        if fluxmodel:      # if a flux model is available, apply it to bp calibrator
            print 'Found flux model for BP calibrator. Applying...'
            setjy(bpms2, modimage=fluxmodel, scalebychan=True)
        gaincal(vis=bpms2, caltable=g0name, preavg=1, solint='inf', spw=spw0, refant=refant, minsnr=3, calmode='ap')
        bandpass(vis=bpms2, caltable=b1name, gaintable=g0name, solint='inf', bandtype='BPOLY', combine='scan', degamp=4, degphase=1, refant=refant, minsnr=3, maskedge=5)

if gainstr:    # follow with gaincal
    # fill gain ms
    if bpstr == gainstr:   # if bp and gain scans are the same, use the same ms file
        gainms2 = msfile[:-3] + '_gain_s' + gainstr + '.ms'
        if not os.path.exists(gainms2):
            os.symlink(bpms2, gainms2)

        gainms = msfile[:-3] + '_gain.ms'
        gainms2 = qrtt.asdm2ms(asdmfile, gainms, gainstr, inttime='1s')   # integrate down to 1s during split
    else:
        gainms = msfile[:-3] + '_gain.ms'
        gainms2 = qrtt.asdm2ms(asdmfile, gainms, gainstr, inttime='1s')   # integrate down to 1s during split

    if not os.path.exists(g1name):
        # flag data with quack, rflag, and clip of zeros
        print 'Starting flagging of gain file...'
        flagdata(vis=gainms2, mode='unflag')
        flagdata(vis=gainms2, mode='manual', antenna='ea06')
        flagdata(vis=gainms2, mode='shadow')
        flagdata(vis=gainms2, mode='manual', uvrange='<6klambda')
        flagdata(vis=gainms2, mode='clip', clipzeros=True)
        flagdata(vis=gainms2, mode='rflag')
        flagdata(vis=gainms2, mode='extend', growaround=True, extendpols=True)
        flagdata(vis=gainms2, mode='quack', quackinterval=15)
        flags = flagdata(vis=gainms2, mode='summary')
        print 'Gain flag summary:'
        print flags

        print 'Starting gain cal...'
        if os.path.exists(b1name):    # if bp cal exists, use it
            gaincal(vis=gainms2, caltable=g1name, gaintable=b1name, preavg=1, solint='inf', spw=spw1, refant=refant, minsnr=3, calmode='p')
        else:
            gaincal(vis=gainms2, caltable=g1name, preavg=1, solint='inf', spw=spw1, refant=refant, minsnr=3, calmode='p')

if plotcals:
    sols = applycals.solutions(g1name, fluxfile=fluxfile)
    sols.parsebp(b1name)
    p.figure(1)
    for spw in range(len(sols.gain[0,0])):
        for pol in range(len(sols.gain[0,0,0])):
            p.plot(n.abs(sols.gain[:,:,spw,pol]),'b.')
    p.xlabel('Solution number')
    p.ylabel('Gain amp')
    p.savefig(g1name[:-3] + '_gain.png')
    p.figure(2)
    p.plot(n.abs(sols.bandpass[:,::10,0].transpose()),'b.')
    p.plot(n.angle(sols.bandpass[:,::10,1].transpose()),'r.')
    p.xlabel('Channel')
    p.ylabel('Bandpass amp (1-mean) and angle (rad)')
    p.savefig(g1name[:-3] + '_bp.png')

if imagetest:
# then image middle integrations for both CASA gain cal and telcal solutions
    for scan in range(len(gainstr.split(','))):
        try:
            d = leanpipet.pipe_thread(filename=gainms2, nints=15, nskip=gainints[scan]/fps/2, iterint=15, spw=[0,1], chans=chans, dmarr=[0], fwhmsurvey=0.5, fwhmfield=0.5, selectpol=['RR','LL'], scan=scan, datacol='data', size=0, res=0, sigma_image=threshold, searchtype='imageall', filtershape=None, savecands=False, candsfile='', flagmode=flagmode, gainfile=g1name, bpfile=b1name, fluxfile=fluxfile)
            d = leanpipet.pipe_thread(filename=gainms2, nints=15, nskip=gainints[scan]/fps/2, iterint=15, spw=[0,1], chans=chans, dmarr=[0], fwhmsurvey=0.5, fwhmfield=0.5, selectpol=['RR','LL'], scan=scan, datacol='data', size=0, res=0, sigma_image=threshold, searchtype='imageall', filtershape=None, savecands=False, candsfile='', flagmode=flagmode, telcalfile=telcalfile)
        except IndexError:
            print 'Ran out of ints in that scan?'
            continue

    if fluxmodel:
        d = leanpipet.pipe_thread(filename=bpms2, nints=15, nskip=30, iterint=15, spw=[0,1], chans=chans, dmarr=[0], fwhmsurvey=0.5, fwhmfield=0.5, selectpol=['RR','LL'], scan=0, datacol='data', size=0, res=0, sigma_image=threshold, searchtype='imageall', filtershape=None, savecands=False, candsfile='', flagmode=flagmode, gainfile=g1name, bpfile=b1name, fluxfile=fluxfile)
        print 'Flux scale calibrator flux = %.3f' % (leanpipet.data.real.mean())
