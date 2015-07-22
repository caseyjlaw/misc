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

import leanpipedt, applycals2, parseasdm
import tasklib as tl
import pylab as p
import numpy as n
import sys, os, argparse

parser = argparse.ArgumentParser()
parser.add_argument("asdmfile", help="Input asdm file name")
parser.add_argument("msfile", help="Root ms name for output and cal file name")
parser.add_argument("--gainstr", help="Gain scans (comma-delimited)")
parser.add_argument("--bpstr", help="BP scans (comma-delimited)")
parser.add_argument("--image", help="Image representative data as test", default=1)
parser.add_argument("--plot", help="Plot gain/bandpass solutions", default=1)
args = parser.parse_args(); asdmfile = args.asdmfile; msfile = args.msfile
if args.gainstr:
    gainstr = args.gainstr
else:
    gainstr = ''
if args.gainstr:
    bpstr = args.bpstr
else:
    bpstr = ''

# calibration parameters
refant = ['ea17', 'ea25']
antsel = '' # '!ea08;!ea28;!ea19;!ea26;!ea04;!ea22'
uvrange = '' #'<5klambda'
delaycal = False     # if doing delay cal => can only be applied in casa
fluxmodel = ''
#fluxname = ''    # can use scan instead?
spw0 = '0~1:60~70'  # initial gain solution
spw1 = '0~1:6~122'  # final gain/bp solution
flags1 = ["mode='unflag'", "mode='shadow'", "mode='clip' clipzeros=True", "mode='manual' antenna='ea11,ea19'", "mode='rflag'", "mode='extend' growaround=True extendpols=True", "mode='quack' quackinterval=20", "mode='summary'"]

# file names
antposname = msfile[:-3]+'.antpos'   # antpos
delayname = msfile[:-3]+'.delay'   # delay cal
g0name = msfile[:-3]+'.g0'   # initial gain correction before bp
b1name = msfile[:-3]+'.b1'   # bandpass file
g1name = msfile[:-3]+'.g1'   # gain cal per scan
g2name = msfile[:-3]+'.g2'   # flux scale applied
telcalfile = asdmfile + '.GN'
gainms = msfile[:-3] + '_gain.ms'
bpms = msfile[:-3] + '_bp.ms'
concatms = msfile[:-3] + '_concat.ms'

# imaging parameters
if args.image:
    import qimg_cython
    flagmode = 'ring2medcht1.5badbp2blstd3'
    threshold = 10.
    chans = range(6,122)+range(134,250)       # 128 ch/spw version cutting 5% of channels from each edge plus bad chans at bottom
    fps = 200  # frames per second. needed to find middle int of 1s data based on int number of high fps data.
    npix = 1024
    res = 58

if asdmfile[-1] == '/':
    asdmfile = asdmfile[:-1]
workdir = os.getcwd()
allscans = parseasdm.getscans(asdmfile)

# If gain and bp not provided, find them
if not gainstr:
    gainints = []
    for scan in allscans:
        if 'CALIBRATE_PHASE' in scan[3] and 'BANDPASS' not in scan[3]:
            gainstr = gainstr+str(scan[0])+','
            gainints.append(scan[2])
    gainstr = gainstr[:-1]
    print 'Found gain calibrators in scans %s' % gainstr
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
    print 'Found bandpass calibrators in scans %s' % bpstr

# prepare ms files
bpms2 = parseasdm.asdm2ms(asdmfile, bpms, bpstr, inttime='1')   # integrate down to 1s during split
if bpstr == gainstr:   # if bp and gain scans are the same, use the same ms file
    gainms2 = msfile[:-3] + '_gain_s' + gainstr + '.ms'
    if not os.path.exists(gainms2):
        os.symlink(bpms2, gainms2)
else:
    gainms2 = parseasdm.asdm2ms(asdmfile, gainms, gainstr, inttime='1')   # integrate down to 1s during split

# make concat file
if not os.path.exists(concatms):
    print 'Concatenating bp and gain scans...'
    tl.concat([gainms2, bpms2], concatms)

    # flag data via text file
    flfile = open('flags.txt','w')
    for flag in flags1:
        flfile.write(flag + '\n')
    flfile.close()

    print 'Flagging with these commands:'
    for ff in enumerate(open('flags.txt')): print ff[1].rstrip()

    print 'Starting flagging of bp file...'
    cfg = tl.FlaglistConfig()  # configure split
    cfg.vis = concatms
    cfg.inpfile = "flags.txt"
    tl.flaglist(cfg)  # run task

    # clean up
    os.remove('flags.txt')

# Calibrate!
if fluxmodel:
    if not os.path.exists(g0name):
        print 'Applying flux model for BP calibrator...'
        cfg = tl.SetjyConfig()
        cfg.vis = concatms
        cfg.scan = bpstr
        cfg.modimage = fluxmodel
        tl.setjy(cfg)

        print 'Starting initial gain cal...'
        cfg = tl.GaincalConfig()
        cfg.vis = concatms
        cfg.caltable = g0name
        if delaycal:
            cfg.gaintable = [antposname]   # if doing delay cal
        else:
            cfg.gaintable = []
        cfg.scan = bpstr
        cfg.gaintype = 'G'
        cfg.solint = 'inf'
        cfg.spw = spw0
        cfg.refant = refant
        cfg.minsnr = 5.
        cfg.calmode = 'p'
        cfg.antenna = antsel
        cfg.uvrange = uvrange
        tl.gaincal(cfg)

    if delaycal and not os.path.exists(delayname):   # no delay cal for now
        print 'Starting delay calibration...'
        cfg = tl.GaincalConfig()
        cfg.vis = concatms
        cfg.caltable = delayname
        cfg.gaintable = [antposname, g0name]
        cfg.scan = bpstr
        cfg.gaintype = 'K'
        cfg.solint = 'inf'
        cfg.combine = ['scan']
        cfg.spw = spw1
        cfg.refant = refant
        cfg.minsnr = 5.
        cfg.antenna = antsel
        cfg.uvrange = uvrange
        tl.gaincal(cfg)

    if not os.path.exists(b1name):
        print 'Starting bp cal...'
        cfg = tl.GaincalConfig()
        cfg.vis = concatms
        cfg.caltable = b1name
        if delaycal:
            cfg.gaintable = [antposname, g0name, delayname]    # if doing delay cal
        else:
            cfg.gaintable = [g0name]
        cfg.scan = bpstr
        cfg.spw = spw1
        cfg.gaintype = 'BPOLY'
        cfg.degamp = 5
        cfg.degphase = 2
        cfg.maskedge = 6
        cfg.solint = 'inf'
        cfg.combine = ['scan']
        cfg.solnorm = True
        cfg.refant = refant
        cfg.antenna = antsel
        cfg.uvrange = uvrange
        tl.gaincal(cfg)

    if not os.path.exists(g1name) or not os.path.exists(g2name):
        print 'Starting gain cal...'
        cfg = tl.GaincalConfig()
        cfg.vis = concatms
        cfg.caltable = g1name
        if delaycal:
            cfg.gaintable = [antposname, delayname, b1name]     # if doing delay cal
        else:
            cfg.gaintable = [b1name]
        cfg.scan = gainstr
        cfg.gaintype = 'G'
        cfg.solint = 'inf'
        cfg.spw = spw1
        cfg.refant = refant
        cfg.minsnr = 5.
        cfg.calmode='ap'
        cfg.antenna = antsel
        cfg.uvrange = uvrange
        tl.gaincal(cfg)

        print 'Transferring flux scale...'
        cfg = tl.FluxscaleConfig()
        cfg.vis = concatms
        cfg.caltable = g1name
        cfg.fluxtable = g2name
        cfg.reference = fluxname
        tl.fluxscale(cfg)

else:    # without fluxscale

# gencal not yet supported in casapy-free mode
#       if delaycal:
#        print 'Starting antpos calibration...'
#        cfg = tl.GencalConfig()
#        cfg.vis = concatms
#        cfg.caltable = antposname
#        cfg.scan = bpstr
#        cfg.caltype = 'antpos'
#        tl.gencal(cfg)

    if not os.path.exists(g0name):
        print 'Starting initial gain cal...'
        cfg = tl.GaincalConfig()
        cfg.vis = concatms
        cfg.caltable = g0name
        if delaycal:
            cfg.gaintable = [antposname]   # if doing delay cal
        else:
            cfg.gaintable = []
        cfg.scan = bpstr
        cfg.gaintype = 'G'
        cfg.solint = 'inf'
        cfg.spw = spw0
        cfg.refant = refant
        cfg.minsnr = 5.
        cfg.calmode = 'p'
        cfg.antenna = antsel
        cfg.uvrange = uvrange
        tl.gaincal(cfg)

    if delaycal and not os.path.exists(delayname):   # no delay cal for now
        print 'Starting delay calibration...'
        cfg = tl.GaincalConfig()
        cfg.vis = concatms
        cfg.caltable = delayname
        cfg.gaintable = [antposname, g0name]
        cfg.scan = bpstr
        cfg.gaintype = 'K'
        cfg.solint = 'inf'
        cfg.combine = ['scan']
        cfg.spw = spw1
        cfg.refant = refant
        cfg.minsnr = 5.
        cfg.antenna = antsel
        cfg.uvrange = uvrange
        tl.gaincal(cfg)

    if not os.path.exists(b1name):
        print 'Starting bp cal...'
        cfg = tl.GaincalConfig()
        cfg.vis = concatms
        cfg.caltable = b1name
        if delaycal:
            cfg.gaintable = [antposname, g0name, delayname]    # if doing delay cal
        else:
            cfg.gaintable = [g0name]
        cfg.scan = bpstr
        cfg.spw = spw1
        cfg.gaintype = 'BPOLY'
        cfg.degamp = 6
        cfg.degphase = 3
        cfg.maskedge = 5
        cfg.solint = 'inf'
        cfg.combine = ['scan']
        cfg.solnorm = True
        cfg.refant = refant
        cfg.antenna = antsel
        cfg.uvrange = uvrange
        tl.gaincal(cfg)

    if not os.path.exists(g1name):
        print 'Starting gain cal...'
        cfg = tl.GaincalConfig()
        cfg.vis = concatms
        cfg.caltable = g1name
        if delaycal:
            cfg.gaintable = [antposname, delayname, b1name]     # if doing delay cal
        else:
            cfg.gaintable = [b1name]
        cfg.scan = gainstr
        cfg.gaintype = 'G'
        cfg.solint = 'inf'
        cfg.spw = spw1
        cfg.refant = refant
        cfg.minsnr = 5.
        cfg.calmode='ap'
        cfg.antenna = antsel
        cfg.uvrange = uvrange
        tl.gaincal(cfg)

    if delaycal:        
        print 'Applying calibration...'
        cfg = tl.ApplycalConfig()
        cfg.vis = concatms
        cfg.gaintable = [antposname, delayname, b1name, g1name]     # if doing delay cal
        tl.applycal(cfg)

# plot
if args.plot:
    if fluxmodel:
        finalgain = g2name
    else:
        finalgain = g1name
    sols = applycals2.solutions(finalgain)
    sols.parsebp(b1name)
    p.figure(1)
    for spw in range(len(sols.gain[0,0])):
        for pol in range(len(sols.gain[0,0,0])):
            p.subplot(211)
            p.plot(n.abs(sols.gain[:,:,spw,pol]),'b-')
            p.plot(n.abs(sols.gain[:,:,spw,pol]),'b*')
            p.subplot(212)
            p.plot(n.angle(sols.gain[:,:,spw,pol]),'r-')
            p.plot(n.angle(sols.gain[:,:,spw,pol]),'r*')
    p.subplot(211)
    p.ylabel('Gain amp')
    p.subplot(212)
    p.xlabel('Solution number')
    p.ylabel('Gain phase')
    p.savefig(finalgain[:-3] + '_gain.png')
    p.figure(2)
    p.plot(n.abs(sols.bandpass[:,::10,0].transpose()),'b.')
    p.plot(n.angle(sols.bandpass[:,::10,1].transpose()),'r.')
    p.xlabel('Channel')
    p.ylabel('Bandpass amp (1-mean) and angle (rad)')
    p.savefig(finalgain[:-3] + '_bp.png')

# optionally image
if args.image:
    if delaycal:   # if doing delay cal, good data in corrected_data column
        datacol = 'corrected_data'
        gainfile = ''
        bpfile = ''
    else:         # else, applying on the fly
        datacol = 'data'
        if fluxmodel:
            gainfile = g2name
        else:
            gainfile = g1name
        bpfile = b1name

    for scan in range(len(gainstr.split(','))):
        try:
            d = leanpipedt.pipe_thread(filename=concatms, nints=15, nskip=gainints[scan]/fps/2, iterint=15, spw=[0,1], chans=chans, dmarr=[0], fwhmsurvey=0.5, fwhmfield=0.5, selectpol=['RR','LL'], scan=scan, datacol=datacol, size=npix*res, res=res, sigma_image=threshold, searchtype='imageall', filtershape=None, savecands=False, candsfile='', flagmode=flagmode, gainfile=gainfile, bpfile=bpfile)
            if os.path.exists(telcalfile):
                d = leanpipedt.pipe_thread(filename=concatms, nints=15, nskip=gainints[scan]/fps/2, iterint=15, spw=[0,1], chans=chans, dmarr=[0], fwhmsurvey=0.5, fwhmfield=0.5, selectpol=['RR','LL'], scan=scan, datacol=datacol, size=npix*res, res=res, sigma_image=threshold, searchtype='imageall', filtershape=None, savecands=False, candsfile='', flagmode=flagmode, telcalfile=telcalfile)
        except IndexError:
            print 'Ran out of ints in that scan?'
            continue

    if fluxmodel:
        print
        print 'Imaging flux calibrator...'
        d = leanpipedt.pipe_thread(filename=bpms2, nints=15, nskip=30, iterint=15, spw=[0,1], chans=chans, dmarr=[0], fwhmsurvey=0.5, fwhmfield=0.5, selectpol=['RR','LL'], scan=0, datacol=datacol, size=npix*res, res=res, sigma_image=threshold, searchtype='imageall', filtershape=None, savecands=False, candsfile='', flagmode=flagmode, gainfile=gainfile, bpfile=bpfile)
        print 'Flux scale calibrator flux = %.3f' % (leanpipedt.data.real.mean())

    print 'Comparing beam to single image...'
    im = qimg_cython.imgonefullxy(n.outer(leanpipedt.u[0], d['freq']/d['freq_orig'][0]), n.outer(leanpipedt.v[0], d['freq']/d['freq_orig'][0]), leanpipedt.data[0], npix*res, npix*res, res)
    beam = qimg_cython.beamonefullxy(n.outer(leanpipedt.u[0], d['freq']/d['freq_orig'][0]), n.outer(leanpipedt.v[0], d['freq']/d['freq_orig'][0]), leanpipedt.data[0], npix*res, npix*res, res)
    p.figure(3, figsize=(12,6))
    p.subplot(121)
    p.imshow(im, interpolation='nearest')
    p.colorbar()
    p.title('Image (1 int)')
    p.subplot(122)
    p.imshow(beam, interpolation='nearest')
    p.colorbar()
    p.title('Beam (1 int)')
    p.savefig(msfile[:-3] + '_beam.png')
