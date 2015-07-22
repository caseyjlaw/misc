#!/usr/bin/env python
#
# claw 2oct2012
#
# Script to plot spectrograms of Sgr A* VLA monitoring

import pylab as p
import numpy as n
import os, pickle

#bbn = 'bb18'
#msname = '12A-339_sb9826039_1.56065.225594965275_' + bbn + '.ms'
#fluxname = 'flux_12A-339_sb9826039_1.56065.225594965275_' + bbn + '.txt'
#pklname = 'flux_12A-339_sb9826039_1.56065.225594965275_' + bbn + '.pkl'
msname = '12A-339_sb9826039_1.56065.225594965275.ms'
fluxname = 'flux_12A-339_sb9826039_1.56065.225594965275.txt'
pklname = 'flux_12A-339_sb9826039_1.56065.225594965275.pkl'
#fluxname = 'flux.txt'
#pklname = 'flux.pkl'

# build lists
field = []; scan = []; spw = []; peak = []; err = []; coord = []
time = []; freq = []

reference_frequencies = n.array([  8.33200000e+09,   8.46000000e+09,   3.74780000e+10,  3.76060000e+10,   3.77340000e+10,   3.78620000e+10, 3.79900000e+10,   3.81180000e+10,   3.82460000e+10,  3.83740000e+10,   2.69680000e+10,   2.70960000e+10,   2.72240000e+10,   2.73520000e+10,   2.74800000e+10,  2.76080000e+10,   2.77360000e+10,   2.78640000e+10,  4.79880000e+10,   4.81160000e+10,   4.82440000e+10,  4.83720000e+10,   4.85000000e+10,   4.86280000e+10,  4.87560000e+10,   4.88840000e+10,   3.90380000e+10, 3.91660000e+10,   3.92940000e+10,   3.94220000e+10,  3.95500000e+10,   3.96780000e+10,   3.98060000e+10,  3.99340000e+10,   2.54880000e+10,   2.56160000e+10,  2.57440000e+10,   2.58720000e+10,   2.60000000e+10,  2.61280000e+10,   2.62560000e+10,   2.63840000e+10,  1.84880000e+10,   1.86160000e+10,   1.87440000e+10,  1.88720000e+10,   1.90000000e+10,   1.91280000e+10,  1.92560000e+10,   1.93840000e+10])

if os.path.exists(pklname):
    print 'Found pkl file...'
    pklfile = open(pklname,'r')
    (field, scan, spw, peak, err, coord, time, freq) = pickle.load(pklfile)
else:
    print 'Generating new pkl file...'
    ms.open(msname)
    scinfo = ms.getscansummary()
    spwinfo = ms.getspectralwindowinfo()

    # parse flux text file
    fl = file.readlines(open(fluxname, 'r'))
    for i in range(len(fl)):
        (name, peakstr, errstr, coord1, coord2) = fl[i].split(',')
        field.append(int(name.split('_')[0][5:]))
        scan.append(int(name.split('_')[1][4:]))
        peak.append(float(peakstr[6:]))
        err.append(float(errstr[5:]))
        coord.append((int(coord1[10:]),int(coord2[:3])))
        time.append(scinfo[str(scan[i])]['0']['BeginTime'])

        spwtmp = name.split('_')[2]
        if '~' in spwtmp:
            spw.append(int(spwtmp[3:].split('~')[0]))  # skip 'spw' then take first sb
            freq.append(spwinfo[str(spw[i])]['Chan1Freq'])
        elif 'bb' in spwtmp:
            bb = int(spwtmp[2:])    # skip 'bb' then take bb
            spw.append(bb)
            if bb == 18: spwn = '45'
            if bb == 25: spwn = '34'
            if bb == 27: spwn = '10'
            if bb == 37: spwn = '2'
            if bb == 39: spwn = '28'
            if bb == 48: spwn = '18'
            freq.append(spwinfo[spwn]['Chan1Freq'])

    field = n.array(field); scan = n.array(scan); spw = n.array(spw); peak = n.array(peak); err = n.array(err); coord = n.array(coord)
    time = n.array(time); freq = n.array(freq)
    pklfile = open(pklname, 'w')
    pickle.dump((field, scan, spw, peak, err, coord, time, freq), pklfile)
    pklfile.close()
    

def pllc(fn, sn, coordtol=4):
    """ Plot lightcurve for field number and spectral window number
    sn can be spw number or freq range tuple in GHz
    coordtol defines box of +-coordtol in peak pixel to plot
    """

    print 'Plotting...'
    mjdref = 56065

    centerpix = n.median(coord[n.where(field == fn)])
    if (type(sn) == type(0)):
        wh = n.where( (field == fn) & (spw == sn) & (peak >= -9999) & (peak <= 9999) & (coord[:,0] >= centerpix - coordtol) & (coord[:,0] <= centerpix + coordtol) & (coord[:,1] >= centerpix - coordtol) & (coord[:,1] <= centerpix + coordtol) )
    elif (type(sn) == type( (0,0) )):
        wh = n.where( (field == fn) & (freq/1e9 >= sn[0]) & (freq/1e9 <= sn[1]) & (peak >= -9999) & (peak <= 9999) & (coord[:,0] >= centerpix - coordtol) & (coord[:,0] <= centerpix + coordtol) & (coord[:,1] >= centerpix - coordtol) & (coord[:,1] <= centerpix + coordtol) )
        
    if len(wh[0]) > 0:
        p.errorbar((time[wh]-mjdref)*24, peak[wh], yerr=err[wh], fmt='.', label=str(freq[wh][0]/1e9)+' GHz')
        p.legend()
        p.xlabel('Time (hrs)')
        p.ylabel('Flux (Jy)')
    else:
        print 'No measurements meet those criteria.'
    

def plsp(fn, mins, maxs, coordtol=4):
    """ Plot spectrum for field number and scan range
    coordtol defines box of +-coordtol in peak pixel to plot
    """

    print 'Plotting...'

    centerpix = n.median(coord[n.where(field == fn)])
    wh = n.where( (field == fn) & (scan >= mins) & (scan <= maxs) & (peak >= -9999) & (peak <= 9999) & (coord[:,0] >= centerpix - coordtol) & (coord[:,0] <= centerpix + coordtol) & (coord[:,1] >= centerpix - coordtol) & (coord[:,1] <= centerpix + coordtol) )

    if len(wh[0]) > 0:
        p.errorbar(freq[wh]/1e9, peak[wh], yerr=err[wh], fmt='.')
        p.xlabel('Freq (GHz)')
        p.ylabel('Flux (Jy)')
    else:
        print 'No measurements meet those criteria.'


def im(fn, vmin=0, vmax=0, coordtol=4):
    """ Plot 2d spectrogram of flux versus time and frequency. Needs to grid values.
    """

    print 'Plotting...'
    mjdref = 56065

    # trim arrays to needed field
    centerpix = n.median(coord[n.where(field==fn)])
    wh = n.where( (field == fn) & (peak >= -9999) & (peak <= 9999) & (coord[:,0] >= centerpix - coordtol) & (coord[:,0] <= centerpix + coordtol) & (coord[:,1] >= centerpix - coordtol) & (coord[:,1] <= centerpix + coordtol) )
    if len(wh[0]) == 0:
        print 'No measurements meet those criteria.'
    else:
        time2 = time[wh]
        freq2 = freq[wh]
        peak2 = peak[wh]
    
        tbins = n.unique( (time2 - mjdref)*24 )
        fbins = n.unique(freq2/1e9)
        peakarr = n.zeros( (len(fbins), len(tbins)) )

        for i in range(len(freq2)):      # iterate over all (trimmed) data
            t = (time2[i] - mjdref)*24    # time in hours
            f = freq2[i]/1e9             # freq in GHz
            peakarr[ n.where(f == fbins)[0][0], n.where(t == tbins)[0][0] ] = peak2[i]

        if (vmax==0) & (vmin==0):
            vmax = n.max(peak2)
            vmin = n.min(peak2)
            print 'vmin, vmax = ', vmin, vmax

        p.clf()
        p.imshow(peakarr, vmin=vmin, vmax=vmax, aspect='auto', origin='lower', interpolation='nearest')
        p.colorbar()
        p.xlabel('Time (hrs)')
        p.ylabel('Freq (GHz)')
        p.xticks(n.arange(len(tbins))[::10], n.round(tbins[::10],2), rotation=25)
        p.yticks(n.arange(len(fbins)), fbins)

def lccomp(sn, coordtol=4):

    pklsrc = 'flux_12A-339_sb9826039_1.56065.225594965275_f6c.pkl'
    pklref = 'flux_12A-339_sb9826039_1.56065.225594965275_f7c.pkl'
    pklname = 'flux_12A-339_sb9826039_1.56065.225594965275.pkl'
    
    if (os.path.exists(pklsrc) & os.path.exists(pklref)):
        print 'Found pkl files...'
        pkl1 = open(pklsrc,'r')
        (sfield, sscan, sspw, speak, serr, scoord, stime, sfreq) = pickle.load(pkl1)
        pkl2 = open(pklref,'r')
        (rfield, rscan, rspw, rpeak, rerr, rcoord, rtime, rfreq) = pickle.load(pkl2)
    elif os.path.exists(pklname):
        print 'Found pkl file...'
        pkl = open(pklname,'r')
        (field, scan, spw, peak, err, coord, time, freq) = pickle.load(pkl)
        f6 = n.where(field==6)
        f7 = n.where(field==7)
        sfield = field[f6]; sscan = scan[f6]; sspw = spw[f6]; speak = peak[f6]; serr = err[f6]; scoord = coord[f6]; stime = time[f6]; sfreq = freq[f6]
        rfield = field[f7]; rscan = scan[f7]; rspw = spw[f7]; rpeak = peak[f7]; rerr = err[f7]; rcoord = coord[f7]; rtime = time[f7]; rfreq = freq[f7]

    else:
        print 'Can\'t find both pkl files...'

    print 'Plotting...'
    mjdref = 56065

    centerpix = n.median(scoord)
    swh = n.where( (sspw == sn) & (speak >= -9999) & (speak <= 9999) & (scoord[:,0] >= centerpix - coordtol) & (scoord[:,0] <= centerpix + coordtol) & (scoord[:,1] >= centerpix - coordtol) & (scoord[:,1] <= centerpix + coordtol) )
    rwh = n.where( (rspw == sn) & (rpeak >= -9999) & (rpeak <= 9999) & (rcoord[:,0] >= centerpix - coordtol) & (rcoord[:,0] <= centerpix + coordtol) & (rcoord[:,1] >= centerpix - coordtol) & (rcoord[:,1] <= centerpix + coordtol) )
    scale = speak[swh].mean()/rpeak[rwh].mean()
    meandt = n.mean([stime[swh][i+1] - stime[swh][i] for i in range(len(stime[swh])-1)])

    # for corrected lc, need to select only times with both ref and src within mean delta t
    sgood = []; rgood = []
    for i in range(len(stime[swh])):
        rg = n.where(n.abs(rtime[rwh]-stime[swh][i]) <= meandt/2)
        if len(rg[0]) == 1:
            sgood.append(i)
            rgood.append(rg[0][0])

    if len(swh[0]) > 0:
        p.figure(1)
        p.subplot(211)
        p.errorbar((stime[swh][sgood]-mjdref)*24, speak[swh][sgood]/(rpeak[rwh][rgood]/rpeak[rwh][rgood].mean()), yerr=serr[swh][sgood]/(rpeak[rwh][rgood]/rpeak[rwh][rgood].mean()), fmt='.', label='Corrected Sgr A*, '  + str(sfreq[swh][0]/1e9)+' GHz')
        p.ylabel('Flux (Jy)')
        p.legend()
        p.subplot(212)
        p.errorbar((stime[swh]-mjdref)*24, speak[swh], yerr=serr[swh], fmt='.', label='Sgr A*, '  + str(sfreq[swh][0]/1e9)+' GHz')
        p.errorbar((rtime[rwh]-mjdref)*24, scale*rpeak[rwh], yerr=scale*rerr[rwh], fmt=None, capsize=0, label='Scaled Ref')
        p.legend()
        p.xlabel('Time (hrs)')
        p.ylabel('Flux (Jy)')
        p.figure(2)
        p.errorbar((stime[swh][sgood]-mjdref)*24, speak[swh][sgood]/(rpeak[rwh][rgood]/rpeak[rwh][rgood].mean()), yerr=serr[swh][sgood]/(rpeak[rwh][rgood]/rpeak[rwh][rgood].mean()), fmt='.', label='Corrected Sgr A*, '  + str(sfreq[swh][0]/1e9)+' GHz')
        p.xlabel('Time (hrs)')
        p.ylabel('Flux (Jy)')
        p.legend()
    else:
        print 'No measurements meet those criteria.'
