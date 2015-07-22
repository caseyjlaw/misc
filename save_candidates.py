#! /usr/bin/env python

"""
save_candidates.py --- example pipeline script using tpipe.py
Searches visibility data for transient using bispectrum method, then image to localize pulse and beamform to generate spectrogram.
"""

import tpipe
import numpy as n

# params
bgwindow = 4
sigma = 5
size = 80000
res = 400

nskip=15000

obs = tpipe.pipe_msdisp('12A-336_1b_j0628.ms',nints=300,nskip=nskip,selectpol=['RR','LL'],spw=[0],datacol='data', scan=2)
obs.make_bispectra(bgwindow=bgwindow)
bispname = 'bisp_j0628_sk'+str(nskip)+'.png'
(basnr, bastd, cands) = obs.detect_bispectra(sigma=sigma, save=bispname)
cands = n.array(cands)
for (int,dmind) in cands.transpose():
    print 'Generating plots for candidate with (int, dmind) = (' + str(int) + ',' + str(dmind) + ').'
    imname = 'im_j0628_sk'+str(nskip)+'_d'+str(dmind)+'i'+str(int)+'.png'
    im=obs.imagetrack(obs.tracksub(dmind, int, bgwindow=bgwindow), int=int, pol='i', size=size, res=res, clean=1, save=imname)
    imsize = size/res
    pixelscale = 1./size
    peakm, peakl = n.where(im == im.max())
    shiftm = imsize/2 - peakm   # indexed starting at top...
    shiftl = peakl - imsize/2   # indexed starting at left.
    obs.phaseshift(-shiftl*pixelscale, -shiftm*pixelscale)
    specname = 'spec_j0628_sk'+str(nskip)+'_d'+str(dmind)+'i'+str(int)+'.png'
    obs.spec(ind=[int], save=specname)
