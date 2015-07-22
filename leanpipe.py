##########################################
# functional style, uses multiprocessing #
##########################################

import numpy as n
import pylab as p
import applytelcal, types
from scipy import signal
from math import ceil
import multiprocessing as mp
import string, os, pickle, ctypes
import time as timestamp
import leanpipe_cython as lib
#import leanpipe_external as lib
import qimg_cython as qimg
#import qimg
import cProfile

# set up libraries for reading and imaging visibility data
try:
    # as a backup, casa stuff can be imported if running casapy
    #        from casac import casac; # Recommended by Sanjay B. (Don't ask me why this is so! :)
    #        ms = casac.ms();
    from casa import ms
    from casa import quanta as qa
    print 'Imported CASA'
except ImportError:
    print 'No CASA or pyrap available. Can\'t read MS data.'

def calc_hexcenters(fwhmsurvey, fwhmfield, show=0):
    """ Tile a large circular area with a small circular areas. sizes are assumed to be fwhm. assumes flat sky.
    """

    large = fwhmsurvey
    small = fwhmfield
    centers = []
    (l0,m0) = (0.,0.)

    centers.append((l0,m0))
    l1 = l0-(small/2.)*n.cos(n.radians(60))
    m1 = m0-(small/2.)*n.sin(n.radians(60))
    ii = 0
    while ( n.sqrt((l1-l0)**2+(m1-m0)**2) < large/2.): 
        l1 = l1+((-1)**ii)*(small/2.)*n.cos(n.radians(60))
        m1 = m1+(small/2.)*n.sin(n.radians(60))
        l2 = l1+small/2
        m2 = m1
        while ( n.sqrt((l2-l0)**2+(m2-m0)**2) < large/2.): 
            centers.append((l2,m2))
            l2 = l2+small/2
        l2 = l1-small/2
        m2 = m1
        while ( n.sqrt((l2-l0)**2+(m2-m0)**2) < large/2.): 
            centers.append((l2,m2))
            l2 = l2-small/2
        ii = ii+1
    l1 = l0
    m1 = m0
    ii = 0
    while ( n.sqrt((l1-l0)**2+(m1-m0)**2) < large/2.): 
        l1 = l1-((-1)**ii)*(small/2.)*n.cos(n.radians(60))
        m1 = m1-(small/2.)*n.sin(n.radians(60))
        l2 = l1
        m2 = m1
        while ( n.sqrt((l2-l0)**2+(m2-m0)**2) < large/2.): 
            centers.append((l2,m2))
            l2 = l2+small/2
        l2 = l1-small/2
        m2 = m1
        while ( n.sqrt((l2-l0)**2+(m2-m0)**2) < large/2.): 
            centers.append((l2,m2))
            l2 = l2-small/2
        ii = ii+1

    delaycenters = n.array(centers)

    if len(delaycenters) == 1:
        plural = ''
    else:
        plural = 's'
    print 'For a search area of %.3f and delay beam of %.3f, we will use %d delay beam%s' % (fwhmsurvey, fwhmfield, len(delaycenters), plural)
    return delaycenters

def numpyview(arr, datatype, shape):
    """ Takes mp.Array and returns numpy array with shape of data in MS.
    """

#    return n.frombuffer(arr.get_obj()).view(n.dtype(datatype)).reshape((iterint, nbl, nchan, npol))
    return n.frombuffer(arr.get_obj(), dtype=n.dtype(datatype)).view(n.dtype(datatype)).reshape(shape)

def filtercands(cands, selectbeam=-1, selectcand=-1):
    """ Function to make dictionary that only contains beams with candidates.
    select is a list of indices that specifies a further subset of candidates.
    """

    cands2 = cands.copy()
    keys = cands2.keys()
    for i in xrange(len(keys)):
        beam = keys[i]
        if len(cands2[beam]) == 0:
            junk = cands2.pop(beam)

    if (selectbeam==-1 and selectcand==-1):
        return cands2
    else:
        keys = cands2.keys()
        for i in xrange(len(keys)):
            beam = keys[i]
            if i != selectbeam:
                junk = cands2.pop(beam)

        return {cands2.keys()[0]: [cands2[cands2.keys()[0]][selectcand]]}
    
def count(cands):
    candcount = 0
    for beam in cands.keys():
        cc = len(cands[beam]) 
        if cc:
            print '\nBeam: ', beam, '. Count:', cc
            print cands[beam]
            candcount += cc
    print '\nCand count: ', candcount
    return candcount

def summarize_cands(pklname):
    """ Reads pickle file and summarizes where candidates are in data space.
    """

    with open(pklname, 'r') as pkl:
        shift = 0
        beams = []
        d = pickle.load(pkl)
        (masterloc, masterprop, masterdata, masternoise) = pickle.load(pkl)

        goodloc = filtercands(masterloc)
        goodprop = filtercands(masterprop)
        print 'Candidate beam, location, properties:'
        for beam in goodloc.keys():
            print beam, goodloc[beam], goodprop[beam]
            beams.append(beam)
            for cand in goodloc[beam]:
                p.text(beam[0]+shift, beam[1]+shift, str(cand[0]) + '-' + str(cand[1]) + '-' + str(cand[2]), horizontalalignment='center', verticalalignment='center', fontsize=9)
                shift += 0.00003

        for beam in beams:
            p.plot([beam[0], beam[0]+shift], [beam[1], beam[1]+shift], 'k-')

    p.show()

def spec(data0, startint, endint):
    """ Function generates spectrogram data from startint to endint
    """

    return data0[startint:endint].mean(axis=1).mean(axis=2).real

def detect_bispectra(ba, d, sigma=5., tol=1.3, Q=0, show=0, save=0, verbose=0):
    """ Function to detect transient in bispectra
    sigma gives the threshold for SNR_bisp (apparent). 
    tol gives the amount of tolerance in the sigma_b cut for point-like sources (rfi filter).
    Q is noise per baseline and can be input. Otherwise estimated from data.
    save=0 is no saving, save=1 is save with default name, save=<string>.png uses custom name (must include .png). 
    """

    # using s=S/Q
    mu = lambda s: 1.  # for bispectra formed from visibilities
    sigbQ3 = lambda s: n.sqrt((1 + 3*mu(s)**2) + 3*(1 + mu(s)**2)*s**2 + 3*s**4)  # from kulkarni 1989, normalized by Q**3, also rogers et al 1995
    s = lambda basnr, ntr: (2.*basnr/n.sqrt(ntr))**(1/3.)  # see rogers et al. 1995 for factor of 2

    # measure SNR_bl==Q from sigma clipped times with normal mean and std of bispectra. put into time,dm order
    bamean = ba.real.mean(axis=1)
    bastd = ba.real.std(axis=1)

    (meanmin,meanmax) = lib.sigma_clip(bamean)  # remove rfi to estimate noise-like parts
    (stdmin,stdmax) = lib.sigma_clip(bastd)
    clipped = n.where((bamean > meanmin) & (bamean < meanmax) & (bastd > stdmin) & (bastd < stdmax) & (bamean != 0.0))[0]  # remove rfi and zeros
    bameanstd = ba[clipped].real.mean(axis=1).std()
    basnr = bamean/bameanstd    # = S**3/(Q**3 / n.sqrt(n_tr)) = s**3 * n.sqrt(n_tr)

    if Q and verbose:
        print 'Using given Q =', Q
    else:
        Q = ((bameanstd/2.)*n.sqrt(d['ntr']))**(1/3.)
        if verbose:
            print 'Estimating noise per baseline from data. Q (per DM) =', Q

    # detect
    cands = n.where( (bastd/Q**3 < tol*sigbQ3(s(basnr, d['ntr']))) & (basnr > sigma) )  # get compact sources with high snr

    if show or save:
        p.figure()
        ax = p.axes()
        p.subplot(211)
        p.title(str(d['nskip']) + ' nskip, ' + str(len(cands))+' candidates', transform = ax.transAxes)
        p.plot(basnr, 'b.')
        if len(cands[0]) > 0:
            p.plot(cands, basnr[cands], 'r*')
            p.ylim(-2*basnr[cands].max(),2*basnr[cands].max())
        p.xlabel('Integration',fontsize=12,fontweight="bold")
        p.ylabel('SNR_b',fontsize=12,fontweight="bold")
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_position(('outward', 20))
        ax.spines['left'].set_position(('outward', 30))
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        p.subplot(212)
        p.plot(bastd/Q**3, basnr, 'b.')

        # plot reference theory lines
        smax = s(basnr.max(), d['nants'])
        sarr = smax*n.arange(0,101)/100.
        p.plot(sigbQ3(sarr), 1/2.*sarr**3*n.sqrt(d['ntr']), 'k')
        p.plot(tol*sigbQ3(sarr), 1/2.*sarr**3*n.sqrt(d['ntr']), 'k--')
        if len(cands[0]) > 0:
            p.plot(bastd[cands]/Q**3, basnr[cands], 'r*')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_position(('outward', 20))
        ax.spines['left'].set_position(('outward', 30))
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')

        if len(cands[0]) > 0:
            p.axis([0, tol*sigbQ3(s(basnr[cands].max(), d['nants'])), -0.5*basnr[cands].max(), 1.1*basnr[cands].max()])

            # show spectral modulation next to each point
            for candint in cands:
                sm = n.single(round(specmod(data, d, candint),1))
                p.text(bastd[candint]/Q**3, basnr[candint], str(sm), horizontalalignment='right', verticalalignment='bottom')
            p.xlabel('sigma_b/Q^3',fontsize=12,fontweight="bold")
            p.ylabel('SNR_b',fontsize=12,fontweight="bold")
            if save:
                if save == 1:
                    savename = d['filename'].split('.')[:-1]
                    savename.append(str(d['nskip']) + '_bisp.png')
                    savename = string.join(savename,'.')
                elif isinstance(save, types.StringType):
                    savename = save
                print 'Saving file as ', savename
                p.savefig(self.pathout+savename)

    return cands[0], basnr, bastd, Q

def estimate_noiseperbl(data0):
    """ Takes large data array and sigma clips it to find noise per bl for input to detect_bispectra.
    Takes mean across pols and channels for now, as in detect_bispectra.
    """
    
    # define noise per baseline for data seen by detect_bispectra
    data0mean = data0.mean(axis=3).mean(axis=2).imag                      # use imaginary part to estimate noise without calibrated, on-axis signal
    (data0meanmin, data0meanmax) = lib.sigma_clip(data0mean.flatten())
    good = n.where( (data0mean>data0meanmin) & (data0mean<data0meanmax) )
    noiseperbl = data0mean[good].std()   # measure single noise for input to detect_bispectra
    print 'Sigma clip of %.3f to %.3f keeps %d%% of data' % (data0meanmin, data0meanmax, (100.*len(good[0]))/len(data0mean.flatten()))
    print 'Estimate of noise per baseline: %.3f' % noiseperbl
    return noiseperbl

def imagecands(data0, u0, v0, d, cands, clean=False, show=0, save=0):
    """ Image list of candidates returned by detect_bispectra.
    Assumes globals data, u, v w already defined.
    """

    cands2 = []
    for ii in cands:
        print
        print '**Imaging candidate at int_rel=%d, int_abs=%d**' % (ii, d['nskip']+d['itercount']+ii)

        # make image of field
        im = imageint(data0, d, ii, size=d['size'], res=d['res'], clean=clean, show=show)
        imsize = float(d['size'])/d['res']
        pixelscale = 1./d['size']
        peakm, peakl = n.where(im == im.max())
        l0 = d['l0']
        m0 = d['m0']
        m1 = -(imsize/2. - peakm[0])*pixelscale   # indexed starting at top...
        l1 = -(peakl[0] - imsize/2.)*pixelscale   # indexed starting at left.

        # logic on whether imaged cand is good
        centerpix = d['size']/d['res']/2
        peak = n.where(n.max(im) == im)
        print 'Image peak of %e at (%d,%d). Center pixel %d' % (n.max(im), peak[0][0], peak[1][0], centerpix)
        tmp = im[centerpix-centerpix/2:centerpix+centerpix/2, centerpix-centerpix/2:centerpix+centerpix/2]
        im_std = tmp[n.where(tmp <= 0.9*tmp.max())].std()   # estimate std in center half without image peak
        print 'Peak/RMS (central half) = %e' % (im.max()/im_std)

        # phase shift to center for spectral analysis
        lib.phaseshift(data0, d, l1, m1, u0, v0)
        sm = specmod(data0, d, ii)
        lib.phaseshift(data0, d, l0, m0, u0, v0)   # shift back to original location

        print 'Spectral modulation = %.1f. Max expected = %.1f' % (sm, d['specmodfilter']*n.sqrt(d['nchan'])/d['sigma_image'])
        if im.max()/im_std > d['sigma_image']:
            if sm < d['specmodfilter']*n.sqrt(d['nchan'])/d['sigma_image']:
                print 'Imaging shows significant point source. Spectral modulation suggests broad band signal. Saving.'
                cands2.append(ii)
            elif not d['specmodfilter']:    # if no application of filter, fall back on imaging filter
                print 'Image has SNR higher than bisp threshold. Saving'
                cands2.append(ii)
                
    # optionally plot
        if show:
            p.figure()
            p.clf()
            ax = p.axes()

            # text description of candidate
            p.subplot(221, axisbg='white')
            p.title('Candidate @ Tier 2')
            p.xticks([])
            p.yticks([])
            p.text(0.1, 0.8, d['filename'], fontname='sans-serif')
            beamra = n.round(beam, 3)[0]
            beamdec = n.round(beam, 3)[1]
            p.text(0.1, 0.6, 'Beam: ' + str(beamra) + ', ' + str(beamdec) + ', dt: ' + str(d['timescales'][cand[0]]), fontname='sans-serif')
            p.text(0.1, 0.4, 'Integration: ' + str(cand[1]) + ', DM: ' + str(d['dmarr'][cand[2]]), fontname='sans-serif')
            p.text(0.1, 0.2, 'SNR: ' + str(n.round(candsnr, 1)), fontname='sans-serif')

            # image of dispersed visibilities for candidate
            p.subplot(222)
            fov = n.degrees(1./d['res'])*3600.       # field of view in arcseconds
            p.imshow(im, aspect='auto', origin='lower', interpolation='nearest', extent=[fov/2, -fov/2, -fov/2, fov/2])
            p.colorbar()
            p.xlabel('RA/l Offset (arcsec)')
            p.ylabel('Dec/m Offset (arcsec)')

            p.subplot(223)
            dataph = spec(data, 0, len(data))
            p.plot(d['freq'], dataph[ii], '.')
            p.text(0.05, 0.05, 'specmod =' + str(sm), transform = ax.transAxes)
            p.xlabel('Frequency (GHz)')
            p.ylabel('Flux density (Jy)')

            p.subplot(224)
            dataph = n.rot90(dataph)
            sh = dataph.shape
            im = p.imshow(dataph, aspect='auto', origin='upper', interpolation='nearest', extent=(0, sh[1], 0, sh[0]))
            p.colorbar()
#            p.plot(self.obs.dmtrack0[dmind][0], self.obs.dmtrack0[dmind][1],'k.')   # reference dispersed track
            p.xlabel('Integration')
            p.ylabel('Channel')
            if save:
                p.savefig(string.join(d['filename'].split('.')[:-1], '.') + '_sc' + str(d['scan']) + 'sp' + str(d['spw']) + 'i' + str(cand[1]) + '_tier2_.png')

    return cands2

def transient_search(data0, d, triples):
    """ Searches for transient in data of dimensions time, baseline, channel, polarization.
    Gets shared memory data from global namespace. Creates local copy as masked array.
    Not threaded.
    """

#######################
# ** start for loops **
#######################
#    for dtind in range(len(self.timescales)):
#        data0 = time_filter(data0, d, d['timescales'][dtind])
    dtind = 0

    noiseperbl = estimate_noiseperbl(data0)
    print 'Measured noise per baseline (mean pol, mean chan): %.3f.' % noiseperbl

    for dm in d['dmarr']:
        lib.dedisperse(data0, d, dm, verbose=1)        # dedisperses data. data0dm=data0 for now
        print 'DM = %d (max %d)' % (dm, d['dmarr'][-1])
        
        for (ra, dec) in d['delaycenters']:
            if d['searchtype'] == 'bispectrum':
                #            print 'Searching new delay center (%.3f, %.3f)' % (ra, dec)
                lib.phaseshift(data0, d, n.radians(ra), n.radians(dec), u, v)    # phase shifts data in place
                bispectra = lib.make_bispectra(data0, triples)
                bispectra = n.ma.masked_array(bispectra, bispectra==0j)
                if d['datadelay'].max() > 0:
                    bispectra[-d['datadelay'].max():] = n.ma.masked
                (bispcands, basnr, bastd, Q) = detect_bispectra(bispectra, d, sigma=d['sigma_bisp'], Q=noiseperbl, verbose=0)

            # if data calibrated, use image and spectral modulation to filter candidates.
                if d['calibratedfilter'] and len(bispcands) > 0:
                    cands = imagecands(data0, u, v, d, bispcands, show=0, clean=False)
                else:
                    cands = bispcands
            else:
                #                (cands, imsnr, imnoise) = imagesearch2(data0, d, u, v, w, size=d['size'], res=d['res'], clean=False)
                datam = data0[...,0].mean(axis=2)
                (imp,imn) = qimg.imgallstat(u[d['iterint']/2], v[d['iterint']/2], datam, size=d['size'], res=d['res'])
                cands = n.where(imp/imn > d['sigma_image'])
                if len(cands[0]) > 0:
                    for cand in cands[0]:
                        print 'Got one! Int=%d, SNR=%.1f, ' % (cand, (imp/imn)[cand])

            # save candidate info. should be generalized as trigger to dump buffer. maybe post-aggregation in dm-t0 plane?
            # could also image search here! feed results into 'master' stuff below.
            if False:
                #            if len(cands) > 0:
                # get existing master objects to store candidate info
                pkl = open(d['masterfile'], 'rb')
                constructiondict = pickle.load(pkl)
                (masterloc, masterprop, masterdata, masternoise) = pickle.load(pkl)
                pkl.close()
                loclist = masterloc[(ra,dec)]
                proplist = masterprop[(ra,dec)]
                datalist = masterdata[(ra,dec)]
                noiselist = masternoise[(ra,dec)]

                for i in range(len(cands)):
                    intmin = max(0, cands[i] - 32)                                   # candidate lightcuve and spectrogram
                    intmax = min(cands[i] + 32, d['iterint'])
                    dataph = spec(data0, intmin, intmax)
                    loclist.append( [dtind, d['nskip']+d['itercount']+cands[i], d['dmarr'].index(dm)] )    # candidate location (dtind=1 for now)
                    if d['searchtype'] == 'bispectrum':
                        datalist.append( (basnr, bastd, dataph) )
                        proplist.append( [basnr[cands[i]], bastd[cands[i]]] )                  # candidate bisp properties
                        noiselist.append(Q)
                    else:
                        datalist.append( (dataph) )
                        proplist.append( (imsnr[cands[i]]) )
                        noiselist.append( (imnoise[cands[i]]) )

                masterloc[(ra,dec)] = loclist
                masterprop[(ra,dec)] = proplist
                masterdata[(ra,dec)] = datalist
                masternoise[(ra,dec)] = noiselist
                print '%d candidates in beam (%.3f, %.3f) of %s for itercount %d' % (len(loclist), ra, dec, d['filename'], d['itercount'])
                print loclist, proplist

                # save candidate info in pickle file
                if d['savecands']:
                    pkl = open(d['masterfile'], 'wb')
                    pickle.dump(constructiondict, pkl)
                    pickle.dump((masterloc, masterprop, masterdata, masternoise), pkl)
                    pkl.close()

#    phaseshift(data0, d, 0., 0.)   # restore data

def time_filter(data0, d, width, show=0):
        """ Replaces data array with filtered version via convolution in time. Note that this has trouble with zeroed data.
        kernel specifies the convolution kernel. 'm' for mexican hat (a.k.a. ricker, effectively does bg subtraction), 'g' for gaussian. 't' for a tophat. 'b' is a tophat with bg subtraction (or square 'm'). 'w' is a tophat with width that varies with channel, as kept in 'twidth[dmbin]'.
        width is the kernel width with length nchan. should be tuned to expected pulse width in each channel.
        bgwindow is used by 'b' only.
        An alternate design for this method would be to make a new data array after filtering, so this can be repeated for many assumed widths without reading data in anew. That would require more memory, so going with repalcement for now.
        """
        kernel = d['filtershape']
        bgwindow = d['bgwindow']

        if not isinstance(width, types.ListType):
            width = [width] * len(d['chans'])

        # time filter by convolution. functions have different normlizations. m has central peak integral=1 and total is 0. others integrate to 1, so they don't do bg subtraction.
        kernelset = {}  # optionally could make set of kernels. one per width needed. (used only by 'w' for now).

        if kernel == 't':
            print 'Applying tophat time filter.'
            for w in n.unique(width):
                kernel = n.zeros(len(data0))                    # tophat.
                onrange = range(len(kernel)/2 - w/2, len(kernel)/2 + int(ceil(w/2.)))
                kernel[onrange] = 1.
                kernelset[w] = kernel/n.where(kernel>0, kernel, 0).sum()         # normalize to have peak integral=1, thus outside integral=-1.
        elif kernel == 'b':
            print 'Applying tophat time filter with bg subtraction (square mexican hat).'
            for w in n.unique(width):
                kernel = n.zeros(len(data0))                    # tophat.
                onrange = range(len(kernel)/2 - w/2, len(kernel)/2 + int(ceil(w/2.)))
                kernel[onrange] = 1.
                offrange = range(len(kernel)/2 - (bgwindow+w)/2, len(kernel)/2-w/2) + range(len(kernel)/2 + int(ceil(w/2.)), len(kernel)/2 + int(ceil((w+bgwindow)/2.)))
                kernel[offrange] = -1.
                posnorm = n.where(kernel>0, kernel, 0).sum()           # find normalization of positive
                negnorm = n.abs(n.where(kernel<0, kernel, 0).sum())    # find normalization of negative
                kernelset[w] = n.where(kernel>0, kernel/posnorm, kernel/negnorm)    # pos and neg both sum to 1/-1, so total integral=0
        elif kernel == 'g':
            print 'Applying gaussian time filter. Note that effective width is much larger than equivalent tophat width.'
            for w in n.unique(width):
                kernel = signal.gaussian(len(data0), w)     # gaussian. peak not quite at 1 for widths less than 3, so it is later renormalized.
                kernelset[w] = kernel / (w * n.sqrt(2*n.pi))           # normalize to pdf, not peak of 1.
        elif kernel == 'w':
            print 'Applying tophat time filter that varies with channel.'
            for w in n.unique(width):
                kernel = n.zeros(len(data0))                    # tophat.
                onrange = range(len(kernel)/2 - w/2, len(kernel)/2 + int(ceil(w/2.)))
                kernel[onrange] = 1.
                kernelset[w] = kernel/n.where(kernel>0, kernel, 0).sum()         # normalize to have peak integral=1, thus outside integral=-1.
        elif kernel == None:
            print 'Applying no time filter.'
            return data0

        if show:
            for kernel in kernelset.values():
                p.plot(kernel,'.')
            p.title('Time filter kernel')
            p.show()

        # take ffts (in time)
        datafft = n.fft.fft(data0, axis=0)
        kernelsetfft = {}
        for w in n.unique(width):
            kernelsetfft[w] = n.fft.fft(n.roll(kernelset[w], len(data0)/2))   # seemingly need to shift kernel to have peak centered near first bin if convolving complex array (but not for real array?)

        # filter by product in fourier space
        for i in range(d['nbl']):    # **can't find matrix product I need, so iterating over nbl, chans, npol**
            for j in range(len(d['chans'])):
                for k in range(d['npol']):
                    datafft[:,i,j,k] = datafft[:,i,j,k]*kernelsetfft[width[j]]    # index fft kernel by twidth

        # ifft to restore time series
#        return n.ma.masked_array(n.fft.ifft(datafft, axis=0), self.flags[:self.nints,:, self.chans,:] == 0)
        return n.fft.ifft(datafft, axis=0)

def time_filter2(data0, d, width, show=0):
    data0c = n.zeros_like(data0)
    for i in range(1,len(data0c)-1):
        data0c[i] = data0[i] - (data0[i-1]+data0[i+1])/2
    return data0c

def specmod(data0, d, ii):
    """Calculate spectral modulation for given track.
    Spectral modulation is basically the standard deviation of a spectrum. 
    This helps quantify whether the flux is located in a narrow number of channels or across all channels.
    Broadband signal has small modulation (<sqrt(nchan)/SNR) while RFI has larger values.
    See Spitler et al 2012 for details.
    """

    bfspec = data0[ii].mean(axis=0).mean(axis=1).real   # mean over bl and pol
    sm = n.sqrt( ((bfspec**2).mean() - bfspec.mean()**2) / bfspec.mean()**2 )

    return sm

def imageint(data0, d, i, mode='split', size=48000, res=500, clean=True, gain=0.01, tol=1e-4, newbeam=0, save=0, show=0):
    """ Use apiy to image integration
    mode defines how frequency dependence is handled. 'split' means separate uv and data points in frequency (but not mfs). 'mean' means mean vis across frequency.
    int is used to select uvw coordinates for track. default is first int.
    pol can be 'i' for a Stokes I image (mean over pol dimension) or a pol index.
    default params size and res are good for 1.4 GHz VLA, C-config image.
    clean determines if image is cleaned and beam corrected. gain/tol are cleaning params.
    newbeam forces the calculation of a new beam for restoring the cleaned image.
    save=0 is no saving, save=1 is save with default name, save=<string>.png uses custom name (must include .png). 
        """

    td = data0[i].mean(axis=2)     # take single integration and mean over pol axis

# apply w phase rotation. generally this is done externally (e.g., by data writing software) and is not needed here.
#        wrot = lambda w: n.exp(-2j*n.pi*n.outer(w, self.freq/self.freq_orig[0]))
#        td = td*wrot(self.w[i])

    # define handling of freq axis
    if mode == 'split':
        td = td.flatten()
        uu = n.outer(u[i], d['freq']/d['freq_orig'][0]).flatten()
        vv = n.outer(v[i], d['freq']/d['freq_orig'][0]).flatten()
        ww = n.outer(w[i], d['freq']/d['freq_orig'][0]).flatten()
#        uu = n.outer(d['u0'][i], d['freq']/d['freq_orig'][0]).flatten()
#        vv = n.outer(d['v0'][i], d['freq']/d['freq_orig'][0]).flatten()
#        ww = n.outer(d['w0'][i], d['freq']/d['freq_orig'][0]).flatten()
    elif mode == 'mean':
        td = td.mean(axis=1)
        uu = u[i]
        vv = v[i]
        ww = w[i]
#        uu = d['u0'][i]
#        vv = d['v0'][i]
#        ww = d['w0'][i]
    else:
        print 'Mode must be \'mean\' or \'split\'.'
        return 0

    fov = n.degrees(1./res)*3600.  # field of view in arcseconds

    # make image
    ai = aipy.img.Img(size=size, res=res)
    uvw_new, td_new = ai.append_hermitian( (uu, vv, ww), td)
    ai.put(uvw_new, td_new)
    image = ai.image(center = (size/res/2, size/res/2))
    image_final = image

    # optionally clean image
    if clean:
        print 'Cleaning image...'
        beam = ai.bm_image()
        beamgain = aipy.img.beam_gain(beam[0])
        (cleanim, dd) = aipy.deconv.clean(image, beam[0], verbose=True, gain=gain, tol=tol)

        try:
            import gaussfitter
            print 'Restoring image with fit to beam shape...'
            beam_centered = ai.bm_image(center=(size/res/2, size/res/2))
            peak = n.where(beam_centered[0] >= 0.1*beam_centered[0].max(), beam_centered[0], 0.)
            beam_params = gaussfitter.gaussfit(peak)
            kernel = n.roll(n.roll(gaussfitter.twodgaussian(beam_params, shape=n.shape(beam[0])), size/res/2, axis=0), size/res/2, axis=1)   # fit to beam at center, then roll to corners for later convolution step
        except ImportError:
            print 'Restoring image with peak of beam...'
            kernel = n.where(beam[0] >= 0.4*beam[0].max(), beam[0], 0.)  # take only peak (gaussian part) pixels of beam image

        restored = aipy.img.convolve2d(cleanim, kernel)
        image_restored = (restored + dd['res']).real/beamgain
        image_final = image_restored

    if show or save:
        p.figure()
        ax = p.axes()
        ax.set_position([0.2,0.2,0.7,0.7])
#        im = p.imshow(image_final, aspect='auto', origin='upper', interpolation='nearest', extent=[-fov/2, fov/2, -fov/2, fov/2])
        im = p.imshow(image_final, aspect='auto', origin='lower', interpolation='nearest', extent=[fov/2, -fov/2, -fov/2, fov/2])
        cb = p.colorbar(im)
        cb.set_label('Flux Density (Jy)',fontsize=12,fontweight="bold")
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_position(('outward', 20))
        ax.spines['left'].set_position(('outward', 30))
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        p.xlabel('RA/l Offset (arcsec)',fontsize=12,fontweight="bold")
        p.ylabel('Dec/m Offset (arcsec)',fontsize=12,fontweight="bold")

        centerpix = size/res/2
        peak = n.where(n.max(image_final) == image_final)
        print 'Image peak of %e at (%d,%d)' % (n.max(image_final), peak[0][0], peak[1][0])
        print 'Center pixel coord: %d' % (centerpix)
        tmp = image_final[centerpix-centerpix/2:centerpix+centerpix/2, centerpix-centerpix/2:centerpix+centerpix/2]
        image_std = tmp[n.where(tmp <= 0.9*tmp.max())].std()   # estimate std in center half without image peak
#        print 'Peak/RMS = %e' % (tmpimage.max()/image_std)
        print 'Peak/RMS (central half) = %e' % (image_final.max()/image_std)

        if save:
            if save == 1:
                savename = d['filename'].split('.')[:-1]
                savename.append(str(d['nskip']) + '_im.png')
                savename = string.join(savename,'.')
            elif isinstance(save, string):
                savename = save
            print 'Saving file as ', savename
            p.savefig(savename)

    return image_final

def imagesearch(data0, d, u, v, w, mode='split', size=48000, res=500, clean=True, gain=0.01, tol=1e-4, newbeam=0, save=0, show=0):
    """ Use apiy to image integration
    mode defines how frequency dependence is handled. 'split' means separate uv and data points in frequency (but not mfs). 'mean' means mean vis across frequency.
    int is used to select uvw coordinates for track. default is first int.
    pol can be 'i' for a Stokes I image (mean over pol dimension) or a pol index.
    default params size and res are good for 1.4 GHz VLA, C-config image.
    clean determines if image is cleaned and beam corrected. gain/tol are cleaning params.
    newbeam forces the calculation of a new beam for restoring the cleaned image.
    save=0 is no saving, save=1 is save with default name, save=<string>.png uses custom name (must include .png). 
        """

    data0m = data0.mean(axis=3)

    cands = []
    imsnr = []
    imnoise = []
    for i in range(len(data0m)):
        print 'Imaging int %d' % i

        td = data0m[i]

        # define handling of freq axis
        if mode == 'split':
            td = td.flatten()
            uu = n.outer(u[i], d['freq']/d['freq_orig'][0]).flatten()
            vv = n.outer(v[i], d['freq']/d['freq_orig'][0]).flatten()
            ww = n.outer(w[i], d['freq']/d['freq_orig'][0]).flatten()
        elif mode == 'mean':
            td = td.mean(axis=1)
            uu = u[i]
            vv = v[i]
            ww = w[i]
        else:
            print 'Mode must be \'mean\' or \'split\'.'
            return 0

        fov = n.degrees(1./res)*3600.  # field of view in arcseconds

        # make image
        useaipy = False
        if useaipy:
            ai = aipy.img.Img(size=size, res=res)
            uvw_new, td_new = ai.append_hermitian( (uu, vv, ww), td)
            ai.put(uvw_new, td_new)
            image = ai.image(center = (size/res/2, size/res/2))
        else:
            image = qimg.imgone(uu, vv, ww, td, size, res)
        im = image

        # optionally clean image
        if clean:
            print 'Cleaning image...'
            beam = ai.bm_image()
            beamgain = aipy.img.beam_gain(beam[0])
            (cleanim, dd) = aipy.deconv.clean(image, beam[0], verbose=True, gain=gain, tol=tol)

            try:
                import gaussfitter
                print 'Restoring image with fit to beam shape...'
                beam_centered = ai.bm_image(center=(size/res/2, size/res/2))
                peak = n.where(beam_centered[0] >= 0.1*beam_centered[0].max(), beam_centered[0], 0.)
                beam_params = gaussfitter.gaussfit(peak)
                kernel = n.roll(n.roll(gaussfitter.twodgaussian(beam_params, shape=n.shape(beam[0])), size/res/2, axis=0), size/res/2, axis=1)   # fit to beam at center, then roll to corners for later convolution step
            except ImportError:
                print 'Restoring image with peak of beam...'
                kernel = n.where(beam[0] >= 0.4*beam[0].max(), beam[0], 0.)  # take only peak (gaussian part) pixels of beam image

            restored = aipy.img.convolve2d(cleanim, kernel)
            image_restored = (restored + dd['res']).real/beamgain
            im = image_restored

        # logic on whether imaged cand is good
        centerpix = d['size']/d['res']/2
        tmp = im[centerpix-centerpix/2:centerpix+centerpix/2, centerpix-centerpix/2:centerpix+centerpix/2]
        im_std = tmp[n.where(tmp <= 0.9*tmp.max())].std()   # estimate std in center half without image peak

        imsnr.append(im.max()/im_std)
        imnoise.append(im_std)

        if im.max()/im_std > d['sigma_image']:
            peakm, peakl = n.where(im == im.max())
            print 'Image peak of %e at (%d,%d). Center pixel %d' % (n.max(im), peakm[0], peakl[0], centerpix)
            print 'Peak/RMS (central half) = %e' % (im.max()/im_std)

            # phase to peak
            imsize = float(d['size'])/d['res']
            pixelscale = 1./d['size']
            l0 = d['l0']
            m0 = d['m0']
            m1 = -(imsize/2. - peakm[0])*pixelscale   # indexed starting at top...
            l1 = -(peakl[0] - imsize/2.)*pixelscale   # indexed starting at left.
            lib.phaseshift(data0, d, l1, m1, u, v)
            sm = specmod(data0, d, i)
            lib.phaseshift(data0, d, l0, m0, u, v)
            print 'Spectral modulation = %.1f. Max expected = %.1f' % (sm, d['specmodfilter']*n.sqrt(d['nchan'])/d['sigma_image'])

            if sm < d['specmodfilter']*n.sqrt(d['nchan'])/d['sigma_image']:
                print 'Imaging shows significant point source. Spectral modulation suggests broad band signal. Saving.'
                cands.append(i)
            elif not d['specmodfilter']:    # if no application of filter, fall back on imaging filter
                print 'Image has high SNR.'
                cands.append(i)
    return cands, imsnr, imnoise

def imagesearch2(data0, d, u, v, w, size=48000, res=500, clean=True, gain=0.01, tol=1e-4, newbeam=0, save=0, show=0):
    """ Use qimg to image all integrations simultaneously with single uvw gridding.
    pol can be 'i' for a Stokes I image (mean over pol dimension) or a pol index.
    default params size and res are good for 1.4 GHz VLA, C-config image.
        """

    cands = []
    imsnr = []
    imnoise = []
    
    # expand uvw to include frequency axis
    midint = len(u)/2
    uu = n.outer(u[midint], d['freq']/d['freq_orig'][0])
    vv = n.outer(v[midint], d['freq']/d['freq_orig'][0])
    ww = n.outer(w[midint], d['freq']/d['freq_orig'][0])
    
    fov = n.degrees(1./res)*3600.  # field of view in arcseconds
    sh = data0.shape
    uuu = n.concatenate([uu,-uu], axis=0).reshape(2*sh[1]*sh[2])
    vvv = n.concatenate([vv,-vv], axis=0).reshape(2*sh[1]*sh[2])
    www = n.concatenate([ww,-ww], axis=0).reshape(2*sh[1]*sh[2])
    data0c = n.concatenate([data0, data0.conjugate()], axis=0).reshape(sh[0],2*sh[1]*sh[2],sh[3])
    images = qimg.imgall(uuu, vvv, www, data0c, size, res)
    
    # logic on whether imaged cand is good
    centerpix = d['size']/d['res']/2
    for i in xrange(len(images)):
        im = images[i]
        tmp = im[centerpix-centerpix/2:centerpix+centerpix/2, centerpix-centerpix/2:centerpix+centerpix/2]
        im_std = tmp[n.where(tmp <= 0.9*tmp.max())].std()   # estimate std in center half without image peak

        imsnr.append(im.max()/im_std)
        imnoise.append(im_std)

        if im.max()/im_std > d['sigma_image']:
            peakm, peakl = n.where(im == im.max())
            print 'Image peak of %e at (%d,%d). Center pixel %d' % (n.max(im), peakm[0], peakl[0], centerpix)
            print 'Peak/RMS (central half) = %e' % (im.max()/im_std)

            # phase to peak
            imsize = float(d['size'])/d['res']
            pixelscale = 1./d['size']
            l0 = d['l0']
            m0 = d['m0']
            m1 = -(imsize/2. - peakm[0])*pixelscale   # indexed starting at top...
            l1 = -(peakl[0] - imsize/2.)*pixelscale   # indexed starting at left.
            lib.phaseshift(data0, d, l1, m1, u, v)
            sm = specmod(data0, d, i)
            lib.phaseshift(data0, d, l0, m0, u, v)
            print 'Spectral modulation = %.1f. Max expected = %.1f' % (sm, d['specmodfilter']*n.sqrt(d['nchan'])/d['sigma_image'])

            if sm < d['specmodfilter']*n.sqrt(d['nchan'])/d['sigma_image']:
                print 'Imaging shows significant point source. Spectral modulation suggests broad band signal. Saving.'
                cands.append(i)
            elif not d['specmodfilter']:    # if no application of filter, fall back on imaging filter
                print 'Image has high SNR.'
                cands.append(i)
    return cands, imsnr, imnoise

def find_candidates(pklfile, d_neighbor=2, island_size_min=1, island_snr_min=5., show=0, save=0):
    """ Visualize and define candidates to follow up from pickle file with search state and candidate dictionary.
    d_neighbor is the dm-time distance over which a neighbor is defined.
    island_size_min is the minimum number of detections requried to call an island a candidate
    island_snr_ min is the minimum snr requried to call an island a candidate
    returns the SNR peak of each island identified with d_neighbor.
    """

    # read in pickle file of candidates
    pkl = open(pklfile, 'rb')
    d = pickle.load(pkl)
    (masterloc, masterprop, masterdata, masternoise) = pickle.load(pkl)
    pkl.close()

    neighbors = {}
    for (ra,dec) in d['delaycenters']:
        neighbors[(ra,dec)] = []

    dmarr = n.array(d['dmarr'])

    # measure number of neighbors
    for beam in masterloc.keys():
        neighbors_lm = []
        cands_lm = masterloc[beam]
        for cand in cands_lm:
            nn = -1    # to avoid self-counting
            diff = n.array(cands_lm) - n.array(cand)
            for dd in diff:
                if ( (n.abs(dd[0]) <= d_neighbor) & (n.abs(dd[1]) <= d_neighbor) & (n.abs(dd[2]) <= d_neighbor) ):  # if dmbin and int are within d_neighbor, then they are neighbors
                    nn = nn+1
            neighbors_lm.append(nn)
        neighbors[beam] = neighbors_lm

    # define islands in DM-time
    islands = {}
    for (ra,dec) in masterloc.keys():
        islands[(ra,dec)] = []
    for beam in masterloc.keys():
        cands_lm = masterloc[beam]
        if len(cands_lm) == 0: continue
        candlist = list(cands_lm)
        islands_lm = [[candlist.pop()]]
        icount = 0
        islands_lm[icount] = n.vstack( (islands_lm[icount][0],islands_lm[icount][0]) )
        while ( (icount < len(islands_lm)) & (len(candlist) > 0) ):
            ccount = 0
            while ccount < len(candlist):
                ww = n.where( (n.abs(n.array(candlist[ccount]) - n.array(islands_lm[icount]))[:, 0] <= d_neighbor) & (n.abs(n.array(candlist[ccount]) - n.array(islands_lm[icount]))[:, 1] <= d_neighbor) & (n.abs(n.array(candlist[ccount]) - n.array(islands_lm[icount]))[:, 2] <= d_neighbor) & ((n.abs(n.array(candlist[ccount]) - n.array(islands_lm[icount]))[:, 0] > 0) | (n.abs(n.array(candlist[ccount]) - n.array(islands_lm[icount]))[:, 1] > 0) | (n.abs(n.array(candlist[ccount]) - n.array(islands_lm[icount]))[:, 2] > 0)) )[0]  # coords within 1, but not both equal to zero (i.e., same candidate)
                if len(ww) > 0:
                    newcand = candlist.pop(ccount)
                    islands_lm[icount] = n.vstack( (islands_lm[icount], newcand) )
                    ccount = 0
                else:
                    ccount = ccount + 1
            if len(candlist) == 0: break
            newcand = candlist.pop()
            islands_lm.append(newcand)
            icount = icount + 1
            islands_lm[icount] = n.vstack( (islands_lm[icount],islands_lm[icount]) )

        # trim off superfluous element in each island
        for i in range(len(islands_lm)):
            islands_lm[i] = islands_lm[i][1:]

        islands[beam] = islands_lm
    print 'Identified DM-time islands.'

    # find snr peaks of islands including filters for island size and snr
    islandmaxloc = {}
    islandmaxsnr = {}
    for beam in islands.keys():
        islandmaxlocl = []; islandmaxsnrl = []
        islands_lm = islands[beam]
        for island in islands_lm:
            if len(island) >= island_size_min:
                islandind = []
                for i in range(len(island)):
                    islandind.append(n.where( (n.array(masterloc[beam])[:,0]==island[i,0]) & (n.array(masterloc[beam])[:,1]==island[i,1]) & (n.array(masterloc[beam])[:,2]==island[i,2]) )[0][0])
                maxsnrind = n.where(n.array(masterprop[beam])[islandind,0] == n.max(n.array(masterprop[beam])[islandind,0]))
                if n.array(masterprop[beam])[islandind][maxsnrind][0][0] > island_snr_min:
                    islandmaxlocl.append(n.array(masterloc[beam])[islandind][maxsnrind][0])
                    islandmaxsnrl.append(n.array(masterprop[beam])[islandind][maxsnrind][0][0])
        islandmaxloc[beam] = n.array(islandmaxlocl)
        islandmaxsnr[beam] = n.array(islandmaxsnrl)

    if show or save:
        cm = p.get_cmap('gist_rainbow')
        p.figure(1, figsize=(12,9))
        tl1 = p.subplot(221)
        tr1 = p.subplot(222)
        bl1 = p.subplot(223)
        br1 = p.subplot(224)

        beamind = 0
        for beam in masterloc.keys():  # iterate over beam candidate groups
            loc = n.array(masterloc[beam])
            prop = n.array(masterprop[beam])
            shiftind = float(beamind)/len(masterloc.keys())   # shift point location (and color) per beam for clarity
            if len(loc) == 0: break
            p.figure(1)
            p.subplot(tl1)
            p.scatter(d['inttime']*loc[:,1], dmarr[loc[:,2]] + 0.5*(dmarr[1]-dmarr[0])*shiftind, s=8*prop[:,0], facecolor='none', color=cm(shiftind), alpha=0.5, clip_on=False)
            for j in range(len(loc)):   # iterate over cands in group to plot one point each
                if neighbors[beam][j] > 1:
                    p.text(d['inttime']*loc[j,1], dmarr[loc[j,2]] + 0.5*(dmarr[1]-dmarr[0])*shiftind, s=str(neighbors[beam][j]), horizontalalignment='center', verticalalignment='center', fontsize=9, color=cm(shiftind), alpha=0.5)
                        
            # plot island peaks
            p.scatter(d['inttime']*islandmaxloc[beam][:,1], dmarr[islandmaxloc[beam][:,2]] + 0.5*(dmarr[1]-dmarr[0])*shiftind, s=30*islandmaxsnr[beam], color=cm(shiftind), alpha=0.8, marker='+', clip_on=False)
            p.xticks([])
            axis_tl1 = p.axis()
            
            p.subplot(tr1)
            dms = d['dmarr']
            dms.append(dms[-1]+(dms[-1]-dms[-2]))
            dmhist = n.array(dms) - (dms[-1]-dms[-2])/2.
            p.hist(dmarr[loc[:,2]], orientation='horizontal', color=cm(shiftind), histtype='step', alpha=0.5, bins=dmhist)
            p.yticks([])
            axis_tr1 = p.axis()
            p.axis([0, 1.1*axis_tr1[1], axis_tl1[2], axis_tl1[3]])

            p.subplot(bl1)
            p.scatter(d['inttime']*loc[:,1], prop[:,0], color=cm(shiftind), alpha=0.3, clip_on=False)
            p.scatter(d['inttime']*islandmaxloc[beam][:,1], islandmaxsnr[beam], color=cm(shiftind), marker='+', s=100, alpha=0.8, clip_on=False)
            axis_bl1 = p.axis()
            axis_bl1 = p.axis([axis_tl1[0], axis_tl1[1], axis_bl1[2], axis_bl1[3]])

            p.subplot(br1)
            p.hist(prop[:,0], orientation='horizontal', color=cm(shiftind), histtype='step', alpha=0.5, bins=len(prop)/2)
            p.yticks([])
            axis_br1 = p.axis()
            p.axis([0, 1.1*axis_br1[1], axis_bl1[2], axis_bl1[3]])

            beamind += 1

        p.figure(1)
#        p.title('Search Summary Plot')
        p.subplot(tl1)
#        p.axis((0,tarr[-1],-0.5,dmarr[-1]+0.5*(dmarr[1]-dmarr[0])))
#        p.xlabel('Time (s)',fontsize=12,fontweight="bold")
        p.ylabel('DM (pc/cm3)',fontsize=14,fontweight="bold")
        p.subplot(bl1)
        p.xlabel('Rel. Time (s)',fontsize=14,fontweight="bold")
        p.ylabel('SNR',fontsize=14,fontweight="bold")
        p.subplot(br1)
#        p.ylabel('SNR',fontsize=12,fontweight="bold")
        p.xlabel('Count',fontsize=14,fontweight="bold")
        p.show()

    return islands, islandmaxloc    # return peaks of islands

def candplot1(pklfile, islands, plotcands, save=0):
    """ Uses results of bispectrum (saved in pklfile and organized by find_candidates) to produce summary plot per candidate.
    This is a "tier 1" plot, since it only uses info from bispectrum search algorithm and saved in real time.
    Requires islands and plotcands are made by find_candidates.
    Candidate should be dictionary with keys of beam location and value of a list of 3-element lists (dtind, int, dmbin).
    save defines whether plots are also saved.
    """

    print 'Building Tier 1 candidate plot.'

    # read in pickle file of candidates
    pkl = open(pklfile, 'rb')
    d = pickle.load(pkl)
    (masterloc, masterprop, masterdata, masternoise) = pickle.load(pkl)
    pkl.close()

    if n.array([len(plotcands[beam]) for beam in plotcands.keys()]).sum() == 0:
        print 'No candidates available...'
        return

    for beam in plotcands.keys():
        if len(plotcands[beam]) == 0:
            continue
        if (n.ndim(plotcands[beam]) != 2) | (n.shape(plotcands[beam])[1] != 3):
            print 'Candidate does not seem to be in proper format. Skipping.'
            continue

        for cand in plotcands[beam]:
            # find island containing candidate
            for island in islands[beam]:
                if n.any([n.all(member == cand) for member in island]):   # if any member of island that has all three coords the same as plotcand, break
                    break

            # next find snr for all candidates in island
            islandind = []
            for i in range(len(island)):
                ww = n.where( (n.array(masterloc[beam])[:,0]==island[i,0]) & (n.array(masterloc[beam])[:,1]==island[i,1]) & (n.array(masterloc[beam])[:,2]==island[i,2]) )
                islandind.append(ww[0][0])
            loc = n.array(masterloc[beam])[islandind]
            snr = n.array(masterprop[beam])[islandind, 0]

            # then plot snr as function of dmind and dtind for island
            fixed_dt = n.where(loc[:,0] == cand[0])
            fixed_dm = n.where(loc[:,2] == cand[2])
            dmdist = n.squeeze(n.array(d['dmarr'])[loc[fixed_dt, 2]])   # seems to have superfluous axis...?
            dmsnr = snr[fixed_dt]
            dtdist = n.squeeze(n.array(d['timescales'])[loc[fixed_dm][:, 0]])   # seems to have superfluous axis...?
            dtsnr = snr[fixed_dm]

            # find master index of candidate
            ww = n.where([n.all(mastercand == cand) for mastercand in masterloc[beam]])[0][0]

            # plot candidate info
            p.figure()
            p.clf()

            # text description of candidate
            p.subplot(221, axisbg='white')
            p.title('Candidate @ Tier 1')
            p.xticks([])
            p.yticks([])
            p.text(0.1, 0.8, d['filename'], fontname='sans-serif')
            beamra = n.round(beam, 3)[0]
            beamdec = n.round(beam, 3)[1]
            p.text(0.1, 0.6, 'Beam: ' + str(beamra) + ', ' + str(beamdec) + ', dt: ' + str(d['timescales'][cand[0]]), fontname='sans-serif')
            p.text(0.1, 0.4, 'Integration: ' + str(cand[1]) + ', DM: ' + str(d['dmarr'][cand[2]]), fontname='sans-serif')
            p.text(0.1, 0.2, 'SNR: ' + str(n.round(masterprop[beam][ww][0], 1)), fontname='sans-serif')

            p.subplot(222)
            p.plot(dmdist, dmsnr, 'b.', clip_on=False, label='DM at peak dt')
            p.xlabel('DM (pc/cm3)')
            p.ylabel('SNR')
            p.twiny()
            p.plot(dtdist, dtsnr, 'r+', clip_on=False, label='dt at peak DM')
            p.xlabel('dt (ints)')
            p.subplot(223)
            p.plot(masterdata[beam][ww][0], 'b.', label='Mean B')
            p.xlabel('Integration')
            p.ylabel('Mean, Std of Bispectra')
            p.twinx()
            p.plot(masterdata[beam][ww][1], 'r.', label='Std B')
            p.subplot(224)
            dataph = n.rot90(masterdata[beam][ww][2])
            sh = dataph.shape
            im = p.imshow(dataph, aspect='auto', origin='upper', interpolation='nearest', extent=(0, sh[1], 0, sh[0]))
            p.colorbar()
#            p.plot(self.obs.dmtrack0[cand[2]][0], self.obs.dmtrack0[cand[2]][1],'k.')  # reference dispersed track
            p.xlabel('Integration')
            p.ylabel('Channel')
            if save:
                p.savefig(string.join(d['filename'].split('.')[:-1], '.') + '_sc' + str(d['scan']) + 'sp' + str(d['spw']) + 'i' + str(cand[1]) + '_tier1_.png')

def candplot2(pklfile, plotcands, twindow=100, show=0, save=0, **kargs):
    """ Builds summary plots from scratch, including phased beam and imaging for each plot candidate.
    This is a "tier 2" plot, since operates in follow-up mode and uses all interferometric info.
    Candidate should be dictionary with keys of beam location and value of a list of 3-element lists (dt, int, dmbin).
    twindow defines the number of integrations to read, centered on the candidate.
    Use kargs to define construction dictionary values to override.
    """

    global data, flag, u, v, w, time

    print 'Building Tier 2 candidate plot.'

    # read in pickle file of search state
    pkl = open(pklfile, 'rb')
    d = pickle.load(pkl)
    pkl.close()

    # override saved constructiondict items
    if len(kargs) > 0:
        for key in kargs:
            print 'Overriding %s from %s to %s' % (key, d[key], kargs[key])
            d[key] = kargs[key]

    if n.array([len(plotcands[beam]) for beam in plotcands.keys()]).sum() == 0:
        print 'No candidates available...'
        return

    # set up data arrays
    data_mem = mp.Array(ctypes.c_double, twindow*d['nbl']*d['nchan']*d['npol']*2)    # x2 to store complex values in double array
    flag_mem = mp.Array(ctypes.c_bool, twindow*d['nbl']*d['nchan']*d['npol'])
    u_mem = mp.Array(ctypes.c_double, twindow*d['nbl'])
    v_mem = mp.Array(ctypes.c_double, twindow*d['nbl'])
    w_mem = mp.Array(ctypes.c_double, twindow*d['nbl'])
    time_mem = mp.Array(ctypes.c_double, twindow)
    data = numpyview(data_mem, 'complex', (twindow, d['nbl'], d['nchan'], d['npol']))
    flag = numpyview(flag_mem, 'bool', (twindow, d['nbl'], d['nchan'], d['npol']))
    u = numpyview(u_mem, 'double', (twindow, d['nbl']))
    v = numpyview(v_mem, 'double', (twindow, d['nbl']))
    w = numpyview(w_mem, 'double', (twindow, d['nbl']))
    time = numpyview(time_mem, 'double', (twindow))

    for beam in plotcands.keys():
        if len(plotcands[beam]) == 0:
            continue
        if (n.ndim(plotcands[beam]) != 2) | (n.shape(plotcands[beam])[1] != 3):
            print 'Candidate does not seem to be in proper format. Skipping.'
            continue

        (ra, dec) = beam
        print 'Generating plots for candidate with delaycenter: (%s, %s)' % (str(ra), str(dec))
        for cand in plotcands[beam]:
            (dtind, ii, dmind) = cand
            print 'Candidate at (dtind, ii, dmind) = (%d, %d, %d)' % (dtind, ii, dmind)

            # set up data for single iteration of size for display. not threaded
            d['nskip'] = ii - twindow/2       # read starting a bit before candidate
            d['iterint'] = twindow + 5             # read a bit extra to trim off later
            d['nints'] = twindow + 5              # read a bit extra to trim off later

            # another chance to override saved constructiondict items
            if len(kargs) > 0:
                for key in kargs:
                    print 'Overriding %s from %s to %s' % (key, d[key], kargs[key])
                    d[key] = kargs[key]

            readprep(d)
            data0, flag0, u0, v0, w0, time0 = readiter(d)
            u[:] = u0[:twindow]; v[:] = v0[:twindow]; w[:] = w0[:twindow]; time[:] = time0[:twindow];
            data[:] = n.ma.masked_array(data0[:twindow], flag0[:twindow] == 1)   # mask of True for flagged data (flags=True in MS)
            ms.iterend()         # tidy up
            ms.close()
            triples = lib.make_triples(d)
            d['ntr'] = len(triples)
            d['datadelay'] = lib.calc_delay(d, 0)   

            # process!
            noiseperbl = estimate_noiseperbl(data)
            lib.dedisperse(data, d, d['dmarr'][dmind], verbose=1)        # dedisperses 'data' in place
            lib.phaseshift(data, d, n.radians(ra), n.radians(dec), u, v)            # shift phase center to candidate beam
#            time_filter(self.timescales[dtind], self.filtershape, bgwindow=self.bgwindow)    # tophat filter (bgwindow not used)
#            noiseperbl = self.obs.data.mean(axis=3).mean(axis=2).real.std()   # measure single noise for input to detect_bispectra

            # reproduce bispectrum results
            bispectra = lib.make_bispectra(data, triples)
            bispectra = n.ma.masked_array(bispectra, bispectra == 0j)   # set flags for zero data
            if d['datadelay'].max() > 0:
                bispectra[-d['datadelay'].max():] = n.ma.masked
            (cands, basnr, bastd, Q) = detect_bispectra(bispectra, d, sigma=d['sigma_bisp'], Q=noiseperbl, verbose=1, show=show)
            print 'Reproduced bispectrum candidate, ', [cands[i] + d['nskip'] for i in range(len(cands))], basnr[cands], bastd[cands]
            candsnr = basnr[cands]

            # make image of field
            im = imageint(data, d, twindow/2, size=d['size'], res=d['res'], clean=0, show=show)
            imsize = float(d['size'])/d['res']
            pixelscale = 1./d['size']
            peakm, peakl = n.where(im == im.max())
            m1 = -(imsize/2. - peakm[0])*pixelscale   # indexed starting at top...
            l1 = -(peakl[0] - imsize/2.)*pixelscale   # indexed starting at left.

            # phase shift to center for spectral analysis
            lib.phaseshift(data, d, l1, m1, u, v)
            sm = n.round(specmod(data, d, twindow/2), 2)

            # now plot
            p.figure()
            p.clf()
            ax = p.axes()

            # text description of candidate
            p.subplot(221, axisbg='white')
            p.title('Candidate @ Tier 2')
            p.xticks([])
            p.yticks([])
            p.text(0.1, 0.8, d['filename'], fontname='sans-serif')
            beamra = n.round(beam, 3)[0]
            beamdec = n.round(beam, 3)[1]
            p.text(0.1, 0.6, 'Beam: ' + str(beamra) + ', ' + str(beamdec) + ', dt: ' + str(d['timescales'][cand[0]]), fontname='sans-serif')
            p.text(0.1, 0.4, 'Integration: ' + str(cand[1]) + ', DM: ' + str(d['dmarr'][cand[2]]), fontname='sans-serif')
            p.text(0.1, 0.2, 'SNR: ' + str(n.round(candsnr, 1)), fontname='sans-serif')

            # image of dispersed visibilities for candidate
            p.subplot(222)
            fov = n.degrees(1./d['res'])*3600.       # field of view in arcseconds
            p.imshow(im, aspect='auto', origin='lower', interpolation='nearest', extent=[fov/2, -fov/2, -fov/2, fov/2])
            p.colorbar()
            p.xlabel('RA/l Offset (arcsec)')
            p.ylabel('Dec/m Offset (arcsec)')

            p.subplot(223)
            dataph = spec(data, 0, len(data))
            p.plot(d['freq'], dataph[twindow/2], '.')
            p.text(0.05, 0.05, 'specmod =' + str(sm), transform = ax.transAxes)
            p.xlabel('Frequency (GHz)')
            p.ylabel('Flux density (Jy)')

            p.subplot(224)
            dataph = n.rot90(dataph)
            sh = dataph.shape
            im = p.imshow(dataph, aspect='auto', origin='upper', interpolation='nearest', extent=(0, sh[1], 0, sh[0]))
            p.colorbar()
#            p.plot(self.obs.dmtrack0[dmind][0], self.obs.dmtrack0[dmind][1],'k.')   # reference dispersed track
            p.xlabel('Integration')
            p.ylabel('Channel')
            if save:
                p.savefig(string.join(d['filename'].split('.')[:-1], '.') + '_sc' + str(d['scan']) + 'sp' + str(d['spw']) + 'i' + str(cand[1]) + '_tier2_.png')

def readprep(d):
    """ Prepare to read data
    """

    filename = d['filename']; spw = d['spw']; iterint = d['iterint']; datacol = d['datacol']; selectpol = d['selectpol']
    scan = d['scan']; d['npol'] = len(d['selectpol']); nints = d['nints']; nskip = d['nskip']

    # read metadata either from pickle or ms file
    pklname = string.join(filename.split('.')[:-1], '.') + '_init.pkl'
    if os.path.exists(pklname):
        print 'Pickle of initializing info found. Loading...'
        pkl = open(pklname, 'r')
        try:
            (d['npol_orig'], d['nbl'], d['blarr'], d['inttime'], spwinfo, scansummary) = pickle.load(pkl)
        except EOFError:
            print 'Bad pickle file. Exiting...'
            return 1
        scanlist = sorted(scansummary.keys())
        starttime_mjd = scansummary[scanlist[scan]]['0']['BeginTime']
    else:
        print 'No pickle of initializing info found. Making anew...'
        pkl = open(pklname, 'wb')
        ms.open(filename)
        spwinfo = ms.getspectralwindowinfo()
        scansummary = ms.getscansummary()
        ms.selectinit(datadescid=0)  # reset select params for later data selection
        ms.selectpolarization(selectpol)

        scanlist = sorted(scansummary.keys())
        starttime_mjd = scansummary[scanlist[scan]]['0']['BeginTime']
        d['inttime'] = scansummary[scanlist[scan]]['0']['IntegrationTime']
        print 'Initializing integration time (s):', d['inttime']

        ms.iterinit(['TIME'], iterint*d['inttime'])
        ms.iterorigin()
        da = ms.getdata([datacol, 'axis_info'], ifraxis=True)
        ms.close()

        d['nbl'] = da[datacol].shape[2]
        bls = da['axis_info']['ifr_axis']['ifr_shortname']
        d['blarr'] = n.array([[int(bls[i].split('-')[0]),int(bls[i].split('-')[1])] for i in xrange(len(bls))])
        d['npol'] = len(selectpol)
        d['npol_orig'] = da[datacol].shape[0]
        print 'Initializing %d polarizations' % (d['npol'])

        pickle.dump((d['npol_orig'], d['nbl'], d['blarr'], d['inttime'], spwinfo, scansummary), pkl)
        pkl.close()

    # set ants
    d['ants'] = n.unique(d['blarr'])
    d['nants'] = len(n.unique(d['blarr']))
    print 'Initializing nants:', d['nants']
    print 'Initializing nbl:', d['nbl']

    # define list of spw keys (may not be in order!)
    freqs = []
    for i in spwinfo.keys():
        freqs.append(spwinfo[i]['Chan1Freq'])
    d['spwlist'] = n.array(sorted(zip(freqs, spwinfo.keys())))[:,1][spw].astype(int)  # spwlist defines order of spw to iterate in freq order

    d['freq_orig'] = n.array([])
    for spw in d['spwlist']:
        nch = spwinfo[str(spw)]['NumChan']
        ch0 = spwinfo[str(spw)]['Chan1Freq']
        chw = spwinfo[str(spw)]['ChanWidth']
        d['freq_orig'] = n.concatenate( (d['freq_orig'], (ch0 + chw * n.arange(nch)) * 1e-9) )

    d['freq'] = d['freq_orig'][d['chans']]
    d['nchan'] = len(d['chans'])
    print 'Initializing nchan:', d['nchan']

    # set requested time range based on given parameters
    timeskip = d['inttime']*nskip
    starttime = qa.getvalue(qa.convert(qa.time(qa.quantity(starttime_mjd+timeskip/(24.*60*60),'d'),form=['ymd'], prec=9)[0], 's'))[0]
    stoptime = qa.getvalue(qa.convert(qa.time(qa.quantity(starttime_mjd+(timeskip+(nints+1)*d['inttime'])/(24.*60*60), 'd'), form=['ymd'], prec=9)[0], 's'))[0]  # nints+1 to be avoid buffer running out and stalling iteration
    print 'First integration of scan:', qa.time(qa.quantity(starttime_mjd,'d'),form=['ymd'],prec=9)[0]
    print
    print 'Reading scan', str(scanlist[scan]) ,'for times', qa.time(qa.quantity(starttime_mjd+timeskip/(24.*60*60),'d'),form=['hms'], prec=9)[0], 'to', qa.time(qa.quantity(starttime_mjd+(timeskip+(nints+1)*d['inttime'])/(24.*60*60), 'd'), form=['hms'], prec=9)[0]

    # read data into data structure
    ms.open(filename)
    if len(d['spwlist']) == 1:
        ms.selectinit(datadescid=d['spwlist'][0])
    else:
        ms.selectinit(datadescid=0, reset=True)    # reset includes spw in iteration over time
    selection = {'time': [starttime, stoptime]}
    ms.select(items = selection)
    ms.selectpolarization(selectpol)
    ms.iterinit(['TIME'], iterint*d['inttime'], 0, adddefaultsortcolumns=False)   # read with a bit of padding to get at least nints
    iterstatus = ms.iterorigin()
    d['itercount1'] = 0
    d['l0'] = 0.; d['m0'] = 0.

    return iterstatus

def readiter(d):
    """ Iterates over ms.
    Returns everything needed for analysis as tuple.
    """

    da = ms.getdata([d['datacol'],'axis_info','u','v','w','flag','data_desc_id'], ifraxis=True)
#    spws = n.unique(da['data_desc_id'])    # spw in use
#    good = n.where((da['data_desc_id']) == spws[0])[0]   # take first spw
    good = n.where((da['data_desc_id']) == d['spwlist'][0])[0]   # take first spw
    time0 = da['axis_info']['time_axis']['MJDseconds'][good]
    data0 = n.transpose(da[d['datacol']], axes=[3,2,1,0])[good]
    if d['telcalfile']:    # apply telcal solutions
        if len(d['spwlist']) > 1:
            spwbin = d['spwlist'][0]
        else:
            spwbin = 0
        chanfreq = da['axis_info']['freq_axis']['chan_freq'][:,spwbin]
        sols = applytelcal.solutions(d['telcalfile'], chanfreq)
        for i in range(len(d['selectpol'])):
            try:
                sols.setselection(d['telcalcalibrator'], time0[0]/(24*3600), d['selectpol'][i], verbose=1)   # chooses solutions closest in time that match pol and source name
                sols.apply(data0, d['blarr'], i)
            except:
                pass

    flag0 = n.transpose(da['flag'], axes=[3,2,1,0])[good]
    u0 = da['u'].transpose()[good] * d['freq_orig'][0] * (1e9/3e8)
    v0 = da['v'].transpose()[good] * d['freq_orig'][0] * (1e9/3e8)
    w0 = da['w'].transpose()[good] * d['freq_orig'][0] * (1e9/3e8)
    if len(d['spwlist']) > 1:
        for spw in d['spwlist'][1:]:
            good = n.where((da['data_desc_id']) == spw)[0]
            data1 = n.transpose(da[d['datacol']], axes=[3,2,1,0])[good]
            if d['telcalfile']:    # apply telcal solutions
                chanfreq = da['axis_info']['freq_axis']['chan_freq'][:,spw]
                sols = applytelcal.solutions(d['telcalfile'], chanfreq)
                for i in range(len(d['selectpol'])):
                    try:
                        sols.setselection(d['telcalcalibrator'], time0[0]/(24*3600), d['selectpol'][i], verbose=1)   # chooses solutions closest in time that match pol and source name
                        sols.apply(data1, d['blarr'], i)
                    except:
                        pass

            data0 = n.concatenate( (data0, data1), axis=2 )
            flag0 = n.concatenate( (flag0, n.transpose(da['flag'], axes=[3,2,1,0])[good]), axis=2 )

    d['iterstatus1'] = ms.iternext()
    return data0[:,:,d['chans'],:], flag0[:,:,d['chans'],:], u0, v0, w0, time0

def readloop(d, eproc, emove):
    """ Data generating stage of parallel data function.
    data is either read into 'data' buffer, when ready
    this keeps data reading bottleneck to 1x the read time.
    """

    iterint = d['iterint']; nbl = d['nbl']; nchan = d['nchan']; npol = d['npol']
    datacap, flagcap, ucap, vcap, wcap, timecap = readiter(d)    # read "cap", a hack to make sure any single iteration has enough integrations (artifact of irregular inttime)
    print 'Reading first iteration with shape', datacap.shape

    # dynamically set uvgrid
    if ((d['size'] == 0) or (d['res'] == 0)):
        d['size'] = n.round(max(ucap.max() - ucap.min(), vcap.max() - vcap.min())).astype('int')
        d['res'] = n.round(25./(3e-1/d['freq'][len(d['freq'])/2])).astype('int')    # full field of view. assumes freq in GHz
        print 'Set size to %d and res to %d' % (d['size'], d['res'])

    while 1:
#            name = mp.current_process().name
#            print '%s: filling buffer' % name
        datanext, flagnext, unext, vnext, wnext, timenext = readiter(d)
        print 'Reading next %d ints from iter %d' % (len(datanext), d['itercount1']+iterint)
        datanext = n.vstack((datacap,datanext))
        flagnext = n.vstack((flagcap,flagnext))
        unext = n.vstack((ucap,unext))
        vnext = n.vstack((vcap,vnext))
        wnext = n.vstack((wcap,wnext))
        timenext = n.concatenate((timecap,timenext))
        # select just the next iteration's worth of data and metadata. leave rest for next iteration's buffer cap.
        if len(datanext) >= iterint:
            data1 = datanext[:iterint]
            datacap = datanext[iterint:]   # save rest for next iteration
            flag1 = flagnext[:iterint]
            flagcap = flagnext[:iterint]
            u1 = unext[:iterint]
            ucap = unext[:iterint]
            v1 = vnext[:iterint]
            vcap = vnext[:iterint]
            w1 = wnext[:iterint]
            wcap = wnext[:iterint]
            time1 = timenext[:iterint]
            timecap = timenext[:iterint]

            # wait for signal to move everything to processing buffers
            emove.wait()
            emove.clear()
            # move metadata from read buffer "1" to processing buffer "0"
            data[:] = data1[:]
            flag[:] = flag1[:]
            u[:] = u1[:]
            v[:] = v1[:]
            w[:] = w1[:]
            time[:] = time1[:]
            d['itercount'] = d['itercount1']
            d['iterstatus'] = d['iterstatus1']
            d['itercount1'] += iterint
            eproc.set()    # reading buffer filled, start processing
        else:
            print 'End of data (in buffer)'
            ms.iterend()
            ms.close()
            break
        if not d['iterstatus']:
            print 'End of data (iterator)'
            ms.iterend()
            ms.close()
            break

def processloop(d, triples, eproc, emove):
    """ Processing stage of parallel data function. 
    Only processes from data. Assumes a "first in, first out" model, where 'data' defines next buffer to process.
    Event triggered by readloop when 'data' is filled.
    """

    while 1:
        eproc.wait()
        eproc.clear()
#        name = mp.current_process().name
#        print '%s: processing data' % name

        # define imaging params, if needed
#        data = n.ma.masked_array(data, flag == 1)       # define working copy of data as masked array. note that mask does not change, so dedispersion does not move mask in time.
        d['datadelay'] = lib.calc_delay(d, 0)   
#        transient_search(n.ma.masked_array(data, flag == 1), d, triples)
        transient_search(n.ma.masked_array(data, data == 0j), d, triples)    # flag based on data set to 0
        emove.set()
        if not d['iterstatus']:
            print 'End of processloop'
            break

def readloop2(d, eproc, emove):
    """ Profiles readloop 
    """
    cProfile.runctx('readloop(d, eproc, emove)', globals(), locals(), 'readloop.prof')

def processloop2(d, triples, eproc, emove):
    """ Profiles processloop
    """
    cProfile.runctx('processloop(d, triples, eproc, emove)', globals(), locals(), 'processloop.prof')

def pipe_thread(filename, profile='default', nints=256, nskip=0, iterint=128, spw=[0], chans=range(5,59), dmarr=[0.], fwhmsurvey=0.5, fwhmfield=0.5, selectpol=['I'], scan=0, datacol='corrected_data', size=30000, res=100, sigma_bisp=5., sigma_image=5., filtershape=None, calibratedfilter=True, specmodfilter=1.5, searchtype='imageallstat', telcalfile='', telcalcalibrator='', savecands=0):
    """ Threading for parallel data reading and processing.
    Either side can be faster than the other, since data are held for processing in shared buffer.
    size/res define uvgrid parameters. if either set to 0, then they are dynamically set to image full field of view and include all visibilities.
    Search types include 'bispectrum', 'imageallstat', 'imageallpng', 'imageone'.
    """

    # set up thread management and shared memory and metadata
    global data, flag, u, v, w, time

    mgr = mp.Manager()
    d = mgr.dict()
    eproc = mp.Event()      # event signalling to begin processing
    emove = mp.Event()      # event signalling to move data into processing buffers (data, flag, u, v, w, time)

    # define basic shared params
    d['filename'] = filename
    d['spw'] = spw
    d['datacol'] = datacol
    d['dmarr'] = dmarr
    d['scan'] = scan
    d['nskip'] = nskip
    d['nints'] = nints       # total ints to iterate over
    d['iterint'] = iterint    # time step for msiter
    d['timescales'] = [1]
    d['chans'] = chans
    d['nchan'] = len(chans)
    d['selectpol'] = selectpol
    d['npol'] = len(selectpol)
    d['filtershape'] = filtershape
    d['bgwindow'] = 4
    d['sigma_bisp'] = sigma_bisp
    d['sigma_image'] = sigma_image
    d['size'] = size
    d['res'] = res
    d['calibratedfilter'] = calibratedfilter
    d['specmodfilter'] = specmodfilter     # fudge factor for spectral modulation. 1==ideal, 0==do not apply, >1==non-ideal broad-band signal
    d['searchtype'] = searchtype
    d['delaycenters'] = calc_hexcenters(fwhmsurvey, fwhmfield)
    d['telcalfile'] = telcalfile
    d['telcalcalibrator'] = telcalcalibrator
    d['savecands'] = savecands

    # define basic data state
    print 'Preparing to read...'
    d['iterstatus'] = readprep(d)

    # must come after readprep
    d['datadelay'] = lib.calc_delay(d, 0)   
    triples = lib.make_triples(d)
    d['ntr'] = len(triples)

    # initialize fields
    masterloc = {}           # "location" of candidate: dt, integration (over all data in pipe), dmbin
    masterprop = {}          # properties of candidate: snr and std of bispectra for dedispersed trial
    masterdata = {}          # data around the candidate: bisp lightcurve and spectrogram
    masternoise = {}         # noise per bl as measured by detect_bispectra
    for (ra, dec) in d['delaycenters']:
        masterloc[(ra,dec)] = []
        masterprop[(ra,dec)] = []
        masterdata[(ra,dec)] = []
        masternoise[(ra,dec)] = []

    # define candidate file
    if d['savecands']:
        tt = timestamp.localtime()
        print 'Start time: %s_%s_%s:%s:%s:%s' % (tt.tm_year, tt.tm_mon, tt.tm_mday, tt.tm_hour, tt.tm_min, tt.tm_sec)
        timestring = '%s_%s_%s:%s:%s:%s' % (tt.tm_year, tt.tm_mon, tt.tm_mday, tt.tm_hour, tt.tm_min, tt.tm_sec)
        d['masterfile'] = 'mastercands_'+timestring+'.pkl'
        picklabledict = d.copy()
        pkl = open(d['masterfile'], 'wb')
        pickle.dump(picklabledict, pkl)
        pickle.dump((masterloc, masterprop, masterdata, masternoise), pkl)
        pkl.close()

    # create shared data arrays
    data_mem = mp.Array(ctypes.c_double, iterint*d['nbl']*d['nchan']*d['npol']*2)    # x2 to store complex values in double array
    flag_mem = mp.Array(ctypes.c_bool, iterint*d['nbl']*d['nchan']*d['npol'])
    u_mem = mp.Array(ctypes.c_double, iterint*d['nbl'])
    v_mem = mp.Array(ctypes.c_double, iterint*d['nbl'])
    w_mem = mp.Array(ctypes.c_double, iterint*d['nbl'])
    time_mem = mp.Array(ctypes.c_double, iterint)
# new way is to convert later
    data = numpyview(data_mem, 'complex', (iterint, d['nbl'], d['nchan'], d['npol']))
    flag = numpyview(flag_mem, 'bool', (iterint, d['nbl'], d['nchan'], d['npol']))
    u = numpyview(u_mem, 'double', (iterint, d['nbl']))
    v = numpyview(v_mem, 'double', (iterint, d['nbl']))
    w = numpyview(w_mem, 'double', (iterint, d['nbl']))
    time = numpyview(time_mem, 'double', (iterint))

    try:
        # start processes
        pread = mp.Process(target=readloop, args=(d,eproc,emove))
        pread.start()
        pproc = mp.Process(target=processloop, args=(d,triples,eproc,emove))
        pproc.start()

        # trigger events to allow moving data to working area
        emove.set()

        # wait for threads to end (when read iteration runs out of data)
        pread.join()
        pproc.join()

    except KeyboardInterrupt:
        print 'Ctrl-C received. Shutting down threads...'
        pread.terminate()
        pproc.terminate()
        pread.join()
        pproc.join()

    return d.copy()
