import numpy as np
cimport numpy as np
cimport cython
from cython.view cimport array as cvarray

DTYPE = np.complex64
ctypedef np.complex64_t DTYPE_t

def sigma_clip(np.ndarray[np.float32_t, ndim=1] arr, float sigma=3):
    """ Function takes 1d array of values and returns the sigma-clipped min and max scaled by value "sigma".
    """

    assert arr.dtype == np.float32

    cdef np.ndarray[np.int_t, ndim=1] cliparr = np.arange(len(arr))
    cdef float mean
    cdef float std

    arr = np.append(arr,[np.float32(1)])    # append superfluous item to trigger loop
    while len(cliparr) != len(arr):
        arr = arr[cliparr]
        mean = arr.mean()
        std = arr.std()
        cliparr = np.where((arr < mean + sigma*std) & (arr > mean - sigma*std) & (arr != 0) )[0]
#        print 'Clipping %d from array of length %d' % (len(arr) - len(cliparr), len(arr))
    return mean - sigma*std, mean + sigma*std

cpdef make_triples(d):
    """ Calculates and returns data indexes (i,j,k) for all closed triples.
    """

    ants = d['ants']
    nants = d['nants']
    blarr = d['blarr']

    cdef int t
    cdef int ant1
    cdef int ant2
    cdef int ant3
    cdef np.ndarray[np.int_t, ndim=2] triples = np.zeros((nants*(nants-1)*(nants-2)/6, 3), dtype='int')

    # first make triples indexes in antenna numbering
    anttrips = []
    for i in ants:
        for j in ants[list(ants).index(i)+1:]:
            for k in ants[list(ants).index(j)+1:]:
                anttrips.append([i,j,k])
                
    # next return data indexes for triples
    for t in xrange(len(anttrips)):
        ant1 = anttrips[t][0]
        ant2 = anttrips[t][1]
        ant3 = anttrips[t][2]
        try:
            bl1 = np.where( (blarr[:,0] == ant1) & (blarr[:,1] == ant2) )[0][0]
            bl2 = np.where( (blarr[:,0] == ant2) & (blarr[:,1] == ant3) )[0][0]
            bl3 = np.where( (blarr[:,0] == ant1) & (blarr[:,1] == ant3) )[0][0]
            triples[t,0] = bl1
            triples[t,1] = bl2
            triples[t,2] = bl3
        except IndexError:
            continue

    return triples

@cython.profile(False)
cdef np.ndarray[np.float32_t, ndim=2] fringe_rotation(float dl, float dm, np.ndarray[np.float32_t, ndim=1] u, np.ndarray[np.float32_t, ndim=1] v, np.ndarray[np.float32_t, ndim=1] freq): 
    return np.exp(-2j*3.1415*(dl*np.outer(u,freq) + dm*np.outer(v,freq)))

@cython.wraparound(False)
@cython.boundscheck(False)
cpdef phaseshift(np.ndarray[DTYPE_t, ndim=4] data0, d, float l1, float m1, np.ndarray[np.float32_t, ndim=2] u, np.ndarray[np.float32_t, ndim=2] v, verbose=0):
    """ Shift phase center to (l1, m1).
    Assumes single uv over all times in data0. Reasonable for up to a second or so of data.
    """

    cdef np.ndarray[complex, ndim=2] frot
    cdef np.ndarray[float, ndim=1] freq = d['freq']
    cdef np.ndarray[float, ndim=1] freq_orig = d['freq_orig']
    cdef float dl = l1 - d['l0']
    cdef float dm = m1 - d['m0']
    cdef unsigned int i
    cdef unsigned int j
    cdef unsigned int k
    cdef unsigned int l

    shape = np.shape(data0)
    cdef unsigned int len0 = shape[0]
    cdef unsigned int len1 = shape[1]
    cdef unsigned int len2 = shape[2]
    cdef unsigned int len3 = shape[3]

    if (dl != 0.) or (dm != 0.):
        frot = fringe_rotation(dl, dm, u[len0/2], v[len0/2], freq/freq_orig[0])
        for j in xrange(len1):
            for i in xrange(len0):
                for l in xrange(len3):    # iterate over pols
                    for k in xrange(len2):
                        data0[i,j,k,l] = data0[i,j,k,l] * frot[j,k]
    else:
        if verbose:
            print 'No phase rotation needed'

    d['l0'] = l1
    d['m0'] = m1

@cython.wraparound(False)
@cython.boundscheck(False)
cpdef phaseshift_threaded(np.ndarray[DTYPE_t, ndim=4] data0, d, float l1, float m1, np.ndarray[np.float32_t, ndim=1] u, np.ndarray[np.float32_t, ndim=1] v, verbose=0):
    """ Shift phase center to (l1, m1).
    Assumes single uv over all times in data0. Reasonable for up to a second or so of data.
    """

    cdef np.ndarray[complex, ndim=2] frot
    cdef np.ndarray[float, ndim=1] freq = d['freq']
    cdef np.ndarray[float, ndim=1] freq_orig = d['freq_orig']
    cdef float dl = l1 - d['l0']
    cdef float dm = m1 - d['m0']
    cdef unsigned int i
    cdef unsigned int j
    cdef unsigned int k
    cdef unsigned int l

    shape = np.shape(data0)
    cdef unsigned int len0 = shape[0]
    cdef unsigned int len1 = shape[1]
    cdef unsigned int len2 = shape[2]
    cdef unsigned int len3 = shape[3]

    if (dl != 0.) or (dm != 0.):
        frot = fringe_rotation(dl, dm, u, v, freq/freq_orig[0])
        for j in xrange(len1):
            for i in xrange(len0):
                for l in xrange(len3):    # iterate over pols
                    for k in xrange(len2):
                        data0[i,j,k,l] = data0[i,j,k,l] * frot[j,k]
    else:
        if verbose:
            print 'No phase rotation needed'

cpdef np.ndarray[np.int_t, ndim=1] calc_delay(np.ndarray[float, ndim=1] freq, float inttime, float dm):
    """ Function to calculate delay for each channel, given dm
    freq in GHz, inttime in s.
    """

    cdef float freqref = freq[len(freq)-1]

    return np.round((4.2e-3 * dm * (1/(freq*freq) - 1/(freqref*freqref)))/inttime,0).astype(np.int16)

cpdef calc_resample(float chanwidth, float midfreq, float dm, float inttime):
    """ Function to calculate resmapling factor.
    freqs in GHz. inttime in s. returns intrachannel smearing by number of integrations.
    """

    return np.round(np.sqrt( (8.3e-3 * dm * chanwidth / midfreq**3)**2 + inttime**2)/inttime, 0).astype(np.int16)

@cython.wraparound(False)
@cython.boundscheck(False)
cpdef dedisperse(np.ndarray[DTYPE_t, ndim=4, mode='c'] data0, d, float dm, int verbose=0):
    """ dedisperse the data in place
    replaces data0 in place. dm algorithm on only accurate if moving "up" in dm space.
    """

    cdef unsigned int i
    cdef unsigned int iprime
    cdef unsigned int j
    cdef unsigned int k
    cdef unsigned int l
    cdef unsigned int indmin
    cdef unsigned int indmax
    cdef int shift
    cdef int len0
    cdef unsigned int len1
    cdef unsigned int len2
    cdef unsigned int len3
# making copy is slower, but returns actual roll of original data0. without, we can only move "up" in dm space.
#    cdef np.ndarray[DTYPE_t, ndim=4] data1 = np.empty_like(data0)

    # calc relative delay per channel. only shift minimally
    cdef np.ndarray[short, ndim=1] newdelay = calc_delay(d['freq'], d['inttime'], dm)
    cdef np.ndarray[short, ndim=1] relativedelay = newdelay - d['datadelay']

    shape = np.shape(data0)
    len0 = shape[0]
    len1 = shape[1]
    len2 = shape[2]
    len3 = shape[3]

#    print relativedelay
#    for shift in np.unique(relativedelay):
    for j in xrange(len1):
        for k in xrange(len2):
            shift = relativedelay[k]
            if shift > 0:
                for l in xrange(len3):
# following is cython shortcut for 'iprime = np.mod(i+shift, len0)'
                    for i in xrange(len0):
                        iprime = i+shift
#                    print i, iprime
                        if iprime >= 0 and iprime < len0:    # ignore edge cases
                            data0[i,j,k,l] = data0[iprime,j,k,l]
                        elif iprime >= len0:    # set nonsense shifted data to zero
                            data0[i,j,k,l] = 0j
            elif shift < 0:
                print 'negative delay found. this dedispersion algorithm only works for positive delays. ignoring...'

# alternatives
#                            data1[i,j,k,l] = data0[iprime,j,k,l]
#            data0[:,:,indmin:indmax+1,:] = np.roll(data0.take(range(indmin,indmax+1), axis=2), -1*shift, axis=0)

    if verbose != 0:
        print 'Dedispersed for DM=%d' % dm

    # new delay values
    d['datadelay'] = newdelay
#    return data1        # if doing proper roll, need to return copy

@cython.wraparound(False)
@cython.boundscheck(False)
cpdef dedisperse_threaded(np.ndarray[DTYPE_t, ndim=4, mode='c'] data0, d, float dm, int verbose=0):
    """ dedisperse the data in place
    replaces data0 in place. dm algorithm on only accurate if moving "up" in dm space.
    assumes unshifted data.
    """

    cdef unsigned int i
    cdef unsigned int iprime
    cdef unsigned int j
    cdef unsigned int k
    cdef unsigned int l
    cdef unsigned int indmin
    cdef unsigned int indmax
    cdef int shift
    shape = np.shape(data0)
    cdef unsigned int len0 = shape[0]
    cdef unsigned int len1 = shape[1]
    cdef unsigned int len2 = shape[2]
    cdef unsigned int len3 = shape[3]
# making copy is slower, but returns actual roll of original data0. without, we can only move "up" in dm space.
#    cdef np.ndarray[DTYPE_t, ndim=4] data1 = np.empty_like(data0)

    # calc relative delay per channel. only shift minimally
    cdef np.ndarray[short, ndim=1] relativedelay = calc_delay(d['freq'], d['inttime'], dm)

#    print relativedelay
#    for shift in np.unique(relativedelay):
    for j in xrange(len1):
        for k in xrange(len2):
            shift = relativedelay[k]
            if shift > 0:
                for l in xrange(len3):
# following is cython shortcut for 'iprime = np.mod(i+shift, len0)'
                    for i in xrange(len0):
                        iprime = i+shift
#                    print i, iprime
                        if iprime >= 0 and iprime < len0:    # ignore edge cases
                            data0[i,j,k,l] = data0[iprime,j,k,l]
                        elif iprime >= len0:    # set nonsense shifted data to zero
                            data0[i,j,k,l] = 0j
            elif shift < 0:
                print 'negative delay found. this dedispersion algorithm only works for positive delays. ignoring...'

# alternatives
#                            data1[i,j,k,l] = data0[iprime,j,k,l]
#            data0[:,:,indmin:indmax+1,:] = np.roll(data0.take(range(indmin,indmax+1), axis=2), -1*shift, axis=0)

    if verbose != 0:
        print 'Dedispersed for DM=%d' % dm

@cython.wraparound(False)
@cython.boundscheck(False)
cpdef dedisperse_resample(np.ndarray[DTYPE_t, ndim=4, mode='c'] data0, np.ndarray[float, ndim=1] freq, float inttime, float dm, unsigned int resample, int verbose=0):
    """ dedisperse the data and resample in place. only fraction of array is useful data.
    dm algorithm on only accurate if moving "up" in dm space.
    assumes unshifted data.
    only does resampling by dt. no dm resampling for now.
    """

    cdef unsigned int i
    cdef unsigned int iprime
    cdef unsigned int j
    cdef unsigned int k
    cdef unsigned int l
    cdef unsigned int r
    cdef unsigned int indmin
    cdef unsigned int indmax
    cdef int shift
    shape = np.shape(data0)
    cdef unsigned int len0 = shape[0]
    cdef unsigned int len1 = shape[1]
    cdef unsigned int len2 = shape[2]
    cdef unsigned int len3 = shape[3]
    # calc relative delay per channel. only shift minimally
    cdef np.ndarray[short, ndim=1] relativedelay = calc_delay(freq, inttime, dm)
    cdef unsigned int newlen0 = len0/resample

    lostsamples = np.mod(len0, resample)
    if lostsamples:
        print 'Warning: resampling drops last %d int(s).' % (lostsamples)

    for j in xrange(len1):
        for l in xrange(len3):
            for k in xrange(len2):
                shift = relativedelay[k]
                for i in xrange(newlen0):
                    iprime = i*resample+shift
                    if iprime >= 0 and iprime < len0-(resample-1):    # if within bounds of unshifted data with resample stepping
                        data0[i,j,k,l] = data0[iprime,j,k,l]
                        if resample > 1:
                            for r in xrange(1,resample):
                                data0[i,j,k,l] = data0[i,j,k,l] + data0[iprime+r,j,k,l]
                            data0[i,j,k,l] = data0[i,j,k,l]/resample
                    elif iprime >= len0-(resample):    # set nonsense shifted data to zero
                        data0[i,j,k,l] = 0j

    if verbose != 0:
        print 'Dedispersed for DM=%d' % dm

@cython.wraparound(False)
@cython.boundscheck(False)
cpdef dedisperse_resample2(np.ndarray[DTYPE_t, ndim=4, mode='c'] data0, np.ndarray[float, ndim=1] freq, float inttime, float dm, unsigned int resample, int verbose=0):
    """ dedisperse the data and resample with fixed value in place. only fraction of array is useful data.
    dm algorithm on only accurate if moving "up" in dm space.
    assumes unshifted data.
    """

    cdef unsigned int i
    cdef unsigned int iprime
    cdef unsigned int j
    cdef unsigned int k
    cdef unsigned int l
    cdef unsigned int r
    cdef unsigned int indmin
    cdef unsigned int indmax
    cdef int shift
    shape = np.shape(data0)
    cdef unsigned int len0 = shape[0]
    cdef unsigned int len1 = shape[1]
    cdef unsigned int len2 = shape[2]
    cdef unsigned int len3 = shape[3]
    cdef float chanwidth = freq[1] - freq[0]
    cdef float midfreq = freq[len(freq)/2]
    cdef unsigned int newlen0 = len0/resample

    if np.mod(len0, resample):
        print 'Warning: data length is not an integer number of resamples. Last int(s) will be lost.'

    # calc relative delay per channel. only shift minimally
    cdef np.ndarray[short, ndim=1] relativedelay = calc_delay(freq, inttime, dm)

    for j in xrange(len1):
        for l in xrange(len3):
            for k in xrange(len2):
                shift = relativedelay[k]
                for i in xrange(newlen0):
                    iprime = i*resample+shift
                    if iprime >= 0 and iprime < len0-(resample-1):    # if within bounds of unshifted data with resample stepping
                        data0[i,j,k,l] = data0[iprime,j,k,l]
                        if resample > 1:
                            for r in xrange(1,resample):
                                data0[i,j,k,l] = data0[i,j,k,l] + data0[iprime+r,j,k,l]
                            data0[i,j,k,l] = data0[i,j,k,l]/resample
                    elif iprime >= len0-resample:    # set nonsense shifted data to zero
                        data0[i,j,k,l] = 0j

    if verbose != 0:
        print 'Dedispersed for DM=%d' % dm

@cython.wraparound(False)
@cython.boundscheck(False)
cpdef make_bispectra(np.ndarray[DTYPE_t, ndim=4] data0, np.ndarray[np.int_t, ndim=2] triples):
    """ Makes bispectra by integration, which assumes data have been dedispersed.
    """
    shape = np.shape(data0)
    cdef unsigned int len0 = shape[0]
    cdef unsigned int len1 = shape[1]
    cdef unsigned int len2 = shape[2]
    cdef unsigned int len3 = shape[3]
    cdef unsigned int i
    cdef unsigned int j
    cdef unsigned int k
    cdef unsigned int p
    cdef unsigned int t
    cdef unsigned int l
    cdef unsigned int m
    cdef unsigned int n
    cdef complex sum
    cdef np.ndarray[np.complex64_t, ndim=2] bispectra = np.zeros(shape=(len0, len(triples)), dtype='complex64')
    cdef np.ndarray[np.complex64_t, ndim=2] data0m = np.zeros(shape=(len0,len1), dtype='complex64')

    # mean over chans and pols
    for i in xrange(len0):
        for j in xrange(len1):
            sum = 0.
            for p in xrange(len3):
                for k in xrange(len2):                   
                    sum += data0[i,j,k,p]
            data0m[i,j] = sum/(len2*len3)
#    data0m = data0.mean(axis=2).mean(axis=2)

    for t in xrange(len(triples)):
        l = triples[t, 0]
        m = triples[t, 1]
        n = triples[t, 2]
        for i in xrange(len0):
            bispectra[i, t] = data0m[i, l] * data0m[i, m] * data0m[i, n].conjugate()

    return bispectra

cpdef dataflag(datacal, sigma, d, mode='blstd'):
    """ Flagging function
    """

    sh = datacal.shape
    iterint = sh[0]
    nbl = sh[1]
    nchan = sh[2]
    npol = sh[3]
    datacal = np.ma.masked_array(datacal, datacal == 0j)

    chperspw = len(d['freq_orig'])/len(d['spw'])
    for spw in range(len(d['spw'])):
        freqs = d['freq_orig'][spw*chperspw:(spw+1)*chperspw]  # find chans for spw. only works for 2 or more sb
        chans = np.array([i for i in range(len(d['freq'])) if d['freq'][i] in freqs])
        for pol in range(npol):
            if mode == 'blstd':
                blstd = datacal[:,:,chans,pol].std(axis=1)     # high std over baseline as fcn of time,freq selects RFI well. also bright transients off axis...
                blstd = np.ma.masked_where(blstd > np.median(blstd) + sigma*blstd.std(), blstd)   # first mask out high points
                blstd = np.ma.masked_where(blstd > np.median(blstd) + sigma*blstd.std(), blstd)   # repeat with new std
                ww = np.where(blstd.data > np.median(blstd) + sigma*blstd.std())     # then measure points to flag based on a third std threshold
                flagfrac = float(len(ww[0])*nbl)/datacal.size
                print 'Blstd flagging for (spw %d, pol %d), %d sigma flag fraction %3.5f %%' % (spw,pol,sigma,flagfrac*100)
                for i in xrange(len(ww[0])):
                    datacal[ww[0][i],:,chans[ww[1][i]],pol] = 0j     # must be a better way of doing this!?
            elif mode == 'badchan':
                blstdmean = datacal[:,:,chans,pol].std(axis=1).mean(axis=0)
                blstdmean = np.ma.masked_array(blstdmean, blstdmean==0j)
                ww = np.where(blstdmean > blstdmean.mean() + sigma * blstdmean.std())
                flagfrac = float(len(ww[0])*len(chans)*len(datacal))/datacal.shape[2]
                print 'Badchan flagging for (spw %d, pol %d), %d sigma flag fraction %3.5f %%' % (spw,pol,sigma,flagfrac*100)
                datacal[:,:,chans[list(ww[0])],pol] = 0j     # must be a better way of doing this!?
            elif mode == 'medcht':
                meanamp = np.abs(datacal[:,:,chans,pol]).mean(axis=1)
                meanamp = np.ma.masked_array(meanamp, meanamp==0j)
                medall = np.median(meanamp)
                medch = np.median(meanamp, axis=0)
                medt = np.median(meanamp, axis=1)
                badch = np.where(medch > sigma*medall)
                if len(badch[0]):
                    if len(badch[0]) > 0.2*len(chans):
                        print 'More than 20\% bad channels in spw. Flagging it all...'
                        badch = [range(len(chans))]
                    datacal[:,:,chans[badch],pol] = 0j
                badt = np.where(medt > sigma*medall)
                if len(badt[0]):
                    datacal[:,:,chans,pol][badt] = 0j
                flagfrac = nbl*float(len(chans)*len(badt[0]) + len(datacal)*len(badch[0]))/datacal.size
                print 'Median flagging of %d chans and %d ints for (spw %d, pol %d). %.1f sigma flag fraction %3.5f %%' % (len(badch[0]), len(badt[0]), spw,pol,sigma,flagfrac*100)
            elif mode == 'badbp':
                if pol == 1:   # dual-pol algorithm, so skip second pol
                    continue
                bpa = np.abs(datacal[:,:,chans,:]).mean(axis=0).mean(axis=1)
                bpa_ant = np.array([ (bpa[np.where(np.any(d['blarr'] == i, axis=1))[0]]).mean(axis=0) for i in np.unique(d['blarr']) ])
                bpa_ant = np.ma.masked_invalid(bpa_ant)
                ww = np.where(bpa_ant > np.median(bpa_ant) + sigma * bpa_ant.std())
                badants = np.unique(d['blarr'])[ww[0]]
                badpols = ww[1]
                print 'Bad basepol flagging for spw %d on ants/pols %s/%s' % (spw,badants,badpols)
                for badpol in np.unique(badpols):
                    badantpol = badants[np.where(badpol == badpols)]
                    if float(len(badantpol))/d['nants'] > 0.2:
                        badants = np.concatenate( (badants,np.unique(d['blarr'])) )
                        badpols = np.concatenate( (badpols,[badpol]*d['nants']) )
                        print 'More than 0.2 of ants in pol %d hit bad. Flagging all.' % (badpol)
                for badant, badpol in zip(badants, badpols):
                    badantind = np.where([badant in aa for aa in d['blarr']])[0]
                    for chan in chans:
                        datacal[:,badantind,chan,badpol] = 0j
            elif mode == 'ring':
                # find bls with ringing due to strong RFI. sigma=2 seems to work well.
                spfft = np.abs(np.fft.ifft(datacal[...,chans,pol].mean(axis=0), axis=1))   # delay spectrum of mean data in time
                spfft = np.ma.masked_array(spfft, spfft = 0)
                badbls = np.where(spfft[:,3*len(chans)/8:len(chans)/2].mean(axis=1) > sigma*np.median(spfft[:,1:], axis=1))[0]  # find bls with spectral power at large delay. ignore dc in case this is cal scan.
                print 'Ringing flagging for (spw %d, pol %d) on bls %s' % (spw,pol,badbls)
                if len(badbls) > nbl/3:    # if many bls affected, flag all
                    datacal[:,:,chans,pol] = 0j     
                else:
                    for badbl in badbls:
                        datacal[:,badbl,chans,pol] = 0j
            else:
                print 'Not flagging based on data quality.'
