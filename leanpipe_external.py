import numpy as np

def sigma_clip(arr,sigma=3):
    """ Function takes 1d array of values and returns the sigma-clipped min and max scaled by value "sigma".
    """

    cliparr = range(len(arr))  # initialize
    arr = np.append(arr,[1])    # append superfluous item to trigger loop
    while len(cliparr) != len(arr):
        arr = arr[cliparr]
        mean = arr.mean()
        std = arr.std()
        cliparr = np.where((arr < mean + sigma*std) & (arr > mean - sigma*std) & (arr != 0) )[0]
#        print 'Clipping %d from array of length %d' % (len(arr) - len(cliparr), len(arr))
    return mean - sigma*std, mean + sigma*std

def make_triples(d, amin=0, amax=0):
    """ Calculates and returns data indexes (i,j,k) for all closed triples.
    amin and amax define range of antennas (with index, in order). only used if nonzero.
    """

    ants = d['ants']
    nants = d['nants']
    if amax == 0:
        amax = nants
    blarr = d['blarr']

    # first make triples indexes in antenna numbering
    anttrips = []
    for i in ants[amin:amax+1]:
        for j in ants[list(ants).index(i)+1:amax+1]:
            for k in ants[list(ants).index(j)+1:amax+1]:
                anttrips.append([i,j,k])
                
    # next return data indexes for triples
    bltrips = []
    for (ant1, ant2, ant3) in anttrips:
        try:
            bl1 = np.where( (blarr[:,0] == ant1) & (blarr[:,1] == ant2) )[0][0]
            bl2 = np.where( (blarr[:,0] == ant2) & (blarr[:,1] == ant3) )[0][0]
            bl3 = np.where( (blarr[:,0] == ant1) & (blarr[:,1] == ant3) )[0][0]
            bltrips.append([bl1, bl2, bl3])
        except IndexError:
            continue

    return bltrips

def phaseshift(data0, d, l1, m1, u, v, verbose=0):
    """ Shift phase center
    """

    ang = lambda dl,dm,u,v,freq: (dl*np.outer(u,freq/d['freq_orig'][0]) + dm*np.outer(v,freq/d['freq_orig'][0]))
#    ang = lambda dl,dm,u,v,freq: (dl*np.outer(u,freq/d['freq'][0]) + dm*np.outer(v,freq/d['freq'][0]))
    freq = d['freq']
    dl = l1 - d['l0']
    dm = m1 - d['m0']

    if (dl != 0.) or (dm != 0.):
        frot = ang(dl, dm, u[len(u)/2], v[len(v)/2], freq)
        for i in xrange(data0.shape[0]):    # iterate over integrations
#            frot = ang(dl, dm, u[i], v[i], freq)
            for pol in xrange(data0.shape[3]):    # iterate over pols
                data0[i,:,:,pol] = data0[i,:,:,pol] * np.exp(-2j*np.pi*frot)
    else:
        if verbose:
            print 'No phase rotation needed'

    # new phase center
    d['l0'] = l1
    d['m0'] = m1

def calc_delay(d, dm):
    """ Function to calculate delay for each channel, given dm
    """

    return np.round((4.2e-3 * dm * (d['freq']**(-2) - d['freq'][len(d['chans'])-1]**(-2)))/d['inttime'],0).astype(int)

def dedisperse(data0, d, dm, verbose=0):
    """ dedisperse the data in place
    """

    # calc relative delay per channel. only shift minimally
    newdelay = calc_delay(d, dm)
    relativedelay = newdelay - d['datadelay']

    for i in xrange(len(d['chans'])):
        if relativedelay[i] > 0:
            data0[:,:,i,:] = np.roll(data0[:,:,i,:], -relativedelay[i], axis=0)
            data0[-relativedelay[i]:,:,i,:] = 0j   # flag meaningless roll data
    if verbose:
        print 'Dedispersed for DM=%d' % dm

    # new delay values
    d['datadelay'] = newdelay

def make_bispectra(data0, triples):
    """ Makes bispectra by integration, which assumes data have been dedispersed.
    """

    bisp = lambda data0: data0[:,:,0] * data0[:,:,1] * np.conj(data0[:,:,2])    # bispectrum for data referenced by triple (data[:,triples])

    bispectra = bisp(data0.mean(axis=2).mean(axis=2)[:, triples])
    return bispectra
