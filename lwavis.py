#! /usr/bin/env python

"""lwavis.py - visualization of lwa visibilities
modified from evlavis.py
claw, 9 april 2012
"""

import string
from os.path import join
import os
import numpy as n
import pylab as p
from casa import ms
from casa import quanta as qa


def sigma_clip(arr,sigma=3):
    """Function takes 1d array of values and returns the sigma-clipped min and max scaled by value "sigma".
    Useful for clipping arrays with large outliers (like RFI). Used below to make plots a bit easier to interpret.
    """

    cliparr = range(len(arr))  # initialize
    arr = n.append(arr,[1])    # append superfluous item to trigger loop
    while len(cliparr) != len(arr):
        arr = arr[cliparr]
        mean = arr.mean()
        std = arr.std()
        cliparr = n.where((arr < mean + sigma*std) & (arr > mean - sigma*std))[0]
    return mean - sigma*std, mean + sigma*std


class lwa:
    def __init__(self, directory, nints=1, nskip=0, ddid=-1, selectpol=['XX','YY']):
        """Initializes the class "lwa". This creates new object containing data and metadata for a set of files in a directory.
        It also includes functions to manipulate data and do some analysis, like making lightcurves, etc.
        Note that this uses CASA libraries in a way that requires it to be run from within "casapy".
        Default is to read in first file (alphabetically) in the directory. Use nints and nskip to control where to start and number of files to read. this assumes that the alphabetical order is the time order.

        Examples of usage in python/casapy:
        import lwavis
        obs = lwavis.lwa('directory') -- create observation object for first file in a directory.
        print obs.data.shape -- see the structure of data read in. dimensions are (time, baseline, channel, polarization)
        results = obs.bisplc(show=1) -- create a bispectrum lightcurve and show any candidate transients. results are returned in 'return' object
        """

        # critical parameters. need to edit these
        ants = range(256)    # set what antennas to use. default is to use "range" to specify all antennas with range(n_ant+1)
        self.chans = n.array(range(8))  # set what channesl to use. default is to use range to select all channels.
        self.track0 = [n.zeros(len(self.chans)), self.chans]

        # read in all ms files in given directory, according to nints and nskip. 
        # assumes one integration per file
        files = os.listdir(directory)
        msfiles = []
        for f in files:
            if f[-2:] == 'ms':
                msfiles.append(f)
        msfiles = msfiles[nskip:nskip+nints]

        # open first file and read a bit of data to define data structure
        self.file = msfiles[0]
        print 'Reading ', self.file
        ms.open(self.file)
        spwinfo = ms.getspectralwindowinfo()
        summary = ms.summary()

        # read in multiple subbands ("data id" in casa parlance). lwa probably doesn't use this.
        if ddid < 0:
            ddidlist = range(len(spwinfo['spwInfo']))
        else:
            ddidlist = [ddid]

        freq = n.array([])
        for ddid in ddidlist:
            nch = spwinfo['spwInfo'][str(ddid)]['NumChan']
            ch0 = spwinfo['spwInfo'][str(ddid)]['Chan1Freq']
            chw = spwinfo['spwInfo'][str(ddid)]['ChanWidth']
            freq = n.concatenate( (freq, (ch0 + chw * n.arange(nch)) * 1e-9) )

        self.freq = freq[self.chans]
        self.nchan = len(self.freq)

        # read data into data structure. start with subband 0, then iterate over higher ones. lwa probably only has 0.
        ms.selectinit(datadescid=0)  # reset select params for later data selection
        selection = {'antenna1': ants, 'antenna2': ants}
        ms.select(items = selection)
        print 'Reading SB %d, polarization %s...' % (0, selectpol)
        ms.selectpolarization(selectpol)
        da = ms.getdata(['data','axis_info'], ifraxis=True)
        if da == {}:
            print 'No data found.'
            return 1
        newda = n.transpose(da['data'], axes=[3,2,1,0])  # if using multi-pol data.
        if len(ddidlist) > 1:
            for ddid in ddidlist[1:]:
                ms.selectinit(datadescid=ddid)  # reset select params for later data selection
                ms.select(items = selection)
                print 'Reading SB %d, polarization %s...' % (ddid, selectpol)
                ms.selectpolarization(selectpol)
                da = ms.getdata(['data','axis_info'], ifraxis=True)
                newda = n.concatenate( (newda, n.transpose(da['data'], axes=[3,2,1,0])), axis=2 )
        ms.close()
        rawdata = newda  # array for collecting raw data

        # check pol and baseline dimensions of data
        self.npol_orig = da['data'].shape[0]
        self.npol = len(selectpol)
        self.nbl = da['data'].shape[2]
        print 'Initializing %d of %d polarizations' % (self.npol, self.npol_orig)
        print 'Initializing nchan:', self.nchan
        print 'Initializing nbl:', self.nbl
        self.nskip = int(nskip*self.nbl)    # number of iterations to skip (for reading in different parts of buffer)

        # set number of antennas and names of baselines
        bls = da['axis_info']['ifr_axis']['ifr_shortname']
        self.blarr = n.array([[bls[i].split('-')[0],bls[i].split('-')[1]] for i in range(len(bls))])
        self.ants = n.unique(self.blarr)
        self.nants = len(self.ants)
        print 'Initializing nants:', self.nants

        # set integration time and time axis
        ti = da['axis_info']['time_axis']['MJDseconds']
        self.inttime = n.mean([ti[i+1] - ti[i] for i in range(len(ti)-1)])
        print 'Initializing integration time (s):', self.inttime
        self.reltime = ti - ti[0]

        if len(msfiles) > 1:
            for file in msfiles[1:]:

                ms.open(file)
                print 'Reading ', file
                spwinfo = ms.getspectralwindowinfo()

                # read in multiple subbands ("data id" in casa parlance). lwa probably doesn't use this.
                if ddid < 0:
                    ddidlist = range(len(spwinfo['spwInfo']))
                else:
                    ddidlist = [ddid]

                # read data into data structure. start with subband 0, then iterate over higher ones. lwa probably only has 0.
                ms.selectinit(datadescid=0)  # reset select params for later data selection
                ms.select(items = selection)
                print 'Reading SB %d, polarization %s...' % (0, selectpol)
                ms.selectpolarization(selectpol)
                da = ms.getdata(['data','axis_info'], ifraxis=True)
                if da == {}:
                    print 'No data found.'
                    return 1
                newda = n.transpose(da['data'], axes=[3,2,1,0])  # if using multi-pol data.
                if len(ddidlist) > 1:
                    for ddid in ddidlist[1:]:
                        ms.selectinit(datadescid=ddid)  # reset select params for later data selection
                        ms.select(items = selection)
                        print 'Reading SB %d, polarization %s...' % (ddid, selectpol)
                        ms.selectpolarization(selectpol)
                        da = ms.getdata(['data','axis_info'], ifraxis=True)
                        newda = n.concatenate( (newda, n.transpose(da['data'], axes=[3,2,1,0])), axis=2 )
                ms.close()
                rawdata = n.concatenate( (rawdata, newda), axis=0)

        # create data structures
        self.rawdata = rawdata
        self.data = rawdata[:,:,self.chans]   # remove channels ignored earlier
        self.dataph = (self.data.mean(axis=3).mean(axis=1)).real   # create dataph, which is sum over all baselines. akin to calculating tied-array beam (possibly without calibration)
        self.min = self.dataph.min()
        self.max = self.dataph.max()
        print 'Shape of rawdata, data:'
        print self.rawdata.shape, self.data.shape
        print 'Dataph min, max:'
        print self.min, self.max



    def spec(self, ind=[], save=0, pathout='./'):
        """Plot dynamic spectrum of data by summing all cross correlations and polarizations.
        This is akin to calculating tied-array beam (without calibration if the data are uncalibrated).
        A spectrum like this is useful for checking data quality and RFI levels.
        """

        reltime = self.reltime

        abs = self.dataph
        print 'Data mean, std: %f, %f' % (self.dataph.mean(), self.dataph.std())
        (vmin, vmax) = sigma_clip(abs.ravel())

        p.figure(1)
        p.clf()
        ax = p.axes()
        ax.set_position([0.2,0.2,0.7,0.7])
        if len(ind) > 0:
            for i in ind:
                p.subplot(len(ind),1,list(ind).index(i)+1)
                intmin = n.max([0,i-50])
                intmax = n.min([len(self.reltime),i+50])
                im = p.imshow(n.rot90(abs[intmin:intmax]), aspect='auto', origin='upper', interpolation='nearest', extent=(intmin,intmax,0,len(self.chans)), vmin=vmin, vmax=vmax)
                p.subplot(len(ind),1,1)
        else:
            im = p.imshow(n.rot90(abs), aspect='auto', origin='upper', interpolation='nearest', extent=(0,len(self.reltime),0,len(self.chans)), vmin=vmin, vmax=vmax)
        p.title(str(self.nskip/self.nbl) + ' nskip, candidates ' + str(ind))
        cb = p.colorbar(im)
        cb.set_label('Flux Density (Jy)',fontsize=12,fontweight="bold")
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_position(('outward', 20))
        ax.spines['left'].set_position(('outward', 30))
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        p.yticks(n.arange(0,len(self.chans),4), (self.chans[(n.arange(0,len(self.chans), 4))]))
        p.xlabel('Time (integration number)',fontsize=12,fontweight="bold")
        p.ylabel('Frequency Channel',fontsize=12,fontweight="bold")
        if save:
            savename = self.file.split('.')[:-1]
            savename.append(str(self.nskip/self.nbl) + '.spec.png')
            savename = string.join(savename,'.')
            print 'Saving file as ', savename
            p.savefig(pathout+savename)


    def tracksub(self, tbin, bgwindow = 0):
        """Creates a background-subtracted set of visibilities.
        For a given track (i.e., an integration number) and bg window, tracksub subtractes a background in time and returns an array with new data.
        """

        data = self.data
        track_t,track_c = self.track0  # get track time and channel arrays
        trackon = (list(n.array(track_t)+tbin), track_c)   # create new track during integration of interest
        dataon = data[trackon[0], :, trackon[1]]

        # set up bg track
        if bgwindow:
            # measure max width of pulse (to avoid in bgsub)
            bgrange = range(-bgwindow/2 + tbin, tbin) + range(tbin, bgwindow/2 + tbin)
            for k in bgrange:     # build up super track for background subtraction
                if bgrange.index(k) == 0:   # first time through
                    trackoff = (list(n.array(track_t)+k), track_c)
                else:    # then extend arrays by next iterations
                    trackoff = (trackoff[0] + list(n.array(track_t)+k), trackoff[1] + track_c)

            dataoff = data[trackoff[0], :, trackoff[1]]

        datadiffarr = n.zeros((self.nchan, self.nbl, self.npol),dtype='complex')

        # compress time axis, then subtract on and off tracks
        for ch in n.unique(trackon[1]):
            indon = n.where(trackon[1] == ch)

            if bgwindow:
                indoff = n.where(trackoff[1] == ch)
                datadiffarr[ch] = dataon[indon].mean(axis=0) - dataoff[indoff].mean(axis=0)
            else:
                datadiffarr[ch] = dataon[indon].mean(axis=0)

        return n.transpose(datadiffarr, axes=[2,1,0])


    def tripgen(self, amin=0, amax=0):
        """Calculates and returns data indexes (i,j,k) for all closed triples.
        Used for calculating the bispectrum.
        amin and amax define range of antennas (with index, in order). only used if nonzero.
        """

        if amax == 0:
            amax = self.nants
        blarr = self.blarr

        # first make triples indexes in antenna numbering
        anttrips = []

        for i in self.ants[amin:amax+1]:
            for j in self.ants[list(self.ants).index(i)+1:amax+1]:
                for k in self.ants[list(self.ants).index(j)+1:amax+1]:
                    anttrips.append([i,j,k])

        # next return data indexes for triples
        bltrips = []
        for (ant1, ant2, ant3) in anttrips:
            try:
                bl1 = n.where( (blarr[:,0] == ant1) & (blarr[:,1] == ant2) )[0][0]
                bl2 = n.where( (blarr[:,0] == ant2) & (blarr[:,1] == ant3) )[0][0]
                bl3 = n.where( (blarr[:,0] == ant1) & (blarr[:,1] == ant3) )[0][0]
                bltrips.append([bl1, bl2, bl3])
            except IndexError:
                continue

        return n.array(bltrips)

    def bisplc(self, chan=-1, bgwindow=4, show=0, save=0, sigma=5., pathout='./'):
        """Steps in Bispectrum Transient Detection Algorithm

1) Collect visibility spectra for some length of time. In Python, I read data into an array with a shape of (n_int, n_chan, n_bl).

2) Prepare visibility spectra to create bispectra. Optionally, one can form dedispersed spectra for each baseline. A simpler start (and probably more relevant for LOFAR) would be to instead select a single integration from the data array described above. Either way, this step changes the data shape to (n_chan, n_bl).

3) Subtract visibilities in time. If the sky has many (or complex) sources, the bispectrum is hard to interpret. Subtracting neighboring visibilities in time (or a rolling mean, like v_t2 - (v_t1+v_t3)/2) removes most constant emission. The only trick is that this assumes that the array has not rotated much and that gain and other effects have not changed. This should preserve the data shape as (n_chan, n_bl).

4) Calculate mean visibility for each baseline. After subtracting in time, one can measure the mean visibility across the band. This reduces the shape to (n_bl).

5) Form a bispectrum for every closed triple in the array. There are a total of n_a * (n_a-1) * (n_a-2) / 6 possible closed triples in the array, where n_a is the number of antennas. One way to form all bispectra is to iterate over antenna indices like this:
for i in range(0, len(n_a)-2):
  for j in range(i, len(n_a)-1):
    for k in range(k, len(n_a)):
      bl1, bl2, bl3 = ant2bl(i, j, k)
      bisp = vis[bl1] * vis[bl2] * vis[bl3]

As you can see, this loop needs a function to convert antenna triples to baseline triples (I call it "ant2bl" here). That is, for antennas (i, j, k), you need (bl_ij, bl_jk, bl_ki). Note that the order of the last baseline is flipped; this is a way of showing that the way you "close" a loop is by tracing a single line around all three baselines. This step changes the basic data product from a shape of (n_bl) to (n_tr). 

6) Search the set of bispectra for sign of a source. Each bispectrum is complex, but if there is a point source in the (differenced) data, all bispectra will respond in the same way. This happens regardless of the location in the field of view.
The mean of all bispectra will scale with the source brightness to the third power, since it is formed from the product of three visibilities. Oddly, the standard deviation of the bispectra will *also* change with the source brightness, due to something called "self noise". The standard deviation of bispectra in the real-imaginary plane should be sqrt(3) S^2 sigma_bl, where S is the source brightness and sigma_bl is the noise on an individual baseline.
In practice, this search involves plotting the mean bispectrum versus time and searching for large deviations. At the same time, a plot of mean versus standard deviation of bispectra will show whether any significant deviation obeys the expected self-noise scaling. That scaling is only valid for a single point source in the field of view, which is what you expect for a fast transient. Any other behavior would be either noise-like or caused by RFI. In particular, RFI will look like a transient, but since it does not often look like a point source, it can be rejected in the plot of mean vs. standard deviation of bispectra. This is a point that I've demonstrated on a small scale, but would needs more testing, since RFI is so varied.
        """

        if chan == -1:
            chans = self.chans
        else:
            chans = n.array([chan])

        bisp = lambda d,i,j,k: d[:,i] * d[:,j] * n.conj(d[:,k])    # bispectrum for pol data
#        bisp = lambda d,i,j,k: n.complex(d[i] * d[j] * n.conj(d[k]))

        triples = self.tripgen()

        dibi = n.zeros((len(self.data)-bgwindow, len(triples)), dtype='complex')
        for ii in range(bgwindow+1, len(self.data)-(bgwindow+1)):
            diff = self.tracksub(ii, bgwindow=bgwindow)
            if len(n.shape(diff)) == 1:    # no track
                continue
            diffmean = diff.mean(axis=2)

            for trip in range(len(triples)):
                i, j, k = triples[trip]
                dibi[ii, trip] = bisp(diffmean, i, j, k).mean(axis=0)  # Stokes I bispectrum. should be ok for multipol data...

        if show or save:
            good = self.det_bisplc(dibi, save=save, show=show, sigma=sigma, pathout=pathout)
            return good
        else:
            return dibi


    def det_bisplc(self, dibi, sigma=5., tol=1.3, show=0, save=0, pathout='./'):
        """Function to search for a transient in a bispectrum lightcurve.
        Designed to be used by bisplc function or easily fed the output of that function.
        sigma gives the threshold for SNR_bisp (apparent). 
        tol gives the amount of tolerance in the sigma_b cut for point-like sources (rfi filter).
        Returns the SNR and integration number of any candidate events.
        """

        ntr = lambda num: num*(num-1)*(num-2)/6.

        # using s=S/Q
#        mu = lambda s: s/(1+s)  # for independent bispectra, as in kulkarni 1989
        mu = 1.  # for bispectra from visibilities
        sigbQ3 = lambda s: n.sqrt((1 + 3*mu**2) + 3*(1 + mu**2)*s**2 + 3*s**4)  # from kulkarni 1989, normalized by Q**3
        s = lambda dibisnr, num: (dibisnr/n.sqrt(ntr(num)))**(1/3.)

        dibimean = dibi.real.mean(axis=1)
        dibistd = dibi.real.std(axis=1)
        (meanmin,meanmax) = sigma_clip(dibimean)  # remove rfi
        (stdmin,stdmax) = sigma_clip(dibistd)  # remove rfi
        clipped = n.where((dibimean > meanmin) & (dibimean < meanmax) & (dibistd > stdmin) & (dibistd < stdmax) & (dibimean != 0.0))[0]  # remove rf
        dibimeanstd = dibi[clipped].real.mean(axis=1).std()
        dibisnr = dibimean/dibimeanstd
        Q = (dibimeanstd*n.sqrt(ntr(self.nants)))**(1/3.)
#        Q = dibimeanstd**(1/3.)    # this works for b0329 data, too, but includes sqrt(ntr) term
#        Q = n.median( dibistd[clipped]**(1/3.) )              # alternate for Q
#        Q = sigt0toQ(dibimeanstd, self.nants)              # alternate for Q
        print 'Noise per baseline (system units), Q =', Q

        # detect
        cands = n.where( (dibistd/Q**3 < tol*sigbQ3(s(dibisnr, self.nants))) & (dibisnr > sigma) )[0]  # define compact sources with good snr

        # plot snrb lc and expected snr vs. sigb relation
        if show or save:
            p.figure(1)
            ax = p.axes()
            p.subplot(211)
            p.title(str(self.nskip/self.nbl)+' nskip, ' + str(len(cands))+' candidates', transform = ax.transAxes)
            p.plot(dibisnr, 'b.')
            if len(cands) > 0:
                p.plot(cands, dibisnr[cands], 'r*')
                p.ylim(-2*dibisnr[cands].max(),2*dibisnr[cands].max())
            p.xlabel('Integration')
            p.ylabel('SNR$_{bisp}$')
            p.subplot(212)
            p.plot(dibistd/Q**3, dibisnr, 'b.')

            # plot reference theory lines
            smax = s(dibisnr.max(), self.nants)
            sarr = smax*n.arange(0,21)/20.
            p.plot(sigbQ3(sarr), sarr**3*n.sqrt(ntr(self.nants)), 'k')
            p.plot(tol*sigbQ3(sarr), sarr**3*n.sqrt(ntr(self.nants)), 'k--')
            p.plot(dibistd[cands]/Q**3, dibisnr[cands], 'r*')

            if len(cands) > 0:
                p.axis([0, tol*sigbQ3(s(dibisnr[cands].max(), self.nants)), -0.5*dibisnr[cands].max(), 1.1*dibisnr[cands].max()])

                # show spectral modulation next to each point
                for candint in cands:
                    sm = n.single(round(self.specmod(candint),1))
                    p.text(dibistd[candint]/Q**3, dibisnr[candint], str(sm), horizontalalignment='right', verticalalignment='bottom')
            p.xlabel('$\sigma_b/Q^3$')
            p.ylabel('SNR$_{bisp}$')
            if save:
                savename = self.file.split('.')[:-1]
                savename.append(str(self.nskip/self.nbl) + '_' + '.bisplc.png')
                savename = string.join(savename,'.')
                p.savefig(pathout+savename)
            else:
                pass

        return dibisnr[cands], dibistd[cands], cands


    def specmod(self, tbin, bgwindow=4):
        """Calculate spectral modulation for given track.
        Spectral modulation is basically the standard deviation of a spectrum. 
        This helps quantify whether the flux is located in a narrow number of channels or across all channels.
        Narrow RFI has large (>5) modulation, while spectrally broad emission has low modulation.
        See Spitler et al 2012 for details.
        """

        diff = self.tracksub(tbin, bgwindow=bgwindow)
        bfspec = diff.mean(axis=0).real  # should be ok for multipol data...
        sm = n.sqrt( ((bfspec**2).mean() - bfspec.mean()**2) / bfspec.mean()**2 )

        return sm


def pulse_search_bisp(fileroot, pathin, pathout, nints=600, edge=15):
    """Blind search for pulses with bispectrum algorithm.
    May be useful for running the functions above on a large data file and producing many plots automatically.
    If lwa data files have only one integration, this needs to be modified a bit.
    """

    startint = 0
    maxints = 10000
    bgwindow = 4
    filelist = [fileroot]

    for file in filelist:
        for nskip in range(startint, maxints-(nints-edge), nints-edge):
            print
            print 'Starting file %s with nskip %d' % (file, nskip)

            lwa = lwa(pathin + file, nints=nints, nskip=nskip, selectpol=['XX'])
            if ev.min == ev.max: break   # escape if no data
            (dibimean, dibistd, dibi) = ev.bisplc(bgwindow=bgwindow)
            dets = ev.det_bisplc(dibimean, dibistd, sigma=5, save=1, pathout=pathout)

            print 'Candidates:'
            print dets
            
            lwa.spec(n.unique(indagg), save=1, pathout=pathout)
