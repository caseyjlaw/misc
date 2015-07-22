
#! /usr/bin/env python

"""= evlavis.py - visualization of evla fast dump data
& claw
: Unknown
+
 Script to load and visualize evla data
--
"""

import sys, string, os, shutil
#import cProfile
from os.path import join
import pickle
import numpy as n
import pylab as p
import scipy.optimize as opt
import scipy.stats.morestats as morestats
try:
    import aipy
except ImportError:
    print 'No aipy available. Can\'t image ms data.'
#from threading import Thread
#from matplotlib.font_manager import fontManager, FontProperties
if len(sys.argv) != 2:
    try:
        from casa import ms
        from casa import quanta as qa
    except:
        from mirtask import util
        from mirexec import TaskInvert, TaskClean, TaskRestore, TaskImFit, TaskCgDisp, TaskImStat, TaskUVFit
        import miriad
else:
    from mirtask import util
    from mirexec import TaskInvert, TaskClean, TaskRestore, TaskImFit, TaskCgDisp, TaskImStat, TaskUVFit
    import miriad


def sigma_clip(arr,sigma=3):
    """Function takes 1d array of values and returns the sigma-clipped min and max scaled by value "sigma".
    """

    cliparr = range(len(arr))  # initialize
    arr = n.append(arr,[1])    # append superfluous item to trigger loop
    while len(cliparr) != len(arr):
        arr = arr[cliparr]
        mean = arr.mean()
        std = arr.std()
        cliparr = n.where((arr < mean + sigma*std) & (arr > mean - sigma*std))[0]
    return mean - sigma*std, mean + sigma*std


class evla:
    def __init__(self, file, nints=1000, nskip=0, nocal=False, nopass=False, ddid=[-1], selectpol=['RR','LL'], scan=0, datacol='data'):
        """Initializes the class according to the file type. Miriad or Measurement Set?
        Note that CASA stuff only works if this is run from within casapy.
        Scan is zero-based selection based on scan order, not actual scan number.
        """

        # initialize
        self.approxuvw = True      # flag to make template visibility file to speed up writing of dm track data
#        self.dmarr =  [30.,56.8,90.] # [56.8]  # crab
#        self.pulsewidth = 0.006  # pulse width of crab and m31 candidates. later turned into array of len(chans)
        self.dmarr =  [44,88.] # j0628+09
        self.pulsewidth = 0.0
        self.file = file
        self.scan = scan

        if file.split('.')[-1][:3] == 'mir':
            self.data_type = 'mir'
            print 'Loading miriad file.'
            self.miriad_init(file, nints=nints, nskip=nskip, nocal=nocal, nopass=nopass)
        elif file.split('.')[-1][:2] == 'ms':
            self.data_type = 'ms'
            print 'Loading ms file.'
            status = self.casa_init(file, nints=nints, nskip=nskip, ddid=ddid, selectpol=selectpol, scan=scan, datacol=datacol)
            if status:
                print 'Stopping init...'
                return

        # set up ur tracks
        self.dmtrack0 = {}
        self.twidths = {}
        for dmbin in range(len(self.dmarr)):
            self.dmtrack0[dmbin] = self.dmtrack(self.dmarr[dmbin],0)
            self.twidths[dmbin] = 0
            for k in self.dmtrack0[dmbin][1]:
                self.twidths[dmbin] = max(self.twidths[dmbin], len(n.where(n.array(self.dmtrack0[dmbin][1]) == k)[0]))


    def miriad_init(self, file, nints, nskip, nocal, nopass):
        """Reads in Miriad data using miriad-python.
        Seems to have some small time (~1 integration) errors in light curves and spectrograms, 
        as compared to CASA-read data.
        """

        vis = miriad.VisData(self.file,)

#        li = range(1,13) + range(14,23) + range(24,62)
        li = range(64)
        self.chans = n.array(li)

        # read data into python arrays
        i = 0
        for inp, preamble, data, flags in vis.readLowlevel ('dsl3', False, nocal=True, nopass=True):
            # Loop to skip some data and read shifted data into original data arrays
            if i == 0:
                # get few general variables
                self.nants0 = inp.getScalar ('nants', 0)
                self.inttime0 = inp.getScalar ('inttime', 10.0)
                self.nspect0 = inp.getScalar ('nspect', 0)
                self.nwide0 = inp.getScalar ('nwide', 0)
                self.sdf0 = inp.getScalar ('sdf', self.nspect0)
                self.nschan0 = inp.getScalar ('nschan', self.nspect0)
                self.ischan0 = inp.getScalar ('ischan', self.nspect0)
                self.sfreq0 = inp.getScalar ('sfreq', self.nspect0)
                self.restfreq0 = inp.getScalar ('restfreq', self.nspect0)
                self.pol0 = inp.getScalar ('pol')

                self.sfreq = self.sfreq0
                self.sdf = self.sdf0
                self.npol = 1
                self.nchan = len(data)
                print 'Initializing nchan:', self.nchan
                bls = []

            # build complete list of baselines
            bls.append(preamble[4])

            # end here. assume at least one instance of each bl occurs before ~three integrations
            if len(bls) == 3*len(n.unique(bls)):
                blarr = []
                for bl in n.unique(bls):
                    blarr.append(util.decodeBaseline (bl))
                self.blarr = n.array(blarr)
                bldict = dict( zip(n.unique(bls), n.arange(len(blarr))) )
                break

            i = i+1

        # Initialize more stuff...
        self.freq = self.sfreq + self.sdf * self.chans
        self.pulsewidth = self.pulsewidth * n.ones(len(self.chans)) # pulse width of crab and m31 candidates

        # good baselines
        self.nbl = len(self.blarr)
        print 'Initializing nbl:', self.nbl
        self.ants = n.unique(self.blarr)
        self.nants = len(self.ants)
        print 'Initializing nants:', self.nants
        self.nskip = int(nskip*self.nbl)    # number of iterations to skip (for reading in different parts of buffer)
        nskip = int(self.nskip)

        # define data arrays
        da = n.zeros((nints,self.nbl,self.nchan),dtype='complex64')
        fl = n.zeros((nints,self.nbl,self.nchan),dtype='bool')
        pr = n.zeros((nints*self.nbl,5),dtype='float64')

        print
        # go back and read data into arrays
        i = 0
        for inp, preamble, data, flags in vis.readLowlevel ('dsl3', False, nocal=nocal, nopass=nopass):
            # Loop to skip some data and read shifted data into original data arrays

            if i < nskip:
                i = i+1
                continue 

            # assumes ints in order, but may skip. after nbl iterations, it fills next row, regardless of number filled.
            if (i-nskip) < nints*self.nbl:
                da[(i-nskip)//self.nbl,bldict[preamble[4]]] = data
                fl[(i-nskip)//self.nbl,bldict[preamble[4]]] = flags
                pr[i-nskip] = preamble
            else:
                break     # stop at nints

            if not (i % (self.nbl*100)):
                print 'Read spectrum ', str(i)

            i = i+1

        # build final data structures
#        good = n.where ( (self.blarr[:,0] != 5) & (self.blarr[:,1] != 5) & (self.blarr[:,0] != 10) & (self.blarr[:,1] != 10) )[0] # remove bad ants?
        self.rawdata = n.expand_dims(da, 3)  # hack to get superfluous pol axis
        self.flags = n.expand_dims(fl, 3)
        self.data = self.rawdata[:,:,self.chans,:] # [:,good,:,:]  # remove bad ants?
#        self.blarr = self.blarr[good]  # remove bad ants?
        self.preamble = pr
        self.u = (pr[:,0] * self.freq.mean()).reshape(nints, self.nbl)
        self.v = (pr[:,1] * self.freq.mean()).reshape(nints, self.nbl)
        self.w = (pr[:,2] * self.freq.mean()).reshape(nints, self.nbl)
        # could add uvw, too... preamble index 0,1,2 in units of ns
        self.dataph = (self.data.mean(axis=3).mean(axis=1)).real  #dataph is summed and detected to form TP beam at phase center, multi-pol
        time = self.preamble[::self.nbl,3]
        self.reltime = 24*3600*(time - time[0])      # relative time array in seconds. evla times change...?

        # print summary info
        print
        print 'Data read!\n'
        print 'Shape of raw data, flagged, time:'
        print self.rawdata.shape, self.data.shape, self.reltime.shape
        self.min = self.dataph.min()
        self.max = self.dataph.max()
        print 'Dataph min, max:'
        print self.min, self.max


    def casa_init(self, file, nints=1000, nskip=0, ddid=[-1], selectpol=['RR','LL'], scan=0, datacol='data'):
        """Reads in Measurement Set data using CASA.
        ddid is list of subbands. zero-based.
        Scan is zero-based selection based on scan order, not actual scan number.
        """

        # open file and read a bit of data to define data structure
        ants = range(0,28)
        self.dmtrack0 = 0  # initialize ur-track
  
        # get spw info. either load pickled version (if found) or make new one
        pklname = string.join(file.split('.')[:-1], '.') + '_init.pkl'
#        pklname = pklname.split('/')[-1]  # hack to remove path and write locally
        if os.path.exists(pklname):
            print 'Pickle of initializing info found. Loading...'
            pkl = open(pklname, 'r')
            try:
                (self.npol_orig, self.npol, self.nbl, self.blarr, self.ants, self.nants, self.nants0, self.inttime, self.inttime0, spwinfo, scansummary) = pickle.load(pkl)
            except EOFError:
                print 'Bad pickle file. Exiting...'
                return 1
            scanlist = scansummary['summary'].keys()
            starttime_mjd = scansummary['summary'][scanlist[scan]]['0']['BeginTime']
            self.nskip = int(nskip*self.nbl)    # number of iterations to skip (for reading in different parts of buffer)
            self.npol = len(selectpol)
        else:
            print 'No pickle of initializing info found. Making anew...'
            pkl = open(pklname, 'wb')
            ms.open(self.file)
            spwinfo = ms.getspectralwindowinfo()
            scansummary = ms.getscansummary()
            scanlist = scansummary['summary'].keys()

            starttime_mjd = scansummary['summary'][scanlist[scan]]['0']['BeginTime']
            starttime0 = qa.getvalue(qa.convert(qa.time(qa.quantity(starttime_mjd+0/(24.*60*60),'d'),form=['ymd'], prec=9), 's'))
            stoptime0 = qa.getvalue(qa.convert(qa.time(qa.quantity(starttime_mjd+0.5/(24.*60*60), 'd'), form=['ymd'], prec=9), 's'))
            ms.selectinit(datadescid=0)  # initialize to initialize params
            selection = {'time': [starttime0, stoptime0], 'antenna1': ants, 'antenna2': ants}
            ms.select(items = selection)
            da = ms.getdata([datacol,'axis_info'], ifraxis=True)
            ms.close()

            self.npol_orig = da[datacol].shape[0]
            self.npol = len(selectpol)
            self.nbl = da[datacol].shape[2]
            print 'Initializing %d of %d polarizations' % (self.npol, self.npol_orig)
            print 'Initializing nbl:', self.nbl

            # good baselines
            bls = da['axis_info']['ifr_axis']['ifr_shortname']
            self.blarr = n.array([[int(bls[i].split('-')[0]),int(bls[i].split('-')[1])] for i in range(len(bls))])
            self.ants = n.unique(self.blarr)
            self.nants = len(self.ants)
            self.nants0 = len(self.ants)
            print 'Initializing nants:', self.nants
            self.nskip = int(nskip*self.nbl)    # number of iterations to skip (for reading in different parts of buffer)

            # set integration time
            ti0 = da['axis_info']['time_axis']['MJDseconds']
#            self.inttime = n.mean([ti0[i+1] - ti0[i] for i in range(len(ti0)-1)])
            self.inttime = scansummary['summary'][scanlist[scan]]['0']['IntegrationTime']
            self.inttime0 = self.inttime
            print 'Initializing integration time (s):', self.inttime

            pickle.dump((self.npol_orig, self.npol, self.nbl, self.blarr, self.ants, self.nants, self.nants0, self.inttime, self.inttime0, spwinfo, scansummary), pkl)
        pkl.close()

        # set desired spw
        if (len(ddid) == 1) & (ddid[0] == -1):
            ddidlist = range(len(spwinfo['spwInfo']))
        else:
            ddidlist = ddid

        freq = n.array([])
        for ddid in ddidlist:
            nch = spwinfo['spwInfo'][str(ddid)]['NumChan']
            ch0 = spwinfo['spwInfo'][str(ddid)]['Chan1Freq']
            chw = spwinfo['spwInfo'][str(ddid)]['ChanWidth']
            freq = n.concatenate( (freq, (ch0 + chw * n.arange(nch)) * 1e-9) )

#        self.chans = n.array(range(2,62))  # can flag by omitting channels here
        self.chans = n.arange(nch*len(ddidlist))    # default is to take all chans
#        self.chans = self.chans[range(5,31)+range(32,41)+range(42,59)+range(69,123)]
        self.chans = self.chans[range(5,31)+range(32,41)+range(42,59)]
#        self.chans = self.chans[range(69,123)]
        self.freq = freq[self.chans]
        self.nchan = len(self.freq)
        print 'Initializing nchan:', self.nchan

        # set requested time range based on given parameters
        timeskip = self.inttime*nskip
        starttime = qa.getvalue(qa.convert(qa.time(qa.quantity(starttime_mjd+timeskip/(24.*60*60),'d'),form=['ymd'], prec=9), 's'))
        stoptime = qa.getvalue(qa.convert(qa.time(qa.quantity(starttime_mjd+(timeskip+nints*self.inttime)/(24.*60*60), 'd'), form=['ymd'], prec=9), 's'))
        print 'First integration of scan:', qa.time(qa.quantity(starttime_mjd,'d'),form=['ymd'],prec=9)
        print
        print 'Reading scan', str(scanlist[scan]) ,'for times', qa.time(qa.quantity(starttime_mjd+timeskip/(24.*60*60),'d'),form=['hms'], prec=9), 'to', qa.time(qa.quantity(starttime_mjd+(timeskip+nints*self.inttime)/(24.*60*60), 'd'), form=['hms'], prec=9)

        # read data into data structure
        ms.open(self.file)
        ms.selectinit(datadescid=ddidlist[0])  # reset select params for later data selection
        selection = {'time': [starttime, stoptime], 'antenna1': ants, 'antenna2': ants}
        ms.select(items = selection)
        print 'Reading %s column, SB %d, polarization %s...' % (datacol, ddidlist[0], selectpol)
        ms.selectpolarization(selectpol)
        da = ms.getdata([datacol,'axis_info','u','v','w'], ifraxis=True)
        u = da['u']; v = da['v']; w = da['w']
#        da = ms.getdata(['data','axis_info'], ifraxis=False)
        if da == {}:
            print 'No data found.'
            return 1
        newda = n.transpose(da[datacol], axes=[3,2,1,0])  # if using multi-pol data.
#        newda = n.transpose(da['data'], axes=[2,1,0])  # if using multi-pol data.
#        scale=2.3  # hack!
        if len(ddidlist) > 1:
            for ddid in ddidlist[1:]:
                ms.selectinit(datadescid=ddid)  # reset select params for later data selection
                ms.select(items = selection)
                print 'Reading %s column, SB %d, polarization %s...' % (datacol, ddid, selectpol)
                ms.selectpolarization(selectpol)
                da = ms.getdata([datacol,'axis_info'], ifraxis=True)
                newda = n.concatenate( (newda, n.transpose(da[datacol], axes=[3,2,1,0])), axis=2 )
#                da = ms.getdata(['data','axis_info'], ifraxis=False)
#                newda = n.concatenate( (newda, n.transpose(da['data'], axes=[2,1,0])), axis=1 )
        ms.close()

        # Initialize more stuff...
        self.nschan0 = self.nchan
        self.pulsewidth = self.pulsewidth * n.ones(self.nchan) # pulse width of crab and m31 candidates

        # set variables for later writing data **some hacks here**
        self.nspect0 = 1
        self.nwide0 = 0
        self.sdf0 = da['axis_info']['freq_axis']['resolution'][0][0] * 1e-9
        self.sdf = self.sdf0
        self.ischan0 = 1
        self.sfreq0 = da['axis_info']['freq_axis']['chan_freq'][0][0] * 1e-9
        self.sfreq = self.sfreq0
        self.restfreq0 = 0.0
        self.pol0 = -1 # assumes single pol?

#        self.rawdata = newda[len(newda)/2:]  # hack to remove autos
        self.u = u.transpose() * (-self.freq.mean()*1e9/3e8)  # uvw are in m on ground. scale by -wavelenth to get projected lamba uvw (as in miriad?)
        self.v = v.transpose() * (-self.freq.mean()*1e9/3e8)
        self.w = w.transpose() * (-self.freq.mean()*1e9/3e8)
        self.rawdata = newda
        self.data = self.rawdata[:,:,self.chans]
        self.dataph = (self.data.mean(axis=3).mean(axis=1)).real  # multi-pol
        self.min = self.dataph.min()
        self.max = self.dataph.max()
        print 'Shape of rawdata, data:'
        print self.rawdata.shape, self.data.shape
        print 'Dataph min, max:'
        print self.min, self.max

        # set integration time and time axis
        ti = da['axis_info']['time_axis']['MJDseconds']
        self.reltime = ti - ti[0]
#        self.reltime = self.reltime[len(self.reltime)/2:]    # hack to remove autos


    def spec(self, ind=[], save=0, pathout='./'):
        reltime = self.reltime

        abs = self.dataph
        print 'Data mean, std: %f, %f' % (self.dataph.mean(), self.dataph.std())
        (vmin, vmax) = sigma_clip(abs.ravel())

#        p.figure(1,figsize=(11,6),dpi=120)
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
        p.title(str(self.scan) + ' scan, ' +str(self.nskip/self.nbl) + ' nskip, candidates ' + str(ind))
#        cb = p.colorbar(im, use_gridspec=True)
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
            savename.append(str(self.scan) + '_' + str(self.nskip/self.nbl) + '.spec.png')
            savename = string.join(savename,'.')
            print 'Saving file as ', savename
            p.savefig(pathout+savename)


    def drops(self, chan=0, pol=0, show=1):
        """Displays info on missing baselines.
        """

        nints = float(len(self.reltime))
        bllen = []

        if self.data_type == 'mir':
            bls = self.preamble[:,4]
            for bl in n.unique(bls):
                bllen.append(n.shape(n.where(bls == bl))[1])
        elif self.data_type == 'ms':
            for i in range(len(self.blarr)):
                bllen.append(len(n.where(self.data[:,i,chan,pol] != 0.00)[0]))

        bllen = n.array(bllen)

        if show:
            p.clf()
            for i in range(self.nbl):
                p.text(self.blarr[i,0], self.blarr[i,1], s=str(100*(bllen[i]/nints - 1)), horizontalalignment='center', verticalalignment='center', fontsize=9)
            p.axis((0,29,0,29))
            p.plot([0,29],[0,29],'b--')
#            p.xticks(int(self.blarr[:,0]), self.blarr[:,0])
#            p.yticks(int(self.blarr[:,1]), self.blarr[:,1])
            p.xlabel('Ant 1')
            p.ylabel('Ant 2')
            p.title('Drop fraction for chan %d, pol %d' % (chan, pol))
#            p.show()

        return self.blarr,bllen


    def fitspec(self, obsrms=0, save=0, pol=0, pathout='./'):
        """
        Fits a powerlaw to the mean spectrum at the phase center.
        Returns fit parameters.
        """

#        logname = self.file.split('_')[0:3]
#        logname.append('fitsp.txt')
#        logname = string.join(logname,'_')
#        log = open(logname,'a')

        # estimage of vis rms per channel from spread in imag space at phase center
        if obsrms == 0:
            print 'estimating obsrms from imaginary part of data...'
#            obsrms = n.std((((self.data).mean(axis=1)).mean(axis=0)).imag)/n.sqrt(2)  # sqrt(2) scales it to an amplitude error. indep of signal.
            obsrms = n.std((((self.data).mean(axis=1)).mean(axis=0)).imag)      # std of imag part is std of real part
#        spec = n.abs((((self.data).mean(axis=1))).mean(axis=0))
        spec = (((((self.data).mean(axis=3)).mean(axis=1))).mean(axis=0)).real
        print 'obsrms = %.2f' % (obsrms)

        plaw = lambda a, b, x: a * (x/x[0]) ** b
#        fitfunc = lambda p, x, rms:  n.sqrt(plaw(p[0], p[1], x)**2 + rms**2)   # for ricean-biased amplitudes
#        errfunc = lambda p, x, y, rms: ((y - fitfunc(p, x, rms))/rms)**2
        fitfunc = lambda p, x:  plaw(p[0], p[1], x)              # for real part of data
        errfunc = lambda p, x, y, rms: ((y - fitfunc(p, x))/rms)**2

        p0 = [100.,-5.]
        p1, success = opt.leastsq(errfunc, p0[:], args = (self.freq, spec, obsrms))
        print 'Fit results: ', p1
        chisq = errfunc(p1, self.freq, spec, obsrms).sum()/(self.nchan - 2)
        print 'Reduced chisq: ', chisq

        p.figure(2)
        p.errorbar(self.freq, spec, yerr=obsrms*n.ones(len(spec)), fmt='.')
        p.plot(self.freq, fitfunc(p1, self.freq), label='Fit: %.1f, %.2f. Noise: %.1f, chisq: %.1f' % (p1[0], p1[1], obsrms, chisq))
        p.xlabel('Frequency')
        p.ylabel('Flux Density (Jy)')
        p.legend()
        if save == 1:
            savename = self.file.split('.')[:-1]
            savename.append(str(self.nskip/self.nbl) + '.fitsp.png')
            savename = string.join(savename,'.')
            print 'Saving file as ', savename
            p.savefig(pathout+savename)
#            print >> log, savename, 'Fit results: ', p1, '. obsrms: ', obsrms, '. $\Chi^2$: ', chisq
        else:
            pass
#            print >> log, self.file, 'Fit results: ', p1, '. obsrms: ', obsrms, '. $\Chi^2$: ', chisq


    def dmtrack(self, dm = 0., t0 = 0., show=0):
        """Takes dispersion measure in pc/cm3 and time offset from first integration in seconds.
        t0 defined at first (unflagged) channel. Need to correct by flight time from there to freq=0 for true time.
        Returns an array of (timebin, channel) to select from the data array.
        """

        reltime = self.reltime
        chans = self.chans
        tint = 2.0*(self.reltime[1] - self.reltime[0])
#        tint = self.inttime0  # could do this instead...?

        # given freq, dm, dfreq, calculate pulse time and duration
        pulset_firstchan = 4.2e-3 * dm * self.freq[len(self.chans)-1]**(-2)   # used to start dmtrack at highest-freq unflagged channel
        pulset = 4.2e-3 * dm * self.freq**(-2) + t0 - pulset_firstchan  # time in seconds
        pulsedt = n.sqrt( (8.3e-6 * dm * (1000*self.sdf) * self.freq**(-3))**2 + self.pulsewidth**2)   # dtime in seconds

        timebin = []
        chanbin = []

        for ch in range(len(chans)):
            ontime = n.where(((pulset[ch] + pulsedt[ch]/2.) >= reltime - tint/2.) & ((pulset[ch] - pulsedt[ch]/2.) <= reltime + tint/2.))
            timebin = n.concatenate((timebin, ontime[0]))
            chanbin = n.concatenate((chanbin, (ch * n.ones(len(ontime[0]), dtype=int))))

        track = (list(timebin), list(chanbin))
#        print 'timebin, chanbin:  ', timebin, chanbin

        if show:
#            p.plot(track[1], track[0])
            p.plot(track[0], track[1], 'w*')

        return track


    def tracksub(self, dmbin, tbin, bgwindow = 0):
        """Reads data along dmtrack and optionally subtracts background like writetrack method.
        Returns the difference of the data in the on and off tracks as a single integration with all bl and chans.
        Nominally identical to writetrack, but gives visibilities values off at the 0.01 (absolute) level. Good enough for now.
        """

        data = self.data

        trackon = self.dmtrack(dm=self.dmarr[dmbin], t0=self.reltime[tbin], show=0)
        if ((trackon[1][0] != 0) | (trackon[1][len(trackon[1])-1] != len(self.chans)-1)):
            print 'Track does not span all channels. Skipping.'
            return [0]

        dataon = data[trackon[0], :, trackon[1]]

        # set up bg track
        if bgwindow:
            # measure max width of pulse (to avoid in bgsub)
            twidths = [] 
            for k in trackon[1]:
                twidths.append(len(n.array(trackon)[0][list(n.where(n.array(trackon[1]) == k)[0])]))

            bgrange = range(-bgwindow/2 - max(twidths) + tbin, -max(twidths) + tbin) + range(max(twidths) + tbin, max(twidths) + bgwindow/2 + tbin)
            for k in bgrange:     # build up super track for background subtraction
                if bgrange.index(k) == 0:   # first time through
                    trackoff = self.dmtrack(dm=self.dmarr[dmbin], t0=self.reltime[k], show=0)
                else:    # then extend arrays by next iterations
                    tmp = self.dmtrack(dm=self.dmarr[dmbin], t0=self.reltime[k], show=0)
                    trackoff[0].extend(tmp[0])
                    trackoff[1].extend(tmp[1])

            dataoff = data[trackoff[0], :, trackoff[1]]

        # compress time axis, then subtract on and off tracks
        for ch in n.unique(trackon[1]):
            indon = n.where(trackon[1] == ch)

            if bgwindow:
                indoff = n.where(trackoff[1] == ch)
                datadiff = dataon[indon].mean(axis=0) - dataoff[indoff].mean(axis=0)
            else:
                datadiff = dataon[indon].mean(axis=0)

            if ch == 0:
                datadiffarr = [datadiff]
            else:
                datadiffarr = n.append(datadiffarr, [datadiff], axis=0)

        datadiffarr = n.array([datadiffarr.transpose()])

        return datadiffarr[0]  # remove time axis, since it is always length 1


    def tracksub2(self, dmbin, tbin, bgwindow = 0):
        """Trying to speed up tracksub...
        """

        data = self.data
        track0,track1 = self.dmtrack0[dmbin]
        trackon = (list(n.array(track0)+tbin), track1)
        twidth = self.twidths[dmbin]
        dataon = data[trackon[0], :, trackon[1]]

#        if ((trackon[1][0] != 0) | (trackon[1][len(trackon[1])-1] != len(self.chans)-1)):
#            print 'Track does not span all channels. Skipping.'
#            return [0]

        # set up bg track
        if bgwindow:
            # measure max width of pulse (to avoid in bgsub)
            bgrange = range(-bgwindow/2 - twidth + tbin, - twidth + tbin) + range(twidth + tbin, twidth + bgwindow/2 + tbin)
            for k in bgrange:     # build up super track for background subtraction
                if bgrange.index(k) == 0:   # first time through
                    trackoff = (list(n.array(track0)+k), track1)
                else:    # then extend arrays by next iterations
                    trackoff = (trackoff[0] + list(n.array(track0)+k), trackoff[1] + track1)

            dataoff = data[trackoff[0], :, trackoff[1]]

        datadiffarr = n.zeros((self.nchan, self.nbl, self.npol),dtype='complex')
        
        # compress time axis, then subtract on and off tracks
        for ch in n.unique(trackon[1]):
            indon = n.where(trackon[1] == ch)

            if bgwindow:
                indoff = n.where(trackoff[1] == ch)
                meanon = dataon[indon].mean(axis=0)
                meanoff = dataoff[indoff].mean(axis=0)
                datadiffarr[ch] = meanon - meanoff
                zeros = n.where( (meanon == 0j) | (meanoff == 0j) )  # find baselines and pols with zeros for meanon or meanoff
                datadiffarr[ch][zeros] = 0j    # set missing data to zero # hack! but could be ok if we can ignore zeros later...
            else:
                datadiffarr[ch] = dataon[indon].mean(axis=0)

        return n.transpose(datadiffarr, axes=[2,1,0])


    def writetrack(self, dmbin, tbin, tshift=0, bgwindow=0, show=0):
        """Writes data from track out as miriad visibility file.
        Optional background subtraction bl-by-bl over bgwindow integrations. Note that this is bgwindow *dmtracks* so width is bgwindow+track width
        Optional spectrum plot with source and background dmtracks
        """

        # prep data and track
        rawdatatrim = self.rawdata[:,:,self.chans]
        track = self.dmtrack(dm=self.dmarr[dmbin], t0=self.reltime[tbin-tshift], show=0)
        if ((track[1][0] != 0) | (track[1][len(track[1])-1] != len(self.chans)-1)):
#            print 'Track does not span all channels. Skipping.'
            return 0

        if bgwindow > 0:
            bgrange = range(-bgwindow/2 - twidths + tbin - tshift, -twidths + tbin - tshift) + range(twidths + tbin - tshift, twidths + bgwindow/2 + tbin - tshift + 1)
#            bgrange = range(int(-bgwindow/2.) + tbin - tshift, int(bgwindow/2.) + tbin - tshift + 1)
#            bgrange.remove(tbin - tshift); bgrange.remove(tbin - tshift + 1); bgrange.remove(tbin - tshift - 1); bgrange.remove(tbin - tshift + 2); bgrange.remove(tbin - tshift - 2); bgrange.remove(tbin - tshift + 3); bgrange.remove(tbin - tshift - 3)
            for i in bgrange:     # build up super track for background subtraction
                if bgrange.index(i) == 0:   # first time through
                    trackbg = self.dmtrack(dm=self.dmarr[dmbin], t0=self.reltime[i], show=0)
                else:    # then extend arrays by next iterations
                    tmp = self.dmtrack(dm=self.dmarr[dmbin], t0=self.reltime[i], show=0)
                    trackbg[0].extend(tmp[0])
                    trackbg[1].extend(tmp[1])
        else:
            print 'Not doing any background subtraction.'

        if show:
            # show source and background tracks on spectrum
            p.figure(1)
            p.plot(self.reltime[track[0]], track[1], 'w.')
            if bgwindow > 0:
                p.plot(self.reltime[trackbg[0]], trackbg[1], 'r.')
            self.spec(save=0)

        # define input metadata source and output visibility file names
        outname = string.join(self.file.split('.')[:-1], '.') + '.' + str(self.nskip/self.nbl) + '-' + 'dm' + str(dmbin) + 't' + str(tbin) + '.mir'
        shutil.rmtree(outname, ignore_errors=True)
        vis = miriad.VisData(self.file,)

        i = 0
        int0 = int(self.nskip + (track[0][len(track[0])/2] + tshift) * self.nbl)   # choose integration at center of dispersed track

        for inp, preamble, data, flags in vis.readLowlevel ('dsl3', False, nocal=True, nopass=True):
            if i == 0:
                # create output vis file
                shutil.rmtree(outname, ignore_errors=True)
                out = miriad.VisData(outname)
                dOut = out.open ('c')

                # set variables
                dOut.setPreambleType ('uvw', 'time', 'baseline')
                dOut.writeVarInt ('nants', self.nants0)
                dOut.writeVarFloat ('inttime', self.inttime0)
                dOut.writeVarInt ('nspect', self.nspect0)
                dOut.writeVarDouble ('sdf', self.sdf0)
                dOut.writeVarInt ('nwide', self.nwide0)
                dOut.writeVarInt ('nschan', self.nschan0)
                dOut.writeVarInt ('ischan', self.ischan0)
                dOut.writeVarDouble ('sfreq', self.sfreq0)
                dOut.writeVarDouble ('restfreq', self.restfreq0)
                dOut.writeVarInt ('pol', self.pol0)
#                inp.copyHeader (dOut, 'history')
                inp.initVarsAsInput (' ') # ???
                inp.copyLineVars (dOut)

            if i < int0:  # need to grab only integration at pulse+intoff
                i = i+1
                continue

            elif i < int0 + self.nbl:
                # write out track, if not flagged

                if n.any(flags):
                    bgarr = []
                    for j in range(self.nchan):
                        if j in self.chans:
                            matches = n.where( j == n.array(self.chans[track[1]]) )[0]   # hack, since chans are from 0-64, but track is in trimmed chan space
                            raw = rawdatatrim[track[0], i-int0, track[1]][matches]   # all baselines for the known pulse
                            raw = raw.mean(axis=0)   # create spectrum for each baseline by averaging over time
                            if bgwindow > 0:   # same as above, but for bg
                                matchesbg = n.where( j == n.array(self.chans[trackbg[1]]) )[0]
                                rawbg = rawdatatrim[trackbg[0], i-int0, trackbg[1]][matchesbg]
                                rawbg = rawbg.mean(axis=0)
                                bgarr.append(rawbg)
                                data[j] = raw - rawbg
                            else:
                                data[j] = raw
                        else:
                            flags[j] = False

#                ants = util.decodeBaseline (preamble[4])
#                print preamble[3], ants

#                dOut.write (self.preamble[i], data, flags)   # not working right here...?
                dOut.write (preamble, data, flags)
                i = i+1  # essentially a baseline*int number

            elif i >= int0 + self.nbl:
                break

        dOut.close ()
        return 1


    def writetrack2(self, dmbin, tbin, tshift=0, bgwindow=0, show=0, pol=0):
        """Writes data from track out as miriad visibility file.
        Alternative to writetrack that uses stored, approximate preamble used from start of pulse, not middle.
        Optional background subtraction bl-by-bl over bgwindow integrations. Note that this is bgwindow *dmtracks* so width is bgwindow+track width
        """

        # create bgsub data
        datadiffarr = self.tracksub(dmbin, tbin, bgwindow=bgwindow)
        if n.shape(datadiffarr) == n.shape([0]):    # if track doesn't cross band, ignore this iteration
            return 0

        data = n.zeros(self.nchan, dtype='complex64')  # default data array. gets overwritten.
        data0 = n.zeros(self.nchan, dtype='complex64')  # zero data array for flagged bls
        flags = n.zeros(self.nchan, dtype='bool')

        # define output visibility file names
        outname = string.join(self.file.split('.')[:-1], '.') + '.' + str(self.nskip/self.nbl) + '-' + 'dm' + str(dmbin) + 't' + str(tbin) + '.mir'
        print outname
        vis = miriad.VisData(self.file,)

        int0 = int((tbin + tshift) * self.nbl)
        flags0 = []
        i = 0
        for inp, preamble, data, flags in vis.readLowlevel ('dsl3', False, nocal=True, nopass=True):
            if i == 0:
                # prep for temp output vis file
                shutil.rmtree(outname, ignore_errors=True)
                out = miriad.VisData(outname)
                dOut = out.open ('c')

                # set variables
                dOut.setPreambleType ('uvw', 'time', 'baseline')
                dOut.writeVarInt ('nants', self.nants0)
                dOut.writeVarFloat ('inttime', self.inttime0)
                dOut.writeVarInt ('nspect', self.nspect0)
                dOut.writeVarDouble ('sdf', self.sdf0)
                dOut.writeVarInt ('nwide', self.nwide0)
                dOut.writeVarInt ('nschan', self.nschan0)
                dOut.writeVarInt ('ischan', self.ischan0)
                dOut.writeVarDouble ('sfreq', self.sfreq0)
                dOut.writeVarDouble ('restfreq', self.restfreq0)
                dOut.writeVarInt ('pol', self.pol0)
#                inp.copyHeader (dOut, 'history')
                inp.initVarsAsInput (' ') # ???
                inp.copyLineVars (dOut)
            if i < self.nbl:
                flags0.append(flags.copy())
                i = i+1
            else:
                break

        l = 0
        for i in range(len(flags0)):  # iterate over baselines
            # write out track, if not flagged
            if n.any(flags0[i]):
                k = 0
                for j in range(self.nchan):
                    if j in self.chans:
                        data[j] = datadiffarr[pol, l, k]
#                        flags[j] = flags0[i][j]
                        k = k+1
                    else:
                        data[j] = 0 + 0j
#                        flags[j] = False
                l = l+1
            else:
                data = data0
#                flags = n.zeros(self.nchan, dtype='bool')

            dOut.write (self.preamble[int0 + i], data, flags0[i])

        dOut.close ()
        return 1


    def makedmt0(self):
        """Integrates data at dmtrack for each pair of elements in dmarr, time.
        Not threaded.  Uses dmthread directly.
        Stores mean of detected signal after dmtrack, effectively forming beam at phase center.
        Probably ok for multipol data...
        """

        dmarr = self.dmarr
#        reltime = n.arange(2*len(self.reltime))/2.  # danger!
        reltime = self.reltime
        chans = self.chans

        dmt0arr = n.zeros((len(dmarr),len(reltime)), dtype='float64')

        for i in range(len(dmarr)):
            for j in range(len(reltime)):
                dmtrack = self.dmtrack(dm=dmarr[i], t0=reltime[j])
                if ((dmtrack[1][0] == 0) & (dmtrack[1][len(dmtrack[1])-1] == len(self.chans)-1)):   # use only tracks that span whole band
#                    dmt0arr[i,j] = n.abs((((self.data).mean(axis=1))[dmtrack[0],dmtrack[1]]).mean())
                    dmt0arr[i,j] = ((((self.data).mean(axis=1))[dmtrack[0],dmtrack[1]]).mean()).real    # use real part to detect on axis, but keep gaussian dis'n
            print 'dedispersed for ', dmarr[i]

        self.dmt0arr = dmt0arr


    def plotdmt0(self, save=0, pathout='./'):
        """calculates rms noise in dmt0 space, then plots circles for each significant point
        save=1 means plot to file.
        """
        dmarr = self.dmarr
        arr = self.dmt0arr
        reltime = self.reltime
        peaks = self.peaks
        tbuffer = 7  # number of extra iterations to trim from edge of dmt0 plot

        # Trim data down to where dmt0 array is nonzero
        arreq0 = n.where(arr == 0)
        trimt = arreq0[1].min()
        arr = arr[:,:trimt - tbuffer]
        reltime = reltime[:trimt - tbuffer]
        print 'dmt0arr/time trimmed to new shape:  ',n.shape(arr), n.shape(reltime)

        mean = arr.mean()
        std = arr.std()
        arr = (arr - mean)/std
        peakmax = n.where(arr == arr.max())

        # Plot
#        p.clf()
        ax = p.imshow(arr, aspect='auto', origin='lower', interpolation='nearest', extent=(min(reltime),max(reltime),min(dmarr),max(dmarr)))
        p.colorbar()

        if len(peaks[0]) > 0:
            print 'Peak of %f at DM=%f, t0=%f' % (arr.max(), dmarr[peakmax[0][0]], reltime[peakmax[1][0]])

            for i in range(len(peaks[1])):
                ax = p.imshow(arr, aspect='auto', origin='lower', interpolation='nearest', extent=(min(reltime),max(reltime),min(dmarr),max(dmarr)))
                p.axis((min(reltime),max(reltime),min(dmarr),max(dmarr)))
                p.plot([reltime[peaks[1][i]]], [dmarr[peaks[0][i]]], 'o', markersize=2*arr[peaks[0][i],peaks[1][i]], markerfacecolor='white', markeredgecolor='blue', alpha=0.5)

        p.xlabel('Time (s)')
        p.ylabel('DM (pc/cm3)')
        p.title('Summed Spectra in DM-t0 space')
        if save:
            savename = self.file.split('.')[:-1]
            savename.append(str(self.nskip/self.nbl) + '.dmt0.png')
            savename = string.join(savename,'.')
            p.savefig(pathout+savename)


    def peakdmt0(self, sig=5.):
        """ Method to find peaks in dedispersed data (in dmt0 space).
        Clips noise, also.
        """
        arr = self.dmt0arr
        reltime = self.reltime
        tbuffer = 7  # number of extra iterations to trim from edge of dmt0 plot

        # Trim data down to where dmt0 array is nonzero
        arreq0 = n.where(arr == 0)
        trimt = arreq0[1].min()
        arr = arr[:,:trimt - tbuffer]
        reltime = reltime[:trimt - tbuffer]
        print 'dmt0arr/time trimmed to new shape:  ',n.shape(arr), n.shape(reltime)

        # single iteration of sigma clip to find mean and std, skipping zeros
        mean = arr.mean()
        std = arr.std()
        print 'initial mean, std:  ', mean, std
        min,max = sigma_clip(arr.flatten())
        cliparr = n.where((arr < max) & (arr > min))
        mean = arr[cliparr].mean()
        std = arr[cliparr].std()
        print 'final mean, sig, std:  ', mean, sig, std

        # Recast arr as significance array
#        arr = n.sqrt((arr-mean)**2 - std**2)/std   # PROBABLY WRONG
        arr = (arr-mean)/std   # for real valued trial output (gaussian dis'n)

        # Detect peaks
        self.peaks = n.where(arr > sig)
        peakmax = n.where(arr == arr.max())
        print 'peaks:  ', self.peaks

        return self.peaks,arr[self.peaks]


    def imagedmt0(self, dmbin, t0bin, tshift=0, bgwindow=4, show=0, clean=1, mode='dirty', code='aipy'):
        """ Makes and fits an background subtracted image for a given dmbin and t0bin.
        tshift can shift the actual t0bin earlier to allow reading small chunks of data relative to pickle.
        code can be 'aipy' or 'miriad'
        """

        # set up
        outroot = string.join(self.file.split('.')[:-1], '.') + '.' + str(self.nskip/self.nbl) + '-dm' + str(dmbin) + 't' + str(t0bin)
        shutil.rmtree (outroot+'.map', ignore_errors=True); shutil.rmtree (outroot+'.beam', ignore_errors=True); shutil.rmtree (outroot+'.clean', ignore_errors=True); shutil.rmtree (outroot+'.restor', ignore_errors=True)

        if code == 'aipy':
            tr = self.tracksub2(dmbin, t0bin, bgwindow=bgwindow)[0].mean(axis=1)
            size = 48000; res = 500   # size=16000 => 13" resolution, res=500 => 7' fov
            fov = n.degrees(1./res)*3600.

            # make image
            ai = aipy.img.Img(size=size, res=res)
            ai.put( (self.u[t0bin],self.v[t0bin],self.w[t0bin]), tr)
            image = ai.image(center = (size/res/2, size/res/2))
            # clean image
            beam = ai.bm_image()
            beamgain = aipy.img.beam_gain(beam[0])
            (clean, dd) = aipy.deconv.clean(image, beam[0], verbose=True, gain=0.01, tol=1e-4)  # light cleaning
            kernel = n.where(beam[0] >= 0.4*beam[0].max(), beam[0], 0.)  # take only peak (gaussian part) pixels of beam image
            restored = aipy.img.convolve2d(clean, kernel)
            image_restored = (restored + dd['res']).real/beamgain
            
            if show:
                p.clf()
                ax = p.imshow(image_restored, aspect='auto', origin='upper', interpolation='nearest', extent=[-fov/2, fov/2, -fov/2, fov/2])
                p.colorbar()
                p.xlabel('Offset (arcsec)')
                p.ylabel('Offset (arcsec)')
                peak = n.where(n.max(image_restored) == image_restored)
            print 'Image peak of %e at (%d,%d)' % (n.max(image_restored), peak[0][0], peak[1][0])
            print 'Peak/RMS = %e' % (image_restored.max()/image_restored.std())
            return image_restored

        elif code == 'miriad':
            if self.approxuvw:
                status = self.writetrack2(dmbin, t0bin, tshift=tshift, bgwindow=bgwindow)   # output file at dmbin, trelbin
            else:
                status = self.writetrack(dmbin, t0bin, tshift=tshift, bgwindow=bgwindow)   # output file at dmbin, trelbin

            if not status:  # if writetrack fails, exit this iteration
                return 0

            try:
        # make image, clean, restor, fit point source
                print
                print 'Making dirty image for nskip=%d, dm[%d]=%.1f, and trel[%d] = %.3f.' % (self.nskip/self.nbl, dmbin, self.dmarr[dmbin], t0bin-tshift, self.reltime[t0bin-tshift])
                txt = TaskInvert (vis=outroot+'.mir', map=outroot+'.map', beam=outroot+'.beam', mfs=True, double=True).snarf()  # good for m31 search
                if show:  txt = TaskCgDisp (in_=outroot+'.map', device='/xs', wedge=True, beambl=True, labtyp='hms,dms', region='relpix,boxes(-100,-100,100,100)').snarf () 
                txt = TaskImStat (in_=outroot+'.map').snarf()   # get dirty image stats

                if mode == 'clean':
#                print 'cleaning image'
                    thresh = 2*float(txt[0][10][41:47])       # set thresh to 2*noise level in dirty image. hack! OMG!!
                    txt = TaskClean (beam=outroot+'.beam', map=outroot+'.map', out=outroot+'.clean', cutoff=thresh, region='relpix,boxes(-100,-100,100,100)').snarf ()   # targeted clean
#                txt = TaskClean (beam=outroot+'.beam', map=outroot+'.map', out=outroot+'.clean', cutoff=thresh).snarf ()
                    print 'Cleaned to %.2f Jy after %d iterations' % (thresh, int(txt[0][-4][19:]))
                    txt = TaskRestore (beam=outroot+'.beam', map=outroot+'.map', model=outroot+'.clean', out=outroot+'.restor').snarf () 
                    if show:  txt = TaskCgDisp (in_=outroot+'.restor', device='/xs', wedge=True, beambl=True, labtyp='hms,dms', region='relpix,boxes(-100,-100,100,100)').snarf () 
                    txt = TaskImFit (in_=outroot+'.restor', object='point').snarf () 

                # parse output of imfit
                # print '012345678901234567890123456789012345678901234567890123456789'
                    peak = float(txt[0][16][30:36])     # 16 -> 14 (etc.) for some images?
                    epeak = float(txt[0][16][44:])
                    off_ra = float(txt[0][17][28:39])
                    eoff_ra = float(txt[0][18][30:39])
                    off_dec = float(txt[0][17][40:])
                    eoff_dec = float(txt[0][18][40:])
                    print 'Fit cleaned image peak %.2f +- %.2f' % (peak, epeak)
                    return peak, epeak, off_ra, eoff_ra, off_dec, eoff_dec
                elif mode == 'dirty':
#                print 'stats of dirty image'
                    peak = float(txt[0][10][50:57])       # get peak of dirty image
                    epeak = float(txt[0][10][41:47])       # note that epeak is biased by any true flux
                    print 'Individual dirty image peak %.2f +- %.2f' % (peak, epeak)
                    if clean:
                        shutil.rmtree (outroot + '.mir', ignore_errors=True)
#                    shutil.rmtree (outroot+'.map', ignore_errors=True) 
                        shutil.rmtree (outroot+'.beam', ignore_errors=True); shutil.rmtree (outroot+'.clean', ignore_errors=True); shutil.rmtree (outroot+'.restor', ignore_errors=True)
                    return peak, epeak
            except:
                print 'Something broke with imaging!'
                return 0
        else:
            print 'code = %s not known' % (code)

    def imsearch(self, dmind, tind, nints, sig=5., show=0, edge=0, mode='dirty'):
        """
        Reproduce search result of pulse_search_image.
        """

        bgwindow = 10  # where bg subtraction is made

        if mode == 'dirty':
            # define typical dirty image noise level for this dm
            print 'For DM = %.1f, measuring median image noise level' % (self.dmarr[dmind])
            bgpeak = []; bgepeak = []
            for bgi in range(bgwindow, nints-bgwindow, nints/15):
                print 'Measuring noise in integration %d' % (bgi)
                outname = string.join(self.file.split('.')[:-1], '.') + '.' + str(self.nskip/self.nbl) + '-' + 'dm' + str(dmind) + 't' + str(bgi) + '.mir'
                shutil.rmtree (outname, ignore_errors=True); shutil.rmtree (outname+'.map', ignore_errors=True); shutil.rmtree (outname+'.beam', ignore_errors=True)
                status = self.writetrack2(dmind, bgi, bgwindow=bgwindow)   # output file at dmbin, trelbin
                try:
                    txt = TaskInvert (vis=outname, map=outname+'.map', beam=outname+'.beam', mfs=True, double=True, cell=80, imsize=250).snarf()
                    txt = TaskImStat (in_=outname+'.map').snarf()   # get dirty image stats
                    bgpeak.append(float(txt[0][10][51:61]))       # get peak of dirty image
                    bgepeak.append(float(txt[0][10][41:51]))       # note that epeak is biased by any true flux
                    shutil.rmtree (outname, ignore_errors=True)
                    shutil.rmtree (outname+'.map', ignore_errors=True)
                    shutil.rmtree (outname+'.beam', ignore_errors=True)
                except:
                    pass
                
            print 'Dirty image noises and their median', bgepeak, n.median(bgepeak)

            # now make dirty image
            results = self.imagedmt0(dmind, tind, show=show, bgwindow=bgwindow, clean=1, mode=mode)
            if mode == 'clean':
                peak, epeak, off_ra, eoff_ra, off_dec, eoff_dec = results
            elif mode == 'dirty':  # need to develop more... needs to be ~10ms processing per int!
                peak, epeak = results
                epeak = n.median(bgepeak)
            if peak/epeak >= sig:
                print '\tDetection!'
            if mode == 'clean':
                print self.nskip/self.nbl, nints, (dmind, tind), 'Peak, (sig),  RA, Dec: ', peak, epeak, '(', peak/epeak, ')  ', off_ra, eoff_ra, off_dec, eoff_dec
            elif mode == 'dirty':
                print self.nskip/self.nbl, nints, (dmind, tind), 'Peak, (sig): ', peak, epeak, '(', peak/epeak, ')'


    def uvfitdmt0(self, dmbin, t0bin, bgwindow=4, tshift=0, show=1, mode='fit'):
        """ Makes and fits a point source to background subtracted visibilities for a given dmbin and t0bin.
        tshift can shift the actual t0bin earlier to allow reading small chunks of data relative to pickle.
        """
        
        # set up
        outroot = string.join(self.file.split('.')[:-1], '.') + '.' + str(self.nskip/self.nbl) + '-dm' + str(dmbin) + 't' + str(t0bin)
        shutil.rmtree (outroot + '.mir', ignore_errors=True)

        if self.approxuvw:
            status = self.writetrack2(dmbin, t0bin, tshift=tshift, bgwindow=bgwindow)   # output file at dmbin, trelbin
        else:
            status = self.writetrack(dmbin, t0bin, tshift=tshift, bgwindow=bgwindow)   # output file at dmbin, trelbin

        if not status:  # if writetrack fails, exit this iteration
            return 0

        if mode == 'fit':
            try:
                # fit point source model to visibilities
                print
                print 'UVfit for nskip=%d, dm[%d] = %.1f, and trel[%d] = %.3f.' % (self.nskip/self.nbl, dmbin, self.dmarr[dmbin], t0bin-tshift, self.reltime[t0bin-tshift])
                txt = TaskUVFit (vis=outroot+'.mir', object='point', select='-auto').snarf()

                # parse output of imfit
                # print '012345678901234567890123456789012345678901234567890123456789'
                peak = float(txt[0][8][30:38])
                epeak = float(txt[0][8][46:])
                off_ra = float(txt[0][9][30:38])
                eoff_ra = float(txt[0][10][31:42])
                off_dec = float(txt[0][9][40:])
                eoff_dec = float(txt[0][10][42:])
                print 'Fit peak %.2f +- %.2f' % (peak, epeak)
                shutil.rmtree (outroot + '.mir', ignore_errors=True)
                return peak, epeak, off_ra, eoff_ra, off_dec, eoff_dec
            except:
                print 'Something broke in/after uvfit!'
                shutil.rmtree (outroot + '.mir', ignore_errors=True)
                return 0

        elif mode == 'grid':
            llen = 150
            mlen = 150
            linspl = n.linspace(-4000,4000,llen)
            linspm = n.linspace(-4000,4000,mlen)
            snrarr = n.zeros((llen, mlen))
            
            for i in range(llen):
                print 'Starting loop ', i
                for j in range(mlen):
                    try:
                        txt = TaskUVFit (vis=outroot+'.mir', object='point', select='-auto', spar='0,'+str(linspl[i])+','+str(linspm[j]), fix='xy').snarf()
#                        print txt[0][8:10]
                        peak = float(txt[0][8][30:38])
                        epeak = float(txt[0][8][46:])
                        off_ra = float(txt[0][9][30:38])
                        off_dec = float(txt[0][9][40:])
                        snrarr[i, j] = peak/epeak
                    except:
                        print 'Something broke in/after uvfit!'

            peak = n.where(snrarr == snrarr.max())
            print 'Peak: ', snrarr[peak]
            print 'Location: ', peak, (linspl[peak[0]], linspm[peak[1]])
            print
            print 'Center SNR: ', snrarr[llen/2,llen/2]

            log = open('log.txt','a')
            print >> log, outroot, snrarr[peak], peak, (linspl[peak[0]], linspm[peak[1]])
            log.close()

            ax = p.imshow(snrarr, aspect='auto', origin='lower', interpolation='nearest')
            p.colorbar()
            p.savefig(outroot + '.png')


    def uvfitdmt02 (self, dmbin, t0bin, bgwindow=4, tshift=0, show=1, mode='default'):
        """Experimental alternative function to do uvfitting of visibilities.
        Reads in data from dispersion track, then fits uv model of point source in python.
        mode defines the optimization approach.
        """

        datadiffarr = self.tracksub(dmbin, t0bin, bgwindow=bgwindow)
        int0 = int((t0bin + tshift) * self.nbl)
        meanfreq = n.mean(self.sfreq + self.sdf * self.chans )

        datadiffarr2 = self.tracksub(dmbin, t0bin+bgwindow-20, bgwindow=bgwindow)   # arbitrary offset to measure typical noise in bg
        datadiffarr3 = self.tracksub(dmbin, t0bin+bgwindow-15, bgwindow=bgwindow)   # arbitrary offset to measure typical noise in bg
        datadiffarr4 = self.tracksub(dmbin, t0bin+bgwindow+15, bgwindow=bgwindow)   # arbitrary offset to measure typical noise in bg
        datadiffarr5 = self.tracksub(dmbin, t0bin+bgwindow+20, bgwindow=bgwindow)   # arbitrary offset to measure typical noise in bg
        datadiffarr6 = self.tracksub(dmbin, t0bin+bgwindow+25, bgwindow=bgwindow)   # arbitrary offset to measure typical noise in bg
        obsrms2 = n.std((datadiffarr2.mean(axis=1)).imag)      # std of imag part is std of real part
        obsrms3 = n.std((datadiffarr3.mean(axis=1)).imag)      # std of imag part is std of real part
        obsrms4 = n.std((datadiffarr4.mean(axis=1)).imag)      # std of imag part is std of real part
        obsrms5 = n.std((datadiffarr5.mean(axis=1)).imag)      # std of imag part is std of real part
        obsrms6 = n.std((datadiffarr6.mean(axis=1)).imag)      # std of imag part is std of real part
        obsrms = n.median([obsrms2, obsrms3, obsrms4, obsrms5, obsrms6])
        print 'obsrms:', obsrms
        p.plot(datadiffarr.mean(axis=0).mean(axis=1).real, datadiffarr.mean(axis=1).imag, '.')  # should be ok for multipol data...
        p.show()

        # get flags to help select good ants
        vis = miriad.VisData(self.file,)
        flags0 = []
        i = 0
        for inp, preamble, data, flags in vis.readLowlevel ('dsl3', False, nocal=True, nopass=True):
            if i < self.nbl:
                flags0.append(flags.copy())
                i = i+1
            else:
                break
        goodants = n.where( n.any(n.array(flags0), axis=1) == True)

        u = []; v = []; w = []
        for i in range(self.nbl):  
            u.append(self.preamble[int0 + i][0])  # in units of ns
            v.append(self.preamble[int0 + i][1])  # in units of ns
            w.append(self.preamble[int0 + i][2])  # in units of ns
        u = n.array(u); v = n.array(v); w = n.array(w)
        u = meanfreq * u[goodants]; v = meanfreq * v[goodants]; w = meanfreq * w[goodants]
        data = datadiffarr.mean(axis=1)

        print
        print 'UVfit2 for nskip=%d, dm[%d] = %.1f, and trel[%d] = %.3f.' % (self.nskip/self.nbl, dmbin, self.dmarr[dmbin], t0bin-tshift, self.reltime[t0bin-tshift])

        vi = lambda a, l, m, u, v, w: a * n.exp(-2j * n.pi * (u*l + v*m)) # + w*n.sqrt(1 - l**2 - m**2)))  # ignoring w term
        def errfunc (p, u, v, w, y, obsrms):
            a, l, m = p
            err = (y - vi(a, l, m, u, v, w))/obsrms
            return err.real + err.imag

        def errfunc2 (p, l, m, u, v, w, y):
            a = p
            err = (y - vi(a, l, m, u, v, w))
            return err.real + err.imag

        def ltoa (l):
            a = n.arcsin(l) * 180 / n.pi * 3600
            return a

        if mode == 'default':
            p0 = [100.,0.,0.]
            out = opt.leastsq(errfunc, p0, args = (u, v, w, data, obsrms), full_output=1)
            p1 = out[0]
            covar = out[1]
            print 'Fit results (Jy, arcsec, arcsec):', p1[0], ltoa(p1[1]), ltoa(p1[2])
            print 'Change in sumsq: %.2f to %.2f' % (n.sum(errfunc(p0, u, v, w, data, obsrms)**2), n.sum(errfunc(p1, u, v, w, data, obsrms)**2))

            print 'here\'s some stuff...'
            print covar
            print n.sqrt(covar[0][0])
            print n.sqrt(covar[1][1])
        elif mode == 'map':
            p0 = [100.,0.,0.]
            out = opt.leastsq(errfunc, p0, args = (u, v, w, data, obsrms), full_output=1)
            p1 = out[0]
            covar = out[1]

            llen = 20
            mlen = 20
            sumsq = n.zeros((llen, mlen))
            linspl = n.linspace(-0.005,0.005,llen)   # dl=0.026 => 1.5 deg (half ra width of andromeda)
            linspm = n.linspace(-0.005,0.005,mlen)   # dl=0.01 => 0.57 deg (half dec width of andromeda)
            for i in range(llen):
                for j in range(mlen):
                    sumsq[i, j] = n.sum(errfunc([p1[0], linspl[i], linspm[j]], u, v, w, data, obsrms)**2)

            mins = n.where(sumsq < 1.01 * sumsq.min())
            print 'Best fit: ', p1
            print 'Red. chisq: ', sumsq[mins][0]/(9-3)
            print 'Location: ', mins, (ltoa(linspl[mins[0]])[0], ltoa(linspm[mins[1]])[0])
            p.plot(mins[1], mins[0], 'w.')
            p.imshow(sumsq)
            p.colorbar()
            p.show()
        elif mode == 'lsgrid':
            llen = 30
            mlen = 30
            linspl = n.linspace(-0.005,0.005,llen)   # dl=0.026 => 1.5 deg (half ra width of andromeda)
            linspm = n.linspace(-0.005,0.005,mlen)   # dl=0.01 => 0.57 deg (half dec width of andromeda)
            amparr = n.zeros((llen, mlen))
            p0 = [100.]
            
            for i in range(llen):
                for j in range(mlen):
                    out = opt.leastsq(errfunc2, p0, args = (linspl[i], linspm[j], u, v, w, data), full_output=1)
                    amparr[i, j] = out[0]

            maxs = n.where(amparr >= 0.9*amparr.max())
            print 'Peak: ', amparr.max()
            print 'Location: ', maxs, (ltoa(linspl[maxs[0]]), ltoa(linspm[maxs[1]]))
            print 'SNR: ', amparr.max()/(obsrms/n.sqrt(9))
            print
            print 'Center: ', amparr[15,15]
            print 'Center SNR: ', amparr[15,15]/(obsrms/n.sqrt(9))

            ax = p.imshow(amparr, aspect='auto', origin='lower', interpolation='nearest')
            p.colorbar()
            p.show()


    def tripgen(self, amin=0, amax=0):
        """Calculates and returns data indexes (i,j,k) for all closed triples.
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

    def bisplc(self, chan=-1, dmbin=0, bgwindow=4, show=0, save=0, sigma=5., pathout='./'):
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

        twidth = n.round(self.twidths[dmbin])
        dmwidth = int(n.round(n.max(self.dmtrack0[dmbin][0]) - n.min(self.dmtrack0[dmbin][0])))

        bisp = lambda d,i,j,k: d[:,i] * d[:,j] * n.conj(d[:,k])    # bispectrum for pol data
#        bisp = lambda d,i,j,k: n.complex(d[i] * d[j] * n.conj(d[k]))

        triples = self.tripgen()

        # set up arrays for filling or weighting data
        dibi = n.zeros((len(self.data)-(bgwindow+2*twidth+dmwidth), len(triples)), dtype='complex')
        truearr = n.ones( (self.npol, self.nbl, self.nchan))
        falsearr = n.zeros( (self.npol, self.nbl, self.nchan))

        for ii in range(len(self.data)-(bgwindow+2*twidth+dmwidth)):
            diff = self.tracksub2(dmbin, ii+bgwindow/2+twidth, bgwindow=bgwindow)
            if len(n.shape(diff)) == 1:    # no track
                continue
            weightarr = n.where(diff != 0j, truearr, falsearr)  # ignore zeros in mean across channels # bit of a hack
            try:
                diffmean = n.average(diff, axis=2, weights=weightarr)
            except ZeroDivisionError:
                diffmean = n.mean(diff, axis=2)    # if all zeros, just make mean # bit of a hack

            for trip in range(len(triples)):
                i, j, k = triples[trip]
                dibi[ii, trip] = bisp(diffmean, i, j, k).mean(axis=0)  # Stokes I bispectrum. should be ok for multipol data...

        if show or save:
            good = self.det_bisplc(dibi, save=save, show=show, sigma=sigma, pathout=pathout, dmbin=dmbin)
            return good
        else:
            return dibi


    def det_bisplc(self, dibi, sigma=5., tol=1.3, show=0, save=0, pathout='./', dmbin=0):
        """Function to detect source in bisplc.
        sigma gives the threshold for SNR_bisp (apparent). 
        tol gives the amount of tolerance in the sigma_b cut for point-like sources (rfi filter).
        """

#        ntr = lambda num: num*(num-1)*(num-2)/6.   # assuming all triples are present
        ntr = lambda num: n.mean([len(n.where(dibi[i] != 0j)[0]) for i in range(len(dibi))])   # consider possibility of zeros in data and take mean number of good triples over all times

        # using s=S/Q
#        mu = lambda s: s/(1+s)  # for independent, calibrated bispectra (0 -> 1)
#        mu = lambda s: (s+1/3.)/(s+1)       # for dependent, calibrated bispectra (1/3 -> 1)
        mu = lambda s: 1.  # for bispectra at high S/N from visibilities?
        sigbQ3 = lambda s: n.sqrt((1 + 3*mu(s)**2) + 3*(1 + mu(s)**2)*s**2 + 3*s**4)  # from kulkarni 1989, normalized by Q**3, also rogers et al 1995
        s = lambda dibisnr, num: (2.*dibisnr/n.sqrt(ntr(num)))**(1/3.)

        # measure SNR_bl==Q from sigma clipped times with normal mean and std of bispectra
        dibimean = dibi.real.mean(axis=1)
        dibistd = dibi.real.std(axis=1)
        (meanmin,meanmax) = sigma_clip(dibimean)  # remove rfi
        (stdmin,stdmax) = sigma_clip(dibistd)  # remove rfi
        clipped = n.where((dibimean > meanmin) & (dibimean < meanmax) & (dibistd > stdmin) & (dibistd < stdmax) & (dibimean != 0.0))[0]  # remove rfi
        dibimeanstd = dibi[clipped].real.mean(axis=1).std()
#        dibimeanstd = dibi[clipped].std()  # this catches stuff in b0329 data, but is not quite right
        dibisnr = dibimean/dibimeanstd    # = S**3/(Q**3 / n.sqrt(n_tr)) = s**3 * n.sqrt(n_tr)
        Q = ((dibimeanstd/2.)*n.sqrt(ntr(self.nants)))**(1/3.)
#        Q = (dibimeanstd*n.sqrt(ntr(self.nants)/(self.nants-2)))**(1/3.)    # includes mu^(i) = 1/3 for baseline-based bispectra?
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
            p.title(str(self.scan) + ' scan, ' + str(self.nskip/self.nbl) + ' nskip, ' + str(dmbin) + ' dmbin, ' + str(len(cands))+' candidates', transform = ax.transAxes)
            p.plot(dibisnr, 'b.')
            if len(cands) > 0:
                p.plot(cands, dibisnr[cands], 'r*')
                p.ylim(-2*dibisnr[cands].max(),2*dibisnr[cands].max())
            p.xlabel('Integration')
            p.ylabel('SNR_b')
            p.subplot(212)
            p.plot(dibistd/Q**3, dibisnr, 'b.')

            # plot reference theory lines
            smax = s(dibisnr.max(), self.nants)
            sarr = smax*n.arange(0,51)/50.
            p.plot(sigbQ3(sarr), 1/2.*sarr**3*n.sqrt(ntr(self.nants)), 'k')
            p.plot(tol*sigbQ3(sarr), 1/2.*sarr**3*n.sqrt(ntr(self.nants)), 'k--')
            p.plot(dibistd[cands]/Q**3, dibisnr[cands], 'r*')

            if len(cands) > 0:
                p.axis([0, tol*sigbQ3(s(dibisnr[cands].max(), self.nants)), -0.5*dibisnr[cands].max(), 1.1*dibisnr[cands].max()])

                # show spectral modulation next to each point
                for candint in cands:
                    sm = n.single(round(self.specmod(dmbin,candint),1))
                    p.text(dibistd[candint]/Q**3, dibisnr[candint], str(sm), horizontalalignment='right', verticalalignment='bottom')
            p.xlabel('sigma_b/Q^3')
            p.ylabel('SNR_b')
            if save:
                savename = self.file.split('.')[:-1]
                savename.append(str(self.scan) + '_' + str(self.nskip/self.nbl) + '_' + str(dmbin) + '.bisplc.png')
                savename = string.join(savename,'.')
                p.savefig(pathout+savename)
            else:
                pass

        return dibisnr[cands], dibistd[cands], cands


    def specmod(self, dmbin, tbin, bgwindow=4):
        """Calculate spectral modulation for given dmtrack.
        Narrow RFI has large (>5) modulation, while spectrally broad emission has low modulation.
        See Spitler et al 2012 for details.
        """

#        smarr = n.zeros(len(self.dataph))  # uncomment to do specmod lightcurve
#        for int in range(len(self.dataph)-bgwindow):
        diff = self.tracksub2(dmbin, tbin, bgwindow=bgwindow)
        bfspec = diff.mean(axis=0).real  # should be ok for multipol data...
        sm = n.sqrt( ((bfspec**2).mean() - bfspec.mean()**2) / bfspec.mean()**2 )

        return sm


    def phaseshift(self, l, m):
        """Function to apply phase shift to (l,m) coordinates of data array.
        Should return new data array.
        This should be used before .mean(axis=bl) step is done to produce a new spectrogram.
        """

        newdata = n.zeros(shape=self.data.shape, dtype='complex')
        ang = lambda l,m,u,v,freq: l*n.outer(u,freq/freq.mean()) + m*n.outer(v,freq/freq.mean())  # operates on single time of u,v

        print 'Shifting phase center by (l,m) = (%e,%e) = (%e,%e) arcsec' % (l, m, n.degrees(l)*3600, n.degrees(m)*3600)

        for int in range(len(newdata)):
            for pol in range(self.npol):
                newdata[int,:,:,pol] = self.data[int,:,:,pol] * n.exp(-2j*n.pi*ang(l, m, self.u[int], self.v[int], self.freq))
    
        return newdata


    def dmlc(self, dmbin, tbin, nints = 50):
        """Plots lc for DM bin over range of timebins.
        In principle, should produce lightcurve as if it is a slice across dmt0 plot.
        Actually designed to test writetrack function.
        """
        pass


    def normalreim(self, prob=2.3e-4, bgwindow=0):
        """Calculates p-value of normality for real-imaginary distribution of visibilities (uses real/imag separately).
        Uses baselines and channels separately. prob is the false positive rate (actually non-normal p-value); default is 230/1e6 => 5sigma.
        Returns least likely normal trials
        Probably not ok for multipol data...
        """

        write = 0  # use writetrack to get bgsub visies? takes long time...
        tbuffer = 7  # number of extra iterations to trim from edge of dmt0 plot

        dmarr = self.dmarr
        reltime = self.reltime
        chans = self.chans

        dmt0arr = n.zeros((len(dmarr),len(reltime)), dtype='float64')

        for i in range(len(dmarr)):
            for j in range(len(reltime)):
                if write:
                    # use writetrack to get bgsub visibilities
                    status = self.writetrack(i, j, tshift=0, bgwindow=bgwindow)
                    if status:
                        newfile = string.join(self.file.split('.')[:-1], '.') + '.' + str(self.nskip/self.nbl) + '-' + 'dm' + str(i) + 't' + str(j) + '.mir'
                        print 'Loading file', newfile
                        ev2 = evla(newfile, nints=1)
                        length = ev2.data.shape[1] * ev2.data.shape[2]
                        da = (ev2.data[0]).reshape((1,length))[0]
                        shutil.rmtree(newfile, ignore_errors=True)
                        dmt0arr[i,j] = min(morestats.shapiro(da.real)[1], morestats.shapiro(da.imag)[1])
                else:
                    datadiff = self.tracksub(i, j, bgwindow=bgwindow)
                    if len(n.shape(datadiff)) == 3:
                        length = datadiff.shape[1] * datadiff.shape[2]
                        datadiff = datadiff.reshape((1,length))[0]
                        dmt0arr[i,j] = min(morestats.shapiro(datadiff.real)[1], morestats.shapiro(datadiff.imag)[1])
                    else:
                        continue
            print 'dedispersed for ', dmarr[i]

        # Trim data down to where dmt0 array is nonzero
        arreq0 = n.where(dmt0arr == 0)
        trimt = arreq0[1].min()
        dmt0arr = dmt0arr[:,:trimt - tbuffer]
        reltime = reltime[:trimt - tbuffer]
        print 'dmt0arr/time trimmed to new shape:  ',n.shape(dmt0arr), n.shape(reltime)

        # Detect dips
        self.dmt0arr = dmt0arr
        self.dips = n.where(dmt0arr < prob)
        dipmin = n.where(dmt0arr == dmt0arr.min())
        print 'Least normal re-im distributions: ', self.dips
        print 'Number of trials: ', len(dmarr) * len(reltime)

        return self.dips,dmt0arr[self.dips]


    def stdreim(self, fstd=1.2, bgwindow=0, show=0):
        """Calculates standard deviation of real-imaginary distribution of visibilities. Uses baselines and channels separately. 
        fstd is threshold factor of change in std to make detection.
        Returns trials with large std change.
        Probably not ok for multipol data...
        """

        write = 0  # use writetrack to get bgsub visies? takes long time...
        tbuffer = 7  # number of extra iterations to trim from edge of dmt0 plot

        dmarr = self.dmarr
        reltime = self.reltime
        chans = self.chans

        dmt0arr = n.zeros((len(dmarr),len(reltime)), dtype='float64')

        for i in range(len(dmarr)):
            for j in range(len(reltime)):
                datadiff = self.tracksub(i, j, bgwindow=bgwindow)
                if len(n.shape(datadiff)) == 3:
                    dmt0arr[i,j] = n.std(datadiff)
                else:
                    continue
            print 'dedispersed for ', dmarr[i]

        # Trim data down to where dmt0 array is nonzero
        arreq0 = n.where(dmt0arr == 0)
        trimt = arreq0[1].min()
        dmt0arr = dmt0arr[:,:trimt - tbuffer]
        reltime = reltime[:trimt - tbuffer]
        print 'dmt0arr/time trimmed to new shape:  ',n.shape(dmt0arr), n.shape(reltime)

        # Detect peaks
        self.dmt0arr = dmt0arr
        self.peaks = n.where(dmt0arr > fstd * dmt0arr.mean())
        peakmax = n.where(dmt0arr == dmt0arr.max())
        print 'Least normal re-im distributions: ', self.peaks
        print 'Number of trials: ', len(dmarr) * len(reltime)

        if show:
            for i in range(len(self.peaks[1])):
                ax = p.imshow(dmt0arr, aspect='auto', origin='lower', interpolation='nearest', extent=(min(reltime),max(reltime),min(dmarr),max(dmarr)))
                p.axis((min(reltime),max(reltime),min(dmarr),max(dmarr)))
                p.plot([reltime[self.peaks[1][i]]], [dmarr[self.peaks[0][i]]], 'o', markersize=2*dmt0arr[self.peaks[0][i],self.peaks[1][i]], markerfacecolor='white', markeredgecolor='blue', alpha=0.5)

            p.xlabel('Time (s)')
            p.ylabel('DM (pc/cm3)')
            p.title('Summed Spectra in DM-t0 space')

        return self.peaks,dmt0arr[self.peaks]


    def plotreim(self, save=0, pathout='./'):
        """Plots the visibilities in real-imaginary space. Test of pulse detection concept for uncalibrated data...
        Probably not ok for multipol data...
        """

        da = self.data.mean(axis=0)
        length = da.shape[0] * da.shape[1]
        da = da.reshape((1,length))[0]

        print 'Real center: %.1f +- %.1f ' % (da.real.mean(), da.real.std()/n.sqrt(len(da.real)))
        print 'Imag center: %.1f +- %.1f ' % (da.imag.mean(), da.imag.std()/n.sqrt(len(da.real)))
        print 'Normal p-value (real): ', morestats.shapiro(da.real)[1]
        print 'Normal p-value (imag): ', morestats.shapiro(da.imag)[1]

        if save:
            savename = self.file.split('.')[:-1]
            savename.append(str(self.nskip/self.nbl) + '.reim.png')
            savename = string.join(savename,'.')
            p.savefig(pathout+savename)
        else:
            p.plot(da.real,da.imag, '.')
            p.show()


    def dedisperse2(self):
        """Integrates over data at dmtrack for each pair of elements in dmarr, time.
        Uses threading.  SLOWER than serial.
        """

        dmarr = self.dmarr
        reltime = self.reltime
        chans = self.chans

        self.dmt0arr = n.zeros((len(dmarr),len(reltime)), dtype='float64')
#        accummask = n.zeros(self.dataph.shape, dtype='bool')

        threadlist = []
        for i in range(len(dmarr)):
            for j in range(len(reltime)):
                proc = worker(self, i, j)
                threadlist.append(proc)
                proc.start()
            print 'submitted for dm= ', dmarr[i]
        for proc in threadlist:
            proc.join()


    def dedisperse3(self):
        """Integrates over data at dmtrack for each pair of elements in dmarr, time.
        Uses ipython.
        """
        from IPython.kernel import client

        dmarr = self.dmarr
        reltime = self.reltime
        dmt0arr = n.zeros((len(dmarr),len(reltime)), dtype='float64')

        # initialize engines
        tc = client.TaskClient()
        mec = client.MultiEngineClient()
        ids = mec.get_ids()
        print 'Got engines: ', ids
#        mec.push_function(dict(dmtrack2 = dmtrack2))    # set up function for later
#        mec.push(dict(data=self.dataph, reltime=reltime))
#        mec.execute('import numpy as n')
#        mec.push_function(dict(dmtrack2 = dmtrack2))    # set up function for later

        pr_list = []
        iarr = []; jarr = []
        for i in range(len(dmarr)):
            for j in range(len(reltime)):
#                iarr.append(i)
#                jarr.append(j)
#                k = (j + len(reltime) * i) % len(ids)   # queue each iter amongst targets
                st = client.StringTask('dmt0 = dmtrack2(data, reltime, %f, %f)' % (dmarr[i], reltime[j]), pull='dmt0', push=dict(data = self.dataph,reltime = reltime))
                pr_list.append(tc.run(task=st))
                if len(pr_list) == len(ids):     # wait until node queue is full
                    for l in range(len(pr_list)):
                        tc.barrier(pr_list)         # then wait for all processes to finish
                        dmt0 = tc.get_task_result(pr_list[l])
                        print dmt0
#                        dmt0arr[iarr[l],jarr[l]] = dmt0[l]
#            iarr = []; jarr = []
            print 'dedispersed for ', dmarr[i]

        self.dmt0arr = dmt0arr


def pulse_search_phasecenter(fileroot, pathin, pathout, nints=10000, edge=0):
    """Blind search for pulses at phase center.
    TO DO:  improve handling of edge times and ignoring data without complete DM track?
    """

    maxints = 90000
    filelist = []
#    filelist.append(string.join(fileroot.split('.')[:-1]) + '_rr1.mir')
    filelist.append(string.join(fileroot.split('.')[:-1]) + '_rr2.mir')

    # loop over miriad data and time chunks
    for file in filelist:
        fileout = open(pathout + string.join(file.split('.')[:-1], '.') + '.txt', 'a')

        for nskip in range(0, maxints-(nints-edge), nints-edge):
            print
            print 'Starting file %s with nskip %d' % (file, nskip)

            # load data
            ev = evla(pathin + file, nints=nints, nskip=nskip)

            # searches at phase center  ## TO DO:  need to search over position in primary beam
            ev.makedmt0()
            peaks, peakssig = ev.peakdmt0()
            print >> fileout, file, nskip, nints, peaks

            # save all results (v1.0 pickle format)
            # TO DO:  think of v2.0 of pickle format
            if len(peaks[0]) > 0:
                pklout = open(pathout + string.join(file.split('.')[:-1], '.') + '.' + str(nskip) + '.pkl', 'wb')
                pickle.dump((file, nskip, nints, peaks[0], ev.dmarr[peaks[0][0]], peaks[1], peakssig), pklout)
                pklout.close()

        fileout.close


def pulse_search_bisp(fileroot, pathin, pathout, nints=1000, edge=30):
    """Blind search for pulses with bispectrum algorithm.
    """

    startint = 0
    maxints = 24000
    scans = 4
    bgwindow = 4
    filelist = [fileroot]

    for file in filelist:
        for scan in range(scans):
            for nskip in range(startint, maxints-(nints-edge), nints-edge):
                print
                print 'Starting file %s in scan %s with nskip %d' % (file, scan, nskip)
                fileout = open(pathout + string.join(file.split('.')[:-1], '.') + '.txt', 'a')

                ev = evla(pathin + file, nints=nints, nskip=nskip, selectpol=['RR','LL'], ddid=[0,1], datacol='data', scan=scan)
                if ev.min == ev.max: break   # escape if no data
                results = {}
                for dmbin in range(len(ev.dmarr)):
                    dibi = ev.bisplc(bgwindow=bgwindow, dmbin=dmbin)
                    dets = ev.det_bisplc(dibi, sigma=5, save=1, pathout=pathout, dmbin=dmbin)
                    results[dmbin] = dets

                snragg = [results[dm][0][i] for dm in results.keys() for i in range(len(results[dm][0]))]
                indagg = [results[dm][2][i] for dm in results.keys() for i in range(len(results[dm][2]))]
                dmbinagg = []; dmagg = []
                for dmbin in results.keys():
                    dmbinagg = dmbinagg + [dmbin]*len(results[dmbin][0])
                    dmagg = dmagg + [ev.dmarr[dmbin]]*len(results[dmbin][0])

                print 'Candidates aggregated over DM:'
                print results

                if len(snragg) > 0:
                    print >> fileout, file, scan, nskip, nints, dmagg, indagg, snragg
                    ev.spec(n.unique(indagg), save=1, pathout=pathout)
                    
                    # TO DO:  think of v2.0 of pickle format
                    pklout = open(pathout + string.join(file.split('.')[:-1], '.') + '.' + str(nskip) + '.' + str(scans) + '.pkl', 'wb')
                    pickle.dump((file, scan, nskip, nints, n.array(dmbinagg), n.array(dmagg), n.array(indagg), n.array(snragg)), pklout)   # pkl v1.0
                    pklout.close()
                else:
                    print >> fileout, file, scan, nskip, nints
                    ev.spec(save=1, pathout=pathout)
                fileout.close


def pulse_search_image(fileroot, pathin, pathout, nints=12000, sig=5.0, show=0, edge=0, mode='dirty', dmrange=None, nstart=0, tstop=0):
    """
    Searches for pulses by imaging dedispersed trials.
    dmrange lets outside call limit range of dms to search (good to spread jobs out in parallel).
    nstart and tstop control start integration and duration of run. useful for running on cluster.
    """

    if tstop != 0:
        import time
        t0 = time.time()
        print 'Time limited search starting at ', time.ctime()

    maxints = 131000  # biggest file in integrations
    bgwindow = 10  # where bg subtraction is made

    filelist = []

# for crab 201103
#    for i in [0,1,4,5,6,7,8,9]:
#        filelist.append(string.join(fileroot.split('.')[:-1]) + '_0' + str(i) + '.mir')
#    for i in range(0,8):
#        filelist.append(string.join(fileroot.split('.')[:-1]) + '_1' + str(i) + '.mir')

# for crab 190348
#    for i in [0,1,2,3,4,5,6,7,8,9]:
#        filelist.append(string.join(fileroot.split('.')[:-1]) + '_0' + str(i) + '.mir')
#    for i in [1,3,4,5,6,7,8,9]:
#        filelist.append(string.join(fileroot.split('.')[:-1]) + '_1' + str(i) + '.mir')
#    for i in [0,1,2,3]:
#        filelist.append(string.join(fileroot.split('.')[:-1]) + '_2' + str(i) + '.mir')

# for b0329 173027
#    for i in [0,1,2,3,4,5,6,7,8,9]:
#        filelist.append(string.join(fileroot.split('.')[:-1]) + '_0' + str(i) + '.mir')
#    for i in [0,1,2,6,7,8,9]:
#        filelist.append(string.join(fileroot.split('.')[:-1]) + '_1' + str(i) + '.mir')
#    for i in [0,1,2,3]:
#        filelist.append(string.join(fileroot.split('.')[:-1]) + '_2' + str(i) + '.mir')

# for m31 154202
#    for i in [0,1,2,3,4,5,6,7,8,9]:
#        filelist.append(string.join(fileroot.split('.')[:-1]) + '_0' + str(i) + '.mir')
#    for i in [0,1,2,3,4,5,6]:
#        filelist.append(string.join(fileroot.split('.')[:-1]) + '_1' + str(i) + '.mir')

# hack to search single file for pulses
    filelist = [fileroot]

    print 'Looping over filelist ', filelist, ' with dmrange, ', dmrange
    for file in filelist:
        for nskip in range(nstart, maxints-(nints-edge), nints-edge):
            print 'Starting file %s with nskip %d' % (file, nskip)

            # load data
            ev = evla(pathin + file, nints=nints, nskip=nskip)

            if dmrange == None:
                dmrange = range(len(ev.dmarr))

            fileout = open(pathout + string.join(file.split('.')[:-1], '.') + '_dm' + str(dmrange[0]) + '.txt', 'a')

            # dedisperse
            for i in dmrange:

                if mode == 'dirty':
                # define typical dirty image noise level for this dm
                    print 'For DM = %.1f, measuring median image noise level' % (ev.dmarr[i])
                    bgpeak = []; bgepeak = []
                    for bgi in range(bgwindow, nints-bgwindow, nints/15):
                        print 'Measuring noise in integration %d' % (bgi)
                        outname = string.join(ev.file.split('.')[:-1], '.') + '.' + str(ev.nskip/ev.nbl) + '-' + 'dm' + str(i) + 't' + str(bgi) + '.mir'
                        shutil.rmtree (outname, ignore_errors=True); shutil.rmtree (outname+'.map', ignore_errors=True); shutil.rmtree (outname+'.beam', ignore_errors=True)
                        status = ev.writetrack2(i, bgi, bgwindow=bgwindow)   # output file at dmbin, trelbin
                        try:
                            txt = TaskInvert (vis=outname, map=outname+'.map', beam=outname+'.beam', mfs=True, double=True, cell=80, imsize=250).snarf()
                            txt = TaskImStat (in_=outname+'.map').snarf()   # get dirty image stats
                            bgpeak.append(float(txt[0][10][51:61]))       # get peak of dirty image
                            bgepeak.append(float(txt[0][10][41:51]))       # note that epeak is biased by any true flux
                            shutil.rmtree (outname, ignore_errors=True)
                            shutil.rmtree (outname+'.map', ignore_errors=True)
                            shutil.rmtree (outname+'.beam', ignore_errors=True)
                        except:
                            pass
                    print 'Dirty image noises and their median', bgepeak, n.median(bgepeak)
                # now iterate over integrations
                for j in range(len(ev.reltime)):
                    if (tstop != 0) & ((j % 10) == 0):     # check if time limit exceeded
                        if (time.time() - t0) / 3600 > tstop:
                            print 'Time limit exceeded...'
                            fileout.close
                            return 0
                    try: 
                        results = ev.imagedmt0(i, j, show=show, bgwindow=bgwindow, clean=1, mode=mode)
                        if mode == 'clean':
                            peak, epeak, off_ra, eoff_ra, off_dec, eoff_dec = results
                        elif mode == 'dirty':  # need to develop more... needs to be ~10ms processing per int!
                            peak, epeak = results
                            epeak = n.median(bgepeak)
                        if peak/epeak >= sig:
                            print '\tDetection!'
                            if mode == 'clean':
                                print >> fileout, file, nskip, nints, (i, j), 'Peak, (sig),  RA, Dec: ', peak, epeak, '(', peak/epeak, ')  ', off_ra, eoff_ra, off_dec, eoff_dec
                            elif mode == 'dirty':
                                print >> fileout, file, nskip, nints, (i, j), 'Peak, (sig): ', peak, epeak, '(', peak/epeak, ')'
                            # save all results (v1.0 pickle format)
                            pklout = open(pathout + string.join(file.split('.')[:-1], '.') + '.' + str(nskip) + '-dm' + str(i) + 't' + str(j) + '.pkl', 'wb')
                            pickle.dump((file, nskip, nints, n.array([i]), ev.dmarr[i], n.array([j]), n.array([peak/epeak])), pklout)
                            pklout.close()
                    except:
                        continue
                if mode == 'dirty':
                    print >> fileout, '   Finished ', file, nskip, nints, (i, j), '. Noise = ', n.median(bgepeak)

        fileout.close


def pulse_search_uvfit(fileroot, pathin, pathout, nints=12000, sig=5.0, show=0, edge=0):
    """
    Searches for pulses by imaging dedispersed trials.
    """

    maxints = 131000  # biggest file in integrations
    bgwindow = 10  # where bg subtraction is made

    filelist = []

# for crab 201103
#    for i in [0,1,4,5,6,7,8,9]:
#        filelist.append(string.join(fileroot.split('.')[:-1]) + '_0' + str(i) + '.mir')
#    for i in range(0,8):
#        filelist.append(string.join(fileroot.split('.')[:-1]) + '_1' + str(i) + '.mir')

# for crab 190348
#    for i in [0,1,2,3,4,5,6,7,8,9]:
#        filelist.append(string.join(fileroot.split('.')[:-1]) + '_0' + str(i) + '.mir')
#    for i in [1,3,4,5,6,7,8,9]:
#        filelist.append(string.join(fileroot.split('.')[:-1]) + '_1' + str(i) + '.mir')
#    for i in [0,1,2,3]:
#        filelist.append(string.join(fileroot.split('.')[:-1]) + '_2' + str(i) + '.mir')

# for b0329 173027
    for i in [0,1,2,3,4,5,6,7,8,9]:
        filelist.append(string.join(fileroot.split('.')[:-1]) + '_0' + str(i) + '.mir')
    for i in [0,1,2,6,7,8,9]:
        filelist.append(string.join(fileroot.split('.')[:-1]) + '_1' + str(i) + '.mir')
    for i in [0,1,2,3]:
        filelist.append(string.join(fileroot.split('.')[:-1]) + '_2' + str(i) + '.mir')

# for m31 154202
#    for i in [0,1,2,3,4,5,6,7,8,9]:
#        filelist.append(string.join(fileroot.split('.')[:-1]) + '_0' + str(i) + '.mir')
#    for i in [0,1,2,3,4,5,6]:
#        filelist.append(string.join(fileroot.split('.')[:-1]) + '_1' + str(i) + '.mir')

# hack to search single file for pulses
    filelist = [fileroot]
        
    print 'Looping over filelist: ', filelist
    for file in filelist:
        fileout = open(pathout + string.join(file.split('.')[:-1], '.') + '.txt', 'a')

        for nskip in [0]:
            print 'Starting file %s with nskip %d' % (file, nskip)

            # load data
            ev = evla(pathin + file, nints=nints, nskip=nskip)

            # dedisperse
            for i in range(len(ev.dmarr)):
#                for j in range(1200,1300):
                for j in range(len(ev.reltime)):
                    try: 
                        peak, epeak, off_ra, eoff_ra, off_dec, eoff_dec = ev.uvfitdmt0(i,j, show=show, bgwindow=bgwindow)
#                        print >> fileout, file, nskip, nints, (i, j), 'Peak, RA, Dec: ', peak, epeak, off_ra, eoff_ra, off_dec, eoff_dec

                        if peak/epeak >= sig:
                            print '\tDetection!'
                            print >> fileout, file, nskip, nints, (i, j), 'Peak, (sig),  RA, Dec: ', peak, epeak, '(', peak/epeak, ')  ', off_ra, eoff_ra, off_dec, eoff_dec
                            # save all results (v1.0 pickle format)
                            pklout = open(pathout + string.join(file.split('.')[:-1], '.') + '.' + str(nskip) + '-dm' + str(i) + 't' + str(j) + '.pkl', 'wb')
                            pickle.dump((file, nskip, nints, n.array([i]), ev.dmarr[i], n.array([j]), n.array([peak/epeak])), pklout)
                            pklout.close()
                    except:
                        continue

        fileout.close


def process_pickle(filename, pathin, mode='image'):
    """Processes a pickle file to produce a spectrum of a candidate pulse.
    mode tells whether to produce dmt0 plot ('dmt0'), a spectrogram ('spec'), 
    image the dm track ('image'), or write visibilities to a file ('dump')
    TO DO:  (maybe) modify format of pickle file.
    """

    file = open(filename, 'rb')
    dump = pickle.load(file)
    name = dump[0]
    nints = dump[2]
    nintskip = dump[1]
    dmbinarr = dump[3]
    dmarr = dump[4]
    tbinarr = dump[5]
    snrarr = dump[6]
    if snrarr[0] <= 1:  # reim mode has p-values for normality, which should be small when there is a pulse
        peaktrial = n.where(snrarr == min(snrarr))[0][0]
    else:
        peaktrial = n.where(snrarr == max(snrarr))[0][0]

    bgwindow = 4

    print 'Loaded pickle file for %s plot of %s' % (mode, name)
    print 'Has peaks at DM = ', dmarr, ' with sig ', snrarr
    print 'Using ', snrarr[peaktrial]

    if len(dmbinarr) >= 1:
#    if (len(dmbinarr) >= 1) & (snrarr[peaktrial] > 7.):
        print 'Grabbing %d ints at %d' % (nints, nintskip)
        ev = evla(pathin + name, nints=nints, nskip=nintskip)
        track = ev.dmtrack(dm=ev.dmarr[int(dmbinarr[peaktrial])], t0=ev.reltime[tbinarr[peaktrial]], show=0)

        if mode == 'spec':  # just show spectrum
            # write out bg-subbed track, read back in to fit spectrum
            for trial in range(len(dmbinarr)):
                track = ev.dmtrack(dm=ev.dmarr[int(dmbinarr[trial])], t0=ev.reltime[tbinarr[trial]], show=0)

               # plot track and spectrogram
                p.figure(1)
#                p.plot(ev.reltime[track[0]], track[1], 'w*')
                p.plot(track[0], track[1], 'w*')
            ev.spec(save=0)
            p.show()
        elif mode == 'dmt0':
            ev.makedmt0()
            peaks, peakssig = ev.peakdmt0()
#            p.plot(ev.reltime[bgwindow], ev.dmarr[dmbinarr[peaktrial]], '*' )   # not working?
            ev.plotdmt0(save=1)
        elif mode == 'image':
            immode = 'clean'
            results = ev.imagedmt0(int(dmbinarr[peaktrial]), tbinarr[peaktrial], bgwindow=bgwindow, clean=1, mode=immode, show=1)
            if immode == 'clean':
                peak, epeak, off_ra, eoff_ra, off_dec, eoff_dec = results
                print peak, epeak, off_ra, eoff_ra, off_dec, eoff_dec
            elif immode == 'dirty':
                peak, std = results
                print peak, std
        elif mode == 'uvfit':
            uvmode = 'elsegrid'
            if uvmode == 'grid':
                ev.uvfitdmt02(dmbinarr[peaktrial], tbinarr[peaktrial], mode='lsgrid')
            else:
                uvmode = 'grid'
                results = ev.uvfitdmt0(dmbinarr[peaktrial], tbinarr[peaktrial], mode='grid')
#                peak, epeak, off_ra, eoff_ra, off_dec, eoff_dec = results
#                print peak, epeak, off_ra, eoff_ra, off_dec, eoff_dec
        elif mode == 'dump':
            # write dmtrack data out for imaging
            ev.writetrack(track)
        elif mode == 'reim':
            datasub = ev.tracksub(dmbinarr[peaktrial], tbinarr[peaktrial], bgwindow=bgwindow)
            ev.data = datasub
            ev.plotreim()
        elif mode == 'imsearch':
            ev.imsearch(dmbinarr[peaktrial], tbinarr[peaktrial], nints, sig=7.0)
        else:
            print 'Mode not recognized'
    else:
        print 'No significant detection.  Moving on...'

    file.close()


if __name__ == '__main__':
    """From the command line, evlavis can either load pickle of candidate to interact with data, or
    it will search for pulses blindly.
    """

    print 'Greetings, human.'
    print ''

    fileroot = 'crab_12ms_rr1.mir'
    pathin = './'
    pathout = 'crab_rr1_ms/'
    edge = 10

    if len(sys.argv) == 1:
        # if no args, search for pulses
        print 'Searching for pulses...'
        try:
#            cProfile.run('pulse_search_uvfit(fileroot=fileroot, pathin=pathin, pathout=pathout, nints=2000, edge=edge)')
#            pulse_search_image(fileroot=fileroot, pathin=pathin, pathout=pathout, nints=2000, edge=edge, mode='dirty', sig=7.0)
            pulse_search_phasecenter(fileroot=fileroot, pathin=pathin, pathout=pathout, nints=2000, edge=edge)
#            pulse_search_reim(fileroot=fileroot, pathin=pathin, pathout=pathout, nints=2000, edge=edge)
        except AttributeError:
            exit(0)
    elif len(sys.argv) == 2:
        # if pickle, then plot data or dm search results
        print 'Assuming input file is pickle of candidate...'
        process_pickle(sys.argv[1], pathin=pathin, mode='image')
    elif len(sys.argv) == 7:
        # if pickle, then plot data or dm search results
        print 'Time limited searching for pulses... with %s, %s, %s' % (fileroot, pathin, pathout)
        try:
            dmrange = [int(sys.argv[4])]    # only works for single dm
            nstart = int(sys.argv[5])  # start integration
            tstop = float(sys.argv[6])  # run time in hours
            pulse_search_image(fileroot=sys.argv[1], pathin=sys.argv[2], pathout=sys.argv[3], nints=20000, edge=edge, mode='dirty', sig=7.0, dmrange=dmrange, nstart=nstart, tstop=tstop)
        except AttributeError:
            exit(0)
    elif len(sys.argv) == 6:
        # if full spec of trial, image it
        print 'Imaging DM trial...'
        file = sys.argv[1]
        nskip = int(sys.argv[2])
        nints = int(sys.argv[3])
        dmbin = int(sys.argv[4])
        t0bin = int(sys.argv[5])
        ev = evla(file, nints=nints, nskip=nskip)
        peak, epeak, off_ra, eoff_ra, off_dec, eoff_dec = ev.imagedmt0(dmbin, t0bin, show=1)
        print file, nskip, nints, (dmbin, t0bin), '::Peak, RA, Dec:: ', peak, epeak, off_ra, eoff_ra, off_dec, eoff_dec
