import tpipe
import numpy as n
import pylab as p
import string, os, pickle
from casa import ms
from casa import quanta as qa

class efficient():
    """ Set of efficient functions that override tpipe class methods.
    """

    def __init__(self):
        pass

    def prep(self):
        """ Sets up tracks used to speed up dedispersion code.
        Has the option to delete raw data and flags to save memory.
        """
        print
        print 'Filtering data to masked array...'
        self.data = n.ma.masked_array(self.rawdata[:self.nints,:, self.chans,:], self.flags[:self.nints,:, self.chans,:] == 0)   # mask of True for flagged data (flags=0 in tpipe, which is flags=False in Miriad and flags=True in MS)

        # apply telcal solutions
        if self.telcalfile:    # apply telcal solutions
            print 'telcal only works for a single spw for now. telcal flags not used.'
            sols = applytelcal.solutions(self.telcalfile, 1e9*self.freq_orig)

            for polstr in self.selectpol:
                if polstr == 'RR':
                    pol = 0
                elif polstr == 'LL':
                    pol = 1

            try:
                sols.setselection(telcalibrator, time1[0]/(24*3600), pol, verbose=1)   # chooses solutions closest in time that match pol and source name
                sols.apply(data1, self.blarr, pol)
            except:
                pass

        self.dataph = (self.data.mean(axis=3).mean(axis=1)).real   #dataph is summed and detected to form TP beam at phase center, multi-pol
        self.min = self.dataph.min()
        self.max = self.dataph.max()
        self.reltime = self.reltime[:self.nints]
        print 'Shape of data:'
        print self.data.shape

        # set up triples and arrays for bispectrum considering flagged baselines (only having zeros).
        triples = self.make_triples()
        meanbl = self.data.mean(axis=2).mean(axis=0)   # find bls with no zeros in either pol to ignore in triples
        self.triples = triples[n.all(meanbl[triples][:,0] != 0j, axis=1) & n.all(meanbl[triples][:,1] != 0j, axis=1) & n.all(meanbl[triples][:,2] != 0j, axis=1) == True]   # only take triples if both pols are good. may be smaller than set for an individual pol

        # initialize data state
        self.datadelay = self.calc_delay(0)

        # set up ur tracks (lol)
        self.dmtrack0 = {}
        self.twidths = {}
        for dmbin in xrange(len(self.dmarr)):
            self.dmtrack0[dmbin] = self.dmtrack(self.dmarr[dmbin],0)  # track crosses high-freq channel in first integration
            self.twidths[dmbin] = 0
            for k in self.dmtrack0[dmbin][1]:
                self.twidths[dmbin] = max(self.twidths[dmbin], len(n.where(n.array(self.dmtrack0[dmbin][1]) == k)[0]))

    def calc_delay(self, dm):
        """ Function to calculate delay for each channel, given dm
        """

        return n.round((4.2e-3 * dm * (self.freq**(-2) - self.freq[len(self.chans)-1]**(-2)))/self.inttime,0).astype(int)

    def dedisperse(self, dm):
        """ Dedisperse each channel, but only if data not already shifted.
        """

        newdelay = self.calc_delay(dm)

        relativedelay = newdelay - self.datadelay

        for i in xrange(len(self.chans)):
            if relativedelay[i] != 0:
                self.data[:,:,i,:] = n.roll(self.data[:,:,i,:], -relativedelay[i], axis=0)

        self.datadelay = newdelay

    def make_bispectra(self, stokes=''):
        """ Makes numpy array of bispectra for each integration.
        stokes parameter is ignored.
        """

        bisp = lambda d: d[:,:,0] * d[:,:,1] * n.conj(d[:,:,2])    # bispectrum for data referenced by triple (data[:,triples])

        bispectra = n.zeros((len(self.dmarr), len(self.data), len(self.triples)), dtype='complex')

        # iterate over dm trials
        for dmbin in xrange(len(self.dmarr)):
            self.dedisperse(self.dmarr[dmbin])
            bispectra[dmbin] = bisp(self.data.mean(axis=2).mean(axis=2)[:, self.triples])


            edgeints = range(self.nints - (n.max(self.calc_delay(self.dmarr[dmbin]))-1), self.nints)  # -1 corrects for min of max delay
            bispectra[dmbin, edgeints] = 0j
            print 'dedispersed for ', self.dmarr[dmbin]

        self.dedisperse(0)   # restore data to original locations

        self.bispectra = n.ma.masked_array(bispectra, bispectra == 0j)

    def make_bispectra2(self):
        """ Makes numpy array of bispectra for each integration.
        """

        bisp = lambda d: n.real(d[:,:,0] * d[:,:,1] * n.conj(d[:,:,2]))    # bispectrum for data referenced by triple (data[:,triples])

        return bisp(self.data.mean(axis=2).mean(axis=2)[:, self.triples])

    def phaseshift(self, l1, m1):
        ang = lambda dl, dm, u, v, freq: (dl*n.outer(u, freq/freq[0]) + dm*n.outer(v, freq/freq[0]))

        dl = l1 - self.l0[0]
        dm = m1 - self.m0[0]

        for i in xrange(self.data.shape[0]):    # iterate over integrations
            for pol in xrange(self.data.shape[3]):    # iterate over pols
                self.data[i,:,:,pol] = self.data[i,:,:,pol] * n.exp(-2j*n.pi*ang(dl, dm, self.u[i], self.v[i], self.freq))

        self.l0 = l1 * n.ones(self.nints)
        self.m0 = m1 * n.ones(self.nints)

class pipe_simdisp2(efficient, tpipe.pipe_simdisp2):

    def __init__(self, profile='default', nints=256, inttime=0.01, freq=1.4, bw=0.128, Q=1., array='vla_a', **kargs):
        self.set_profile(profile=profile)
        self.set_params(**kargs)
        self.simulate(nints, inttime, self.chans, freq, bw, Q=Q, array=array)
        self.prep()

    def simulate(self, nints, inttime, chans, freq, bw, Q=1., array='vla10'):
        """ Simulates data
        array is the name of array config, nints is number of ints, inttime is integration duration, chans, freq, bw (in GHz) as normal.
        array can be 'vla_d', 'vla_a', and 'vla10', the latter is the first 10 of vla_d.
        """

        self.file = 'sim'
        self.chans = chans
        self.nints = nints
        self.nchan = len(chans)
        print 'Initializing nchan:', self.nchan
        self.sfreq = freq    # in GHz
        self.sdf = bw/self.nchan
        self.npol = 1
        self.freq = self.sfreq + self.sdf * n.arange(self.nchan)
        self.inttime = inttime   # in seconds
        self.reltime = inttime*n.arange(nints)

        # define relative phase center for each integration
        self.l0 = n.zeros(self.nints)
        self.m0 = n.zeros(self.nints)

        # antennas and baselines
        vla_d = 1e3*n.array([[ 0.00305045,  0.03486681],  [ 0.00893224,  0.10209601],  [ 0.01674565,  0.19140365],  [ 0.02615514,  0.29895461],  [ 0.03696303,  0.42248936],  [ 0.04903413,  0.56046269],  [ 0.06226816,  0.7117283 ],  [ 0.07658673,  0.87539034],  [ 0.09192633,  1.05072281],  [ 0.02867032, -0.02007518],  [ 0.08395162, -0.05878355],  [ 0.1573876 , -0.11020398],  [ 0.24582472, -0.17212832],  [ 0.347405  , -0.2432556 ],  [ 0.46085786, -0.32269615],  [ 0.58524071, -0.40978995],  [ 0.71981691, -0.50402122],  [ 0.86398948, -0.60497195],  [-0.03172077, -0.01479164],  [-0.09288386, -0.04331245],  [-0.17413325, -0.08119967],  [-0.27197986, -0.12682629],  [-0.38436803, -0.17923376],  [-0.509892  , -0.23776654],  [-0.64750886, -0.30193834],  [-0.79640364, -0.37136912],  [-0.95591581, -0.44575086]])
        vla_a = 27*1e3*n.array([[ 0.00305045,  0.03486681],  [ 0.00893224,  0.10209601],  [ 0.01674565,  0.19140365],  [ 0.02615514,  0.29895461],  [ 0.03696303,  0.42248936],  [ 0.04903413,  0.56046269],  [ 0.06226816,  0.7117283 ],  [ 0.07658673,  0.87539034],  [ 0.09192633,  1.05072281],  [ 0.02867032, -0.02007518],  [ 0.08395162, -0.05878355],  [ 0.1573876 , -0.11020398],  [ 0.24582472, -0.17212832],  [ 0.347405  , -0.2432556 ],  [ 0.46085786, -0.32269615],  [ 0.58524071, -0.40978995],  [ 0.71981691, -0.50402122],  [ 0.86398948, -0.60497195],  [-0.03172077, -0.01479164],  [-0.09288386, -0.04331245],  [-0.17413325, -0.08119967],  [-0.27197986, -0.12682629],  [-0.38436803, -0.17923376],  [-0.509892  , -0.23776654],  [-0.64750886, -0.30193834],  [-0.79640364, -0.37136912],  [-0.95591581, -0.44575086]])

        if array == 'vla_d':
            antloc = vla_d
        elif array == 'vla_a':
            antloc = vla_a
        elif array == 'vla10':
            antloc = vla_d[5:15]    # 5:15 choses inner part  of two arms
        elif array == 'mwa':
            antloc=1e3*n.array([[0,0]])
        self.nants = len(antloc)
        print 'Initializing nants:', self.nants
        blarr = []; u = []; v = []; w = []
        for i in range(1, self.nants+1):
            for j in range(i, self.nants+1):
                blarr.append([i,j])
                u.append(antloc[i-1][0] - antloc[j-1][0])   # in meters (like MS, fwiw)
                v.append(antloc[i-1][1] - antloc[j-1][1])
                w.append(0.)
        self.blarr = n.array(blarr)

        self.nbl = len(self.blarr)
        self.u = n.zeros((nints,self.nbl),dtype='float64')
        self.v = n.zeros((nints,self.nbl),dtype='float64')
        self.w = n.zeros((nints,self.nbl),dtype='float64')
        print 'Initializing nbl:', self.nbl
        self.ants = n.unique(self.blarr)
        self.nskip = 0

        # no earth rotation yet
        for i in range(nints):
            self.u[i] = n.array(u) * self.freq[0] * (1e9/3e8)
            self.v[i] = n.array(v) * self.freq[0] * (1e9/3e8)
            self.w[i] = n.array(w) * self.freq[0] * (1e9/3e8)

        # simulate data
        self.rawdata = n.zeros((nints,self.nbl,self.nchan,self.npol),dtype='complex64')
        self.flags = n.ones((nints,self.nbl,self.nchan,self.npol),dtype='bool')
        self.rawdata.real = Q * n.sqrt(self.nchan) * n.random.randn(nints,self.nbl,self.nchan,self.npol)   # normal width=1 after channel mean
        self.rawdata.imag = Q * n.sqrt(self.nchan) * n.random.randn(nints,self.nbl,self.nchan,self.npol)

        # print summary info
        print
        print 'Shape of raw data, time:'
        print self.rawdata.shape, self.reltime.shape

class pipe_msdisp2(efficient, tpipe.pipe_msdisp2):

    def __init__(self, file, profile='default', nints=256, nskip=0, spw=[-1], selectpol=['RR','LL'], scan=0, datacol='corrected_data', telcalfile='', telcalibrator='', **kargs):
        self.telcalfile=telcalfile
        self.telcalibrator=telcalibrator
        self.set_profile(profile=profile)
        self.set_params(**kargs)
        self.read(file=file, nints=nints, nskip=nskip, spw=spw, selectpol=selectpol, scan=scan, datacol=datacol)
        self.prep()

    def read(self, file, nints, nskip, spw, selectpol, scan, datacol):
        """ Reads in Measurement Set data using CASA.
        spw is list of subbands. zero-based.
        Scan is zero-based selection based on scan order, not actual scan number.
        selectpol is list of polarization strings (e.g., ['RR','LL'])
        """
        self.file = file
        self.scan = scan
        self.nints = nints

        # get spw info. either load pickled version (if found) or make new one
        pklname = string.join(file.split('.')[:-1], '.') + '_init.pkl'
#        pklname = pklname.split('/')[-1]  # hack to remove path and write locally
        if os.path.exists(pklname):
            print 'Pickle of initializing info found. Loading...'
            pkl = open(pklname, 'r')
            try:
                (self.npol_orig, self.nbl, self.blarr, self.inttime, spwinfo, scansummary) = pickle.load(pkl)
            except EOFError:
                print 'Bad pickle file. Exiting...'
                return 1
# old way, casa 3.3?
#            scanlist = scansummary['summary'].keys()
#            starttime_mjd = scansummary['summary'][scanlist[scan]]['0']['BeginTime']
# new way, casa 4.0?
            scanlist = scansummary.keys()
            starttime_mjd = scansummary[scanlist[scan]]['0']['BeginTime']
            self.nskip = int(nskip*self.nbl)    # number of iterations to skip (for reading in different parts of buffer)
            self.npol = len(selectpol)
        else:
            print 'No pickle of initializing info found. Making anew...'
            pkl = open(pklname, 'wb')
            ms.open(self.file)
            spwinfo = ms.getspectralwindowinfo()
            scansummary = ms.getscansummary()

# original (general version)
#            scanlist = scansummary['summary'].keys()
#            starttime_mjd = scansummary['summary'][scanlist[scan]]['0']['BeginTime']
#            starttime0 = qa.getvalue(qa.convert(qa.time(qa.quantity(starttime_mjd+0/(24.*60*60),'d'),form=['ymd'], prec=9), 's'))
#            stoptime0 = qa.getvalue(qa.convert(qa.time(qa.quantity(starttime_mjd+0.5/(24.*60*60), 'd'), form=['ymd'], prec=9), 's'))

# for casa 4.0 (?) and later
            scanlist = scansummary.keys()
            # set time info
            self.inttime = scansummary[scanlist[scan]]['0']['IntegrationTime']
            self.inttime0 = self.inttime
            print 'Initializing integration time (s):', self.inttime
            starttime_mjd = scansummary[scanlist[scan]]['0']['BeginTime']
            starttime0 = qa.getvalue(qa.convert(qa.time(qa.quantity(starttime_mjd+0/(24.*60*60),'d'),form=['ymd'], prec=9)[0], 's'))[0]
            stoptime0 = qa.getvalue(qa.convert(qa.time(qa.quantity(starttime_mjd+self.inttime/(24.*60*60), 'd'), form=['ymd'], prec=9)[0], 's'))[0]

            ms.selectinit(datadescid=0)  # initialize to initialize params
            selection = {'time': [starttime0, stoptime0]}
            ms.select(items = selection)
            da = ms.getdata([datacol, 'axis_info'], ifraxis=True)
            ms.close()

            self.npol_orig = da[datacol].shape[0]
            self.nbl = da[datacol].shape[2]
            print 'Initializing nbl:', self.nbl

            # good baselines
            bls = da['axis_info']['ifr_axis']['ifr_shortname']
            self.blarr = n.array([[int(bls[i].split('-')[0]),int(bls[i].split('-')[1])] for i in xrange(len(bls))])
            self.nskip = int(nskip*self.nbl)    # number of iterations to skip (for reading in different parts of buffer)

            pickle.dump((self.npol_orig, self.nbl, self.blarr, self.inttime, spwinfo, scansummary), pkl)
        pkl.close()

        self.ants = n.unique(self.blarr)
        self.nants = len(n.unique(self.blarr))
        self.nants0 = len(n.unique(self.blarr))
        print 'Initializing nants:', self.nants
        self.npol = len(selectpol)
        print 'Initializing %d of %d polarizations' % (self.npol, self.npol_orig)

        # set desired spw
        if (len(spw) == 1) & (spw[0] == -1):
            #            spwlist = spwinfo['spwInfo'].keys()    # old way
            spwlist = spwinfo.keys()    # new way
        else:
            spwlist = spw

        self.freq_orig = n.array([])
        for spw in spwlist:
            # new way
            nch = spwinfo[str(spw)]['NumChan']
            ch0 = spwinfo[str(spw)]['Chan1Freq']
            chw = spwinfo[str(spw)]['ChanWidth']
            self.freq_orig = n.concatenate( (self.freq_orig, (ch0 + chw * n.arange(nch)) * 1e-9) )
# old way
#            nch = spwinfo['spwInfo'][str(spw)]['NumChan']
#            ch0 = spwinfo['spwInfo'][str(spw)]['Chan1Freq']
#            chw = spwinfo['spwInfo'][str(spw)]['ChanWidth']

        self.freq = self.freq_orig[self.chans]
        self.nchan = len(self.freq)
        print 'Initializing nchan:', self.nchan

        # set requested time range based on given parameters
        timeskip = self.inttime*nskip
# new way        
        starttime = qa.getvalue(qa.convert(qa.time(qa.quantity(starttime_mjd+timeskip/(24.*60*60),'d'),form=['ymd'], prec=9)[0], 's'))[0]
        stoptime = qa.getvalue(qa.convert(qa.time(qa.quantity(starttime_mjd+(timeskip+nints*self.inttime)/(24.*60*60), 'd'), form=['ymd'], prec=9)[0], 's'))[0]
        print 'First integration of scan:', qa.time(qa.quantity(starttime_mjd,'d'),form=['ymd'],prec=9)[0]
        print
# new way
        print 'Reading scan', str(scanlist[scan]) ,'for times', qa.time(qa.quantity(starttime_mjd+timeskip/(24.*60*60),'d'),form=['hms'], prec=9)[0], 'to', qa.time(qa.quantity(starttime_mjd+(timeskip+nints*self.inttime)/(24.*60*60), 'd'), form=['hms'], prec=9)[0]

        # read data into data structure
        ms.open(self.file)
        ms.selectinit(datadescid=spwlist[0])  # reset select params for later data selection
        selection = {'time': [starttime, stoptime]}
        ms.select(items = selection)
        print 'Reading %s column, SB %d, polarization %s...' % (datacol, spwlist[0], selectpol)
        ms.selectpolarization(selectpol)
        da = ms.getdata([datacol,'axis_info','u','v','w','flag'], ifraxis=True)
        u = da['u']; v = da['v']; w = da['w']
        if da == {}:
            print 'No data found.'
            return 1
        newda = n.transpose(da[datacol], axes=[3,2,1,0])  # if using multi-pol data.
        flags = n.transpose(da['flag'], axes=[3,2,1,0])
        if len(spwlist) > 1:
            for spw in spwlist[1:]:
                ms.selectinit(datadescid=spw)  # reset select params for later data selection
                ms.select(items = selection)
                print 'Reading %s column, SB %d, polarization %s...' % (datacol, spw, selectpol)
                ms.selectpolarization(selectpol)
                da = ms.getdata([datacol,'axis_info','flag'], ifraxis=True)
                newda = n.concatenate( (newda, n.transpose(da[datacol], axes=[3,2,1,0])), axis=2 )
                flags = n.concatenate( (flags, n.transpose(da['flag'], axes=[3,2,1,0])), axis=2 )
        ms.close()

        # Initialize more stuff...
        self.nschan0 = self.nchan

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

        # Assumes MS files store uvw in meters. Corrects by mean frequency of channels in use.
        self.u = u.transpose() * self.freq_orig[0] * (1e9/3e8)
        self.v = v.transpose() * self.freq_orig[0] * (1e9/3e8)
        self.w = w.transpose() * self.freq_orig[0] * (1e9/3e8)

        # set integration time and time axis
        ti = da['axis_info']['time_axis']['MJDseconds']
        self.reltime = ti - ti[0]

        # define relative phase center for each integration
        self.l0 = n.zeros(self.nints)
        self.m0 = n.zeros(self.nints)

        self.rawdata = newda
        self.flags = n.invert(flags)             # tests show that MS has opposite flag convention as Miriad! using complement of MS flag in tpipe.
        print 'Shape of raw data, time:'
        print self.rawdata.shape, self.reltime.shape


