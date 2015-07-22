#! /usr/bin/env python

"""
bisp_mbsearch.py --- tpipe pipeline to do bispectrum search over multiple dm and delay beams
Summarizes candidates for time segment and optionally images/beamforms to localize and measure spectrum of candidate.
"""

import leanpipe as tpipe
import pickle, sys
import numpy as n
import pylab as p

class search():
    """ Class that defines search process.
    """

    def __init__(self, datatype='', datafile='', startints=[0], nints=128, timescales=[1], filtershape='b', dmarr=[0.], fwhmsurvey=0.5, fwhmfield=0.5, selectpol=['RR','LL'], stokes='prebisp', scan=0, spw=[0], chans=range(5,59), bgwindow=4, sigma=5, size=300000, res=1000):
        """ Initialize search pipeline. Parameters of search defined here, including search domain (ints, dm, dt, lm).
        searches over list of startints of size ints
        selectpol, scan, spw are only used for ms data
        """

        if (datatype == '') | (datafile == ''):
            print 'Must specify data and data type.'
            return

# data reading
        self.datatype = datatype
        self.datafile = datafile
        self.fileroot = self.datafile.split('.')[0]
# search domain
        self.startints = startints
        self.nints = nints    # number of ints to read at a time
        self.scan = scan
        self.spw = spw
        self.chans = chans
        self.timescales = timescales    # number of ints in boxcar filter
        self.filtershape = filtershape
        self.dmarr = dmarr
        self.fwhmsurvey = fwhmsurvey  # field of view to tile (fwhm)
        self.fwhmfield = fwhmfield   # delay beam size (fwhm)
        self.selectpol = selectpol
        self.delaycenters = self.calc_hexcenters()
        self.stokes = stokes
# search params
        self.bgwindow = bgwindow
        self.sigma = sigma
# imaging
        self.size = size   # aipy param determining image pixel scale
        self.res = res  # aipy param determining image field of view

    def description(self):
        if self.datatype=='ms':
            return '**Bispectrum search pipe parameters** \nData: %s.%s \nTime: %s, %s, %s \nSelection: %s, %s, %s, %s, %s, %s \nDM: %s \nBeams: %s, %s, \n %s \nSearch: %s %s \nImaging: %s %s' % (self.fileroot, self.datatype, self.startints, self.nints, self.scan, self.spw, self.chans, self.selectpol, self.stokes, self.timescales, self.filtershape, self.dmarr, self.fwhmsurvey, self.fwhmfield, self.delaycenters, self.bgwindow, self.sigma, self.size, self.res)
        else:
            return '**Bispectrum search pipe parameters** \nData: %s.%s \nTime: %s, %s, %s, %s, %s, %s \nDM: %s \nBeams: %s, %s, \n %s \nSearch: %s %s \nImaging: %s %s' % (self.fileroot, self.datatype, self.startints, self.nints, self.timescales, self.filtershape, self.chans, self.stokes, self.dmarr, self.fwhmsurvey, self.fwhmfield, self.delaycenters, self.bgwindow, self.sigma, self.size, self.res)

    def calc_hexcenters(self, show=0):
        """ Tile a large circular area with a small circular areas. sizes are assumed to be fwhm. assumes flat sky.
        """

        large = self.fwhmsurvey
        small = self.fwhmfield
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

        if show:
            p.plot(delaycenters[:,0], delaycenters[:,1], '.')
            p.xlabel('x Offset (deg)',fontsize=12,fontweight="bold")
            p.ylabel('y Offset (deg)',fontsize=12,fontweight="bold")
            p.title('Delay centers')

        if len(delaycenters) == 1:
            plural = ''
        else:
            plural = 's'
        print 'For a search area of %.2f and delay beam of %.2f, we will use %d delay beam%s' % (self.fwhmsurvey, self.fwhmfield, len(delaycenters), plural)
        return delaycenters

    def do_search(self):
        """ Bispectrum transient search over all delay beams.
        """

        self.masterloc = {}           # "location" of candidate: dt, integration (over all data in pipe), dmbin
        self.masterprop = {}          # properties of candidate: snr and std of bispectra for dedispersed trial
        self.masterdata = {}          # data around the candidate: bisp lightcurve and spectrogram
        self.masternoise = {}         # noise per bl as measured by detect_bispectra
        for (ra,dec) in self.delaycenters:
            self.masterloc[(ra,dec)] = []
            self.masterprop[(ra,dec)] = []
            self.masterdata[(ra,dec)] = []
            self.masternoise[(ra,dec)] = []

        for startint in self.startints:
            if self.datatype == 'mir':
                self.obs = tpipe.pipe_mirdisp2(file=self.datafile, nints=self.nints, nskip=startint, profile='default', dmarr=self.dmarr, chans=self.chans)
            elif self.datatype == 'ms':
                self.obs = tpipe.pipe_msdisp2(file=self.datafile, nints=self.nints, nskip=startint, profile='default', selectpol=self.selectpol, scan=self.scan, spw=self.spw, chans=self.chans, dmarr=self.dmarr, datacol='corrected_data')
            elif self.datatype == 'sim':
                self.obs = tpipe.pipe_simdisp2(array='vla_a', nints=self.nints, profile='default', dmarr=self.dmarr, chans=range(64))
#                self.obs.add_transient(0., 0., 1., n.random.randint(self.nints))
            else:
                print 'Unknown datatype!'
                return

            # accumulate running reltime
            if startint == self.startints[0]:
                self.tarr = self.obs.reltime
            else:
                self.tarr = n.append(self.tarr, self.obs.reltime + self.tarr[-1] + self.obs.inttime)

            for dtind in range(len(self.timescales)):

                # after first pass, need to reset the working data so effect of time_filter is as expected
#                if self.timescales.index(dt) > 0:
                if dtind > 0:
                    print 'Resetting data to unfiltered version...'
                    self.obs.prep()

                self.obs.time_filter(self.timescales[dtind], self.filtershape, bgwindow=self.bgwindow)

                for (ra,dec) in self.delaycenters:
                    print
                    print 'Searching int %d, dt %d, and delay beam (ra, dec) = %.3f, %.3f' % (startint, self.timescales[dtind], ra, dec)
                    l1 = n.radians(ra); m1 = n.radians(dec)
                    if ( (l1 != self.obs.l0[0]) | (m1 != self.obs.m0[0]) ):
                        dl = l1 - self.obs.l0[0]; dm = m1 - self.obs.m0[0]  # find phase shift for new location
                        self.obs.phaseshift(dl, dm)

                    # measure single noise for input to detect_bispectra
                    noiseperbl = self.obs.data.mean(axis=3).mean(axis=2).real.std()

                    # find candidates
                    self.obs.make_bispectra(stokes=self.stokes)
# add a way to cut candidates near edges (due to time_filter)?
                    (candsnr, candstd, cands) = self.obs.detect_bispectra(sigma=self.sigma, Q=noiseperbl, save=0)

                    # build master list of candidate location from bispectrum search
                    loclist = self.masterloc[(ra,dec)]
                    proplist = self.masterprop[(ra,dec)]
                    datalist = self.masterdata[(ra,dec)]
                    for i in range(len(cands)):
                        loclist.append( [dtind, (startint-self.startints[0])+cands[i][0], cands[i][1]] )    # candidate location
                        proplist.append( [candsnr[i], candstd[i]] )                  # candidate bisp properties
                        intmin = max(0, cands[i][0] - 64)                            # candidate lightcuve and spectrogram
                        intmax = min(cands[i][0] + 64, self.nints)
                        bamean = self.obs.bispectra.real.mean(axis=2)[cands[i][1], intmin:intmax]     # trim in time to save space
                        bastd = self.obs.bispectra.real.std(axis=2)[cands[i][1], intmin:intmax]
                        spec = self.obs.dataph[intmin:intmax]
                        datalist.append( (bamean, bastd, spec) )

                    self.masterloc[(ra,dec)] = loclist
                    self.masterprop[(ra,dec)] = proplist
                    self.masterdata[(ra,dec)] = datalist
                    self.masternoise[(ra,dec)] = (dtind, self.obs.Q)
                    print '%d candidates found in beam of %s at int %d' % (len(loclist), self.datafile, startint)
        
    def find_candidates(self, d_neighbor=2, island_size_min=1, island_snr_min=5., show=1, save=0):
        """ Visualize and define candidates to follow up.
        d_neighbor is the dm-time distance over which a neighbor is defined.
        island_size_min is the minimum number of detections requried to call an island a candidate
        island_snr_ min is the minimum snr requried to call an island a candidate
        returns the SNR peak of each island identified with d_neighbor.
        """

        neighbors = {}
        for (ra,dec) in self.delaycenters:
            neighbors[(ra,dec)] = []

        dmarr = n.array(self.dmarr)
        tarr = self.tarr

        # measure number of neighbors
        for beam in self.masterloc.keys():
            neighbors_lm = []
            cands_lm = self.masterloc[beam]
            for cand in cands_lm:
                nn = -1    # to avoid self-counting
                diff = n.array(cands_lm) - n.array(cand)
                for dd in diff:
                    if ( (n.abs(dd[0]) <= d_neighbor) & (n.abs(dd[1]) <= d_neighbor) & (n.abs(dd[2]) <= d_neighbor) ):  # if dmbin and int are within d_neighbor, then they are neighbors
                        nn = nn+1
                neighbors_lm.append(nn)
            neighbors[beam] = neighbors_lm

        # define islands in DM-time
        self.islands = {}
        for (ra,dec) in self.masterloc.keys():
            self.islands[(ra,dec)] = []
        for beam in self.masterloc.keys():
            cands_lm = self.masterloc[beam]
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

            self.islands[beam] = islands_lm
        print 'Identified DM-time islands.'

        # find snr peaks of islands including filters for island size and snr
        self.islandmaxloc = {}
        self.islandmaxsnr = {}
        for beam in self.islands.keys():
            islandmaxloc = []; islandmaxsnr = []
            islands_lm = self.islands[beam]
            for island in islands_lm:
                if len(island) >= island_size_min:
                    islandind = []
                    for i in range(len(island)):
                        islandind.append(n.where( (n.array(self.masterloc[beam])[:,0]==island[i,0]) & (n.array(self.masterloc[beam])[:,1]==island[i,1]) & (n.array(self.masterloc[beam])[:,2]==island[i,2]) )[0][0])
                    maxsnrind = n.where(n.array(self.masterprop[beam])[islandind,0] == n.max(n.array(self.masterprop[beam])[islandind,0]))
                    if n.array(self.masterprop[beam])[islandind][maxsnrind][0][0] > island_snr_min:
                        islandmaxloc.append(n.array(self.masterloc[beam])[islandind][maxsnrind][0])
                        islandmaxsnr.append(n.array(self.masterprop[beam])[islandind][maxsnrind][0][0])
            self.islandmaxloc[beam] = n.array(islandmaxloc)
            self.islandmaxsnr[beam] = n.array(islandmaxsnr)

        if show or save:
            cm = p.get_cmap('gist_rainbow')
            p.figure(1, figsize=(12,9))
            if self.datatype != 'ms':   # casapy doesn't get along with recent matplotlib stuff
                tl1 = p.subplot2grid( (3,3), (0,0), colspan=2, rowspan=2)   # dm vs. time
                tr1 = p.subplot2grid( (3,3), (0,2), rowspan=2)      # dm hist
                bl1 = p.subplot2grid( (3,3), (2,0), colspan=2)      # snr vs. time
                br1 = p.subplot2grid( (3,3), (2,2))           # snr hist
            elif self.datatype == 'ms':
                tl1 = p.subplot(221)
                tr1 = p.subplot(222)
                bl1 = p.subplot(223)
                br1 = p.subplot(224)

            beamind = 0
            for beam in self.masterloc.keys():  # iterate over beam candidate groups
                loc = n.array(self.masterloc[beam])
                prop = n.array(self.masterprop[beam])
                shiftind = float(beamind)/len(self.masterloc.keys())   # shift point location (and color) per beam for clarity
                if len(loc) == 0: break
                p.figure(1)
                p.subplot(tl1)
                p.scatter(tarr[loc[:,1]], dmarr[loc[:,2]] + 0.5*(dmarr[1]-dmarr[0])*shiftind, s=8*prop[:,0], facecolor='none', color=cm(shiftind), alpha=0.5, clip_on=False)
                for j in range(len(loc)):   # iterate over cands in group to plot one point each
                    if neighbors[beam][j] > 1:
                        p.text(tarr[loc[j,1]], dmarr[loc[j,2]] + 0.5*(dmarr[1]-dmarr[0])*shiftind, s=str(neighbors[beam][j]), horizontalalignment='center', verticalalignment='center', fontsize=9, color=cm(shiftind), alpha=0.5)
                        
                # plot island peaks
                p.scatter(tarr[self.islandmaxloc[beam][:,1]], dmarr[self.islandmaxloc[beam][:,2]] + 0.5*(dmarr[1]-dmarr[0])*shiftind, s=30*self.islandmaxsnr[beam], color=cm(shiftind), alpha=0.8, marker='+', clip_on=False)
                p.xticks([])
                axis_tl1 = p.axis()

                p.subplot(tr1)
                dms = self.dmarr
                dms.append(dms[-1]+(dms[-1]-dms[-2]))
                dmhist = n.array(dms) - (dms[-1]-dms[-2])/2.
                p.hist(dmarr[loc[:,2]], orientation='horizontal', color=cm(shiftind), histtype='step', alpha=0.5, bins=dmhist)
                p.yticks([])
                axis_tr1 = p.axis()
                p.axis([0, 1.1*axis_tr1[1], axis_tl1[2], axis_tl1[3]])

                p.subplot(bl1)
                p.scatter(tarr[loc[:,1]], prop[:,0], color=cm(shiftind), alpha=0.3, clip_on=False)
                p.scatter(tarr[self.islandmaxloc[beam][:,1]], self.islandmaxsnr[beam], color=cm(shiftind), marker='+', s=100, alpha=0.8, clip_on=False)
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
            p.xlabel('Time (s)',fontsize=14,fontweight="bold")
            p.ylabel('SNR',fontsize=14,fontweight="bold")
            p.subplot(br1)
#        p.ylabel('SNR',fontsize=12,fontweight="bold")
            p.xlabel('Count',fontsize=14,fontweight="bold")
            if self.datatype == 'mir':   # casapy doesn't get along with recent matplotlib stuff
                p.tight_layout()

            if save:
                p.savefig(self.figfile)
                if show:
                    p.show()
            else:
                if show:
                    p.show()

        return self.islandmaxloc    # return peaks of islands

    def candplot1(self, plotcands={}, save=0):
        """ Uses results of bispectrum results (calculated in do_search and find_candidates) to produce summary plot per candidate.
        This is a "tier 1" plot, since it only uses info available to bispectrum search algorithm.
        Requires do_search and find_candidates to be run first to build "master" objects.
        Candidate should be dictionary with keys of beam location and value of a list of 3-element lists (dtind, int, dmbin).
        save defines whether plots are also saved.
        """

        print 'Building Tier 1 candidate plot.'

        if n.array([len(plotcands[beam]) for beam in plotcands.keys()]).sum() == 0:
            print 'No candidates available...'
            return

        figi = 1
        for beam in plotcands.keys():
            if len(plotcands[beam]) == 0:
                continue
            if (n.ndim(plotcands[beam]) != 2) | (n.shape(plotcands[beam])[1] != 3):
                print 'Candidate does not seem to be in proper format. Skipping.'
                continue

            for cand in plotcands[beam]:
                # find island containing candidate
                for island in self.islands[beam]:
                    if n.any([n.all(member == cand) for member in island]):   # if any member of island that has all three coords the same as plotcand, break
                        break

                # next find snr for all candidates in island
                islandind = []
                for i in range(len(island)):
                    ww = n.where( (n.array(self.masterloc[beam])[:,0]==island[i,0]) & (n.array(self.masterloc[beam])[:,1]==island[i,1]) & (n.array(self.masterloc[beam])[:,2]==island[i,2]) )
                    islandind.append(ww[0][0])
                loc = n.array(self.masterloc[beam])[islandind]
                snr = n.array(self.masterprop[beam])[islandind, 0]

                # then plot snr as function of dmind and dtind for island
                fixed_dt = n.where(loc[:,0] == cand[0])
                fixed_dm = n.where(loc[:,2] == cand[2])
                dmdist = n.squeeze(n.array(self.dmarr)[loc[fixed_dt, 2]])   # seems to have superfluous axis...?
                dmsnr = snr[fixed_dt]
                dtdist = n.squeeze(n.array(self.timescales)[loc[fixed_dm][:, 0]])   # seems to have superfluous axis...?
                dtsnr = snr[fixed_dm]

                # find master index of candidate
                ww = n.where([n.all(mastercand == cand) for mastercand in self.masterloc[beam]])[0][0]

                # plot candidate info
                p.figure(figi)
                p.clf()

                p.subplot(221, axisbg='white')
                p.title('Candidate @ Tier 1')
                p.xticks([])
                p.yticks([])
                p.text(0.1, 0.8, self.datafile, fontname='sans-serif')
                beamra = n.round(beam, 3)[0]
                beamdec = n.round(beam, 3)[1]
                p.text(0.1, 0.6, 'Beam: ' + str(beamra) + ', ' + str(beamdec) + ', dt: ' + str(self.timescales[cand[0]]), fontname='sans-serif')
                p.text(0.1, 0.4, 'Integration: ' + str(cand[1]) + ', DM: ' + str(self.dmarr[cand[2]]), fontname='sans-serif')
                p.text(0.1, 0.2, 'SNR: ' + str(n.round(self.masterprop[beam][ww][0], 1)), fontname='sans-serif')

                p.subplot(222)
                p.plot(dmdist, dmsnr, 'b.', clip_on=False, label='DM at peak dt')
                p.xlabel('DM (pc/cm3)')
                p.ylabel('SNR')
                p.twiny()
                p.plot(dtdist, dtsnr, 'r+', clip_on=False, label='dt at peak DM')
                p.xlabel('dt (ints)')
                p.subplot(223)
                p.plot(self.masterdata[beam][ww][0], 'b.', label='Mean B')
                p.xlabel('Integration')
                p.ylabel('Mean, Std of Bispectra')
                p.twinx()
                p.plot(self.masterdata[beam][ww][1], 'r.', label='Std B')
                p.subplot(224)
                dataph = n.rot90(self.masterdata[beam][ww][2])
                sh = dataph.shape
                im = p.imshow(dataph, aspect='auto', origin='upper', interpolation='nearest', extent=(0, sh[1], 0, sh[0]))
                p.colorbar()
                p.plot(self.obs.dmtrack0[cand[2]][0], self.obs.dmtrack0[cand[2]][1],'k.')  # reference dispersed track
                p.xlabel('Integration')
                p.ylabel('Channel')
                if self.datatype == 'mir':   # casapy doesn't get along with recent matplotlib stuff
                    p.tight_layout()
                if save:
                    p.savefig(self.fileroot + '_sc' + str(self.scan) + 'sp' + str(self.spw[0]) + 'i' + str(self.startints[0]) + '_tier1_' + str(figi) + '.png')
                figi += 1

    def bispectrumplot(self, plotcand={}, show=1, save=0):
        """ Produces plot of bispectrum for given candidate.
        """

        if len(plotcand) != 1:
            print 'plotcands should have exactly one beam.'
            return 0

        # collect candidate params
        (ra, dec) = plotcand.keys()[0]
        cand = plotcand[(ra,dec)][0]
        startint = self.startints[0]

        print 'Producing bispectrum plot for beam', ra, dec, ', candidate, ', cand

        self.obs = tpipe.pipe_msdisp2(file=self.datafile, nints=self.nints, nskip=startint, profile='default', selectpol=self.selectpol, scan=self.scan, spw=self.spw, chans=self.chans, dmarr=self.dmarr, datacol='corrected_data')

        # select data for candidate
        self.obs.time_filter(cand[0], self.filtershape, bgwindow=self.bgwindow)
        noiseperbl = self.obs.data.mean(axis=3).mean(axis=2).real.std()   # measure single noise for input to detect_bispectra

        l1 = n.radians(ra); m1 = n.radians(dec)
        if ( (l1 != self.obs.l0[0]) | (m1 != self.obs.m0[0]) ):
            dl = l1 - self.obs.l0[0]; dm = m1 - self.obs.m0[0]  # find phase shift for new location
            self.obs.phaseshift(dl, dm)

        # find candidates
        self.obs.make_bispectra(stokes=self.stokes)
        (candsnr, candstd, cands) = self.obs.detect_bispectra(sigma=self.sigma, Q=noiseperbl, show=show, save=save)

    def candplot2(self, plotcands={}, save=0):
        """ Builds summary plots, including phased beam and imaging for each plot candidate.
        This is a "tier 2" plot, since it only uses all interferometric info for a given candidate data location.
        Candidate should be dictionary with keys of beam location and value of a list of 3-element lists (dt, int, dmbin).
        Note: this starts from scratch by reading data to produce a set of plots tailored to each candidate.
        No requirement for do_search and find_candidates to be run first.
        """

        print 'Building Tier 2 candidate plot.'

        if n.array([len(plotcands[beam]) for beam in plotcands.keys()]).sum() == 0:
            print 'No candidates available...'
            return

#        twindow = 128     # number of integrations to read for plots
        figi = 1
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

                if self.datatype == 'mir':
# older version centered pulse plot in time
#                    obs = tpipe.pipe_mirdisp2(file=self.datafile, nints=twindow, nskip=self.startints[0]+ii-twindow/2, profile='default', dmarr=self.dmarr, chans=self.chans)
                    self.obs = tpipe.pipe_mirdisp2(file=self.datafile, nints=self.nints, nskip=self.startints[0], profile='default', dmarr=self.dmarr, chans=self.chans)
                elif self.datatype == 'ms':
#                    obs = tpipe.pipe_msdisp2(file=self.datafile, nints=twindow, nskip=self.startints[0]+ii-twindow/2, profile='default', selectpol=self.selectpol, scan=self.scan, spw=self.spw, chans=self.chans, dmarr=self.dmarr, datacol='corrected_data')
                    self.obs = tpipe.pipe_msdisp2(file=self.datafile, nints=self.nints, nskip=self.startints[0], profile='default', selectpol=self.selectpol, scan=self.scan, spw=self.spw, chans=self.chans, dmarr=self.dmarr, datacol='corrected_data')
                elif self.datatype == 'sim':
                    print 'For datatype=\'sim\', we have no original data.'
#                    obs = tpipe.pipe_simdisp2(array='vla_d', nints=twindow, profile='default', dmarr=self.dmarr)
#                    obs.add_transient(0., 0., 1., twindow/2)
                    self.obs = tpipe.pipe_simdisp2(array='vla_d', nints=self.nints, profile='default', dmarr=self.dmarr)
                    self.obs.add_transient(0., 0., 1., ii)

                self.obs.time_filter(self.timescales[dtind], self.filtershape, bgwindow=self.bgwindow)    # tophat filter (bgwindow not used)
                noiseperbl = self.obs.data.mean(axis=3).mean(axis=2).real.std()   # measure single noise for input to detect_bispectra

                # shift phase center to candidate
                l1 = n.radians(ra); m1 = n.radians(dec)
                if ( (l1 != self.obs.l0[0]) | (m1 != self.obs.m0[0]) ):
                    dl = l1 - self.obs.l0[0]; dm = m1 - self.obs.m0[0]  # find phase shift for new location (in degrees)
                    self.obs.phaseshift(dl, dm)

                # reproduce bisp candidates
                self.obs.make_bispectra(stokes=self.stokes)
                (candsnr, candstd, cands) = self.obs.detect_bispectra(sigma=self.sigma, Q=noiseperbl, save=0)
                print 'Reproduced bispectrum candidate, ', cands

                # make image of field
#                im = obs.imagetrack(obs.tracksub(dmind, twindow/2), i=twindow/2, pol='i', size=self.size, res=self.res, clean=1, show=0)
                im = self.obs.imagetrack(self.obs.tracksub(dmind, ii), i=ii, pol='i', size=self.size, res=self.res, clean=1, show=0)
                imsize = float(self.size)/self.res
                pixelscale = 1./self.size
                peakm, peakl = n.where(im == im.max())
                dm = -(imsize/2. - peakm[0])*pixelscale   # indexed starting at top...
                dl = -(peakl[0] - imsize/2.)*pixelscale   # indexed starting at left.
                # phase shift to center for spectral analysis
                self.obs.phaseshift(dl, dm)
#                specmod = obs.specmod(dmind, twindow/2)
                specmod = self.obs.specmod(dmind, ii)

                # now plot
                p.figure(figi)

                # text description of candidate
                p.subplot(221, axisbg='white')
                p.title('Candidate @ Tier 2')
                p.xticks([])
                p.yticks([])
                p.text(0.1, 0.8, self.datafile, fontname='sans-serif')
                beamra = n.round(beam, 3)[0]
                beamdec = n.round(beam, 3)[1]
                p.text(0.1, 0.6, 'Beam: ' + str(beamra) + ', ' + str(beamdec) + ', dt: ' + str(self.timescales[cand[0]]), fontname='sans-serif')
                p.text(0.1, 0.4, 'Integration: ' + str(cand[1]) + ', DM: ' + str(self.dmarr[cand[2]]), fontname='sans-serif')
#                ww = n.where([n.all(mastercand == cand) for mastercand in self.masterloc[beam]])[0][0]
                p.text(0.1, 0.2, 'SNR: ' + str(n.round(candsnr, 1)), fontname='sans-serif')

                # image of dispersed visibilities for candidate
                p.subplot(222)
                fov = n.degrees(1./self.res)*3600.       # field of view in arcseconds
                p.imshow(im, aspect='auto', origin='lower', interpolation='nearest', extent=[fov/2, -fov/2, -fov/2, fov/2])
                p.colorbar()
                p.xlabel('RA/l Offset (arcsec)')
                p.ylabel('Dec/m Offset (arcsec)')

                p.subplot(223)
#                spec = obs.dedisperse(dmind).mean(axis=3).mean(axis=1).real[twindow/2]    # stokes I
                spec = self.obs.dedisperse(dmind).mean(axis=3).mean(axis=1).real[ii]    # stokes I
                p.plot(self.obs.freq, spec, '.')
                p.text(0.1, 0.1, 'specmod =' + str(specmod))
                p.xlabel('Frequency (GHz)')
                p.ylabel('Flux density (Jy)')

                p.subplot(224)
                dataph = n.rot90(self.obs.dataph)
                sh = dataph.shape
                im = p.imshow(dataph, aspect='auto', origin='upper', interpolation='nearest', extent=(0, sh[1], 0, sh[0]))
                p.colorbar()
                p.plot(self.obs.dmtrack0[dmind][0], self.obs.dmtrack0[dmind][1],'k.')   # reference dispersed track
                p.xlabel('Integration')
                p.ylabel('Channel')
                if self.datatype == 'mir':   # casapy doesn't get along with recent matplotlib stuff
                    p.tight_layout()

                if save:
                    p.savefig(self.fileroot + '_sc' + str(self.scan) + 'sp' + str(self.spw[0]) + 'i' + str(self.startints[0]) + '_tier1_' + str(figi) + '.png')
                figi += 1


################################
## high level functions ##
################################

def filtercands(cands):
    """ Function to make dictionary that only contains beams with candidates.
    """
    for beam in cands.keys():
        if len(cands[beam]) == 0:
            junk = cands.pop(beam)
            
    return cands
    
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
        totalcount = 0
        goodcands = []; descriptions = []
        shift = 0
        firsttime=1
        while 1:
            try:
                (description, cdict, noise, cands) = pickle.load(pkl)
                goodcands.append(filtercands(cands))
                descriptions.append(description)
                if firsttime:
                    descstr = description.split('\n')
                    print
                    print '**For setup of (last dump)**'
                    for desc in descstr:
                        print desc
                    firsttime=0
            except EOFError:
                break
        beams = []
        j = 0
        for i in range(len(goodcands)):
            canddict = goodcands[i]
            print descriptions[i].split('\n')[2]
            for beam in canddict:
                beams.append(beam)
                print 'Candidate ' + str(j),  canddict[beam], beam
                j += 1
                for cand in canddict[beam]:
                    p.text(beam[0]+shift, beam[1]+shift, str(cand[0]) + '-' + str(cand[1]) + '-' + str(cand[2]), horizontalalignment='center', verticalalignment='center', fontsize=9)
                    shift += 0.00003

        for beam in beams:
            p.plot([beam[0], beam[0]+shift], [beam[1], beam[1]+shift], 'k-')

    p.show()

def reproduce_cands(name, candnum=-1, plot=0, tier1=0, tier2=1, bisp=0, getpipe=0, save=0):
    """ Read pickle file to identify and/or plot candidates.
    Summarizes parameters of candidates. Optionally reproduces plots for candidates for saving.
    """

    print
    print 'Looking at candidates in %s...' % (name)

    with open(name, 'r') as pkl:
        totalcount = 0
        itot = 0
        goodcands = []; cdicts = []
        firsttime=1
        while 1:
            try:
                (description, cdict, noise, cands) = pickle.load(pkl)
                descstr = description.split('\n')
                if firsttime:
                    print
                    print '**For setup of**'
                    for desc in descstr:
                        print desc
                    firsttime = 0
                else:
                    print descstr[2]
                goodcands.append(filtercands(cands))
                cdicts.append(cdict)
            except EOFError:
                break

        j = 0
        for i in range(len(goodcands)):
            canddict = goodcands[i]
            cdict = cdicts[i]
            for beam in canddict.keys():
                print 'Found %d candidates in this line from pickle file.' % (len(canddict.keys()))
                if (plot and ((j == candnum) or (candnum == -1))):    # plot if asked and candnum is this iteration or default value
                    sp = search(datatype=cdict['datatype'], datafile=cdict['datafile'], startints=cdict['startints'], nints=cdict['nints'], timescales=cdict['timescales'], filtershape=cdict['filtershape'], dmarr=cdict['dmarr'], fwhmsurvey=cdict['fwhmsurvey'], fwhmfield=cdict['fwhmfield'], selectpol=cdict['selectpol'], stokes=cdict['stokes'], scan=cdict['scan'], spw=cdict['spw'], chans=cdict['chans'], bgwindow=cdict['bgwindow'], sigma=cdict['sigma'], size=cdict['size'], res=cdict['res'])
                            
                    sp.delaycenters = [beam]     # only search in beams with candidates
                    print sp.description()

                    if tier1:     # optionally can replot tier1 figure
                        print
                        print 'Searching over delaycenters:', sp.delaycenters
                        sp.do_search()
                        plotcands = sp.find_candidates(show=0, island_snr_min=sp.sigma)
                        print 'Candidates per delay beam:', [(beam, len(plotcands[beam])) for beam in plotcands.keys()]
                        sp.candplot1(plotcands=plotcands, save=save)
                    if tier2:
                        sp.candplot2(plotcands={beam: canddict[beam]}, save=save)
                    if bisp:
                        sp.bispectrumplot(plotcand={beam: canddict[beam]}, save=save, show=1)

                j += 1
        print 'Total cand count in pickle: ', j
        if (plot and getpipe):
            return sp

def mbsearch():
    """ Defines search parameters for large scale search (over multiple scans, times, beams, dms, timescales).
    All results, including noise levels, saved in pickle.
    """

    print 'Starting multibeam bispectrum search...'
    import time, sys

# data
    datatype = 'ms'
#    datafile = '12A-336_1b_j0628.ms'
    datafile = '12A-336_J1911_s4~7.ms'
#    datafile = '12A-336_J1851.ms'
# search domain
#    startints = range(10000,18000,110)    # dm=150 delay span is 15 ints for vla rrat data plus 3 for one side of time filter
#    startints = range(2000,12000,110)
    startints = range(2000,5000,238)
    nints = 256    # number of ints to read at a time
    timescales = [1]    # number of ints in filter
    filtershape = None   # 't', None, 'b'
    dmarr = range(0,175,33)
#    dmarr = range(0,150,22)
#    fwhmsurvey = 0.233   # field of view to tile (fwhm)  # 0.233 is actual PMPS fwhm
    fwhmsurvey = 0.03   # field of view to tile (fwhm)  # 0.233 is actual PMPS fwhm
    fwhmfield = 0.008   # delay beam size (fwhm)
#    fwhmsurvey = 0.5   # field of view to tile (fwhm)  # 0.233 is actual PMPS fwhm
#    fwhmfield = 0.5   # delay beam size (fwhm)
    selectpol = ['I']
#    selectpol = ['RR','LL']
    stokes = 'prebisp'
    scans = range(1)
    spw = [0]
    chans = range(5,50) + range(51,59)
#    chans = range(5,24) + range(25,50) + range(51,59)
#    chans = range(64)
# search params
    bgwindow = 4
    sigma = 5
# imaging
    size = 350000   # vla-a
    res = 500  # vla

#    log = open('search.txt', 'a')    # redefine stdout for automatic search logging
#    temp = sys.stdout
#    sys.stdout = log
    tt = time.localtime()
    print 'Start time: %s_%s_%s:%s:%s:%s' % (tt.tm_year, tt.tm_mon, tt.tm_mday, tt.tm_hour, tt.tm_min, tt.tm_sec)
    timestring = '%s_%s_%s:%s:%s:%s' % (tt.tm_year, tt.tm_mon, tt.tm_mday, tt.tm_hour, tt.tm_min, tt.tm_sec)
    pkl = open('mbsearch_'+timestring+'.pkl', 'w')

    for scan in scans:
        for startint in startints:
            sp = search(datatype=datatype, datafile=datafile, startints=[startint], nints=nints, timescales=timescales, filtershape=filtershape, dmarr=dmarr, fwhmsurvey=fwhmsurvey, fwhmfield=fwhmfield, selectpol=selectpol, stokes=stokes, scan=scan, spw=spw, chans=chans, bgwindow=bgwindow, sigma=sigma, size=size, res=res)
            description = sp.description()
            print description
#            pipegenstring = 'sp = bisp_mbsearch.search(datatype=\'%s\', datafile=\'%s\', startints=[%s], nints=%s, timescales=%s, filtershape=\'%s\', dmarr=%s, fwhmsurvey=%s, fwhmfield=%s, selectpol=%s, stokes=\'%s\', scan=%s, spw=%s, chans=%s, bgwindow=%s, sigma=%s, size=%s, res=%s)' % (datatype, datafile, startint, nints, timescales, filtershape, dmarr, fwhmsurvey, fwhmfield, selectpol, stokes, scan, spw, chans, bgwindow, sigma, size, res)
            constructiondict = {'datatype':datatype, 'datafile':datafile, 'startints':[startint], 'nints':nints, 'timescales':timescales, 'filtershape':filtershape, 'dmarr':dmarr, 'fwhmsurvey':fwhmsurvey, 'fwhmfield':fwhmfield, 'selectpol':selectpol, 'stokes':stokes, 'scan':scan, 'spw':spw, 'chans':chans, 'bgwindow':bgwindow, 'sigma':sigma, 'size':size, 'res':res}   # dictionary used to store pipe construction parameters
            sp.do_search()
            plotcands = sp.find_candidates(show=0, island_snr_min=sp.sigma)
            print 'Candidates per delay beam:', [(beam, len(plotcands[beam])) for beam in plotcands.keys()]
            sp.candplot1(plotcands, save=1)
            pickle.dump((description, constructiondict, sp.masternoise, plotcands), pkl)

    tt = time.localtime()
    print 'Stop time: %s_%s_%s:%s:%s:%s' % (tt.tm_year, tt.tm_mon, tt.tm_mday, tt.tm_hour, tt.tm_min, tt.tm_sec)
    pkl.close()
#    log.close()
#    sys.stdout = temp  # restore stdout

