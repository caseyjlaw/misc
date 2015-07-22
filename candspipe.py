# Functions process candidate pkl files from leanpipedt
# Assumes structure of FRB search: multiple pkls per scan/ms file
# Batch output writes candidtes per chunk of ints, but after merging, there is only one entry for all cands.
# Also, a "noise" file is named "noise" instead of "cands" and contains one entry per chunk of ints.

import numpy as n
#import pylab as p    # requires display before saving fig...
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
import multiprocessing as mp
import string, ctypes, glob, types, os
import cPickle as pickle
import leanpipedt as leanpipedt
import leanpipedt_cython as lib
import scipy.special as sp

def mergecands(pkls, outname=''):
    """ Converts output of batch processing into single pkl file for visualization below.
    Takes string to select list or list of pkl file names and merges them into a single pkl file with two entries instead of one per scan with cands.
    Assumes all pkl files are part of same basic dmarr.
    In the future, could also merge across scans.
    """

    if isinstance(pkls, types.ListType):
        pkllist = pkls
    elif isinstance(pkls, types.StringType):
        pkllist = glob.glob(pkls)

    pkllist = [pkllist[i] for i in range(len(pkllist)) if 'merge' not in pkllist[i]]
    pkllist.sort(key=lambda i: int(i.rstrip('.pkl').split('_')[3][1:]))  # assumes filename structure

    if len(pkllist) == 0:
        return

    if not outname:
        outname = pkllist[0].split('_s')[0] + '_merge.pkl'

    if os.path.exists(outname):
        print 'Outname %s exists. Not overwriting.' % outname
        return

    print 'Iterating over list %s. Writing to %s' % (pkllist, outname)

    for pklfile in pkllist:
        if 'merge' in pklfile:
            return

    # collect construction dicts and merge chunks in time
    dd = {}
    masterdmarr = []; mastercandloc = []; mastercandprop = []
    for i in range(len(pkllist)):
        pkl = open(pkllist[i], 'r')
        dd[i] = pickle.load(pkl)
        candloc = []; candprop = []
        for cands in pickle_iter(pkl):
            loc, prop = cands
            candloc = candloc+loc
            candprop = candprop+prop
        pkl.close()
        print 'Merged file %s over %d/%d (loc/prop) chunks.' % (pkllist[i], len(candloc), len(candprop))
        masterdmarr = masterdmarr + list(dd[i]['dmarr'])
        dd['loc'+str(i)] = candloc
        dd['prop'+str(i)] = candprop

    # regrid dmarr
    masterdmarr.sort()
    masterdmarr = list(n.unique(masterdmarr))
    masterdataloc = []
    for i in range(len(pkllist)):
        dmbin0 = masterdmarr.index(dd[i]['dmarr'][0])
        scani = int(pkllist[i].rstrip('.pkl').split('_')[3][1:])   # assumes scan name structure
        masterdataloc = masterdataloc + [scani]*len(dd['loc'+str(i)])
        mastercandloc = mastercandloc + [(cand[0], cand[1], cand[2], cand[3]+dmbin0) for cand in dd['loc'+str(i)]]
        mastercandprop = mastercandprop + dd['prop'+str(i)]
#    print 'Master dmarr grid', masterdmarr

    # define master dictionary
    masterd = dd[i]
    masterd['pkllist'] = pkllist
    masterd['dmarr'] = n.array(masterdmarr)
    outputpkl = open(outname, 'wb')
    pickle.dump(masterd, outputpkl)
    pickle.dump((masterdataloc,mastercandloc,mastercandprop), outputpkl)   # merge file has "dataloc" an index to d['pkllist'] to say which file the cand was in

def mergecands_plotall(rootname, saveroot=''):
    """ Wrapper for mergecands that takes argument 'cands*pkl' and makes all merge files, less any extant merge files.
    Then generates summary plots (image, dmt) for all candidates.
    """

    pkllist = glob.glob(rootname)
    scanlist = n.unique([pkllist[i].rstrip('.pkl').split('_')[3] for i in range(len(pkllist))])

    for scan in scanlist:
        sublist = [pklfile for pklfile in glob.glob('cands*_'+scan+'_*pkl') if 'merge' not in pklfile]   # filter out merge files
        mergecands(sublist)

    if saveroot:
        # make summary plots
        cands = filter_candidates('cands*merge*pkl', saveroot=saveroot)    # summary plot for all
        if len(cands) == 2:
            candsloc, candsprop = cands
        elif len(cands) == 3:
            dataloc, candsloc, candsprop = cands

        if len(candsprop[0]) > 1:   # for long form of cand info
            compresscands(glob.glob('cands*merge*pkl'), rootname.split('*')[0]+'_summary.pkl')
        # plot all candidates
            for pklfile in glob.glob('cands*merge*pkl'):
                candsloc, candsprop = filter_candidates(pklfile, saveroot='', threshold=7.)    # only plot cands for snr>7
                for i in range(len(candsloc)):
                    singleplot(pklfile, candsloc, candsprop, i, save=1)

def mergeremove(mergefile, outfile, remove={}):
    """ Function writes new merge file without specified ranges of scan/ints (containing RFI)
    remove is a dict to define data ranges to excise from candidate set: remove = {scan0: [int0,int1,int2,int3]}
    where int0-int1 defines first range in that scan and int2-int3 is next, etc...
    """

    pkl = open(mergefile,'r')
    d = pickle.load(pkl)
    cands = pickle.load(pkl)

    if len(remove) == 0:
        print 'No times to remove.'
    else:
        print 'Removing %s' % remove

    dataloc, candloc, candprop = cands
    dataloc2 = []; candloc2 = []; candprop2 = []

    for i in xrange(len(dataloc)):
        if dataloc[i] in remove.keys():
            nranges = len(remove[dataloc[i]])
            escape = 0
            for first in range(0,nranges,2):
                badrange0 = remove[dataloc[i]][first]
                badrange1 = remove[dataloc[i]][first+1]
                if ((candloc[i][2] > badrange0) & (candloc[i][2] < badrange1)):
                    escape = 1
            if escape:
                continue
        dataloc2.append(dataloc[i])
        candloc2.append(candloc[i])
        candprop2.append(candprop[i])

    # if the scan has any candidates, add nints to count
    goodintcount = 0
    for scani in n.unique(dataloc):
        if len(n.where(dataloc == scani)[0]):   
            goodintcount += d['nints']

    # correct int count by removed range
    for scan in remove.keys():
        goodintcount -= remove[scan][1] - remove[scan][0]

    # write filtered output
    d['remove'] = remove
    d['goodintcount'] = goodintcount
    outputpkl = open(outfile, 'wb')
    pickle.dump(d, outputpkl)
    pickle.dump((dataloc2,candloc2,candprop2), outputpkl)

def pickle_iter(infile):
    while 1:
        try:
            yield pickle.load(infile)
        except EOFError:
            break

def filter_candidates(pklroot, threshold=-999, saveroot=''):
    """ Quick ingestion and filtering of candidate pkl file.
    Returns candidate location and property arrays.
    """

    locs = []; snrs = []; props = []
    # read in pickle file of candidates
    pkllist = glob.glob(pklroot)
    for pklfile in pkllist:
        print 'Filling from first data in %s...' % pklfile
        pkl = open(pklfile, 'rb')
        d = pickle.load(pkl)
        cands = pickle.load(pkl)
        if len(cands) == 2:
            loc, prop = cands
        elif len(cands) == 3:
            dataloc, loc, prop = cands

        snrcol = parsepropinfo(prop)

        pkl.close()
        locs = locs + loc
        snrs = snrs + [prop[i][snrcol] for i in range(len(prop))]
        props = props + [prop[i] for i in range(len(prop))]

    candsloc = n.array(locs).astype(int)
    snrs = n.array(snrs)
    ww = n.where(snrs>threshold)[0].astype(int)
    candsloc = candsloc[ww]
    candsprop = [props[i] for i in ww]

    print 'Found %d candidates above threshold %.1f' % (len(candsloc), threshold)

    if saveroot and len(candsprop):
        dmtplot(d, candsloc, snrs, saveroot=saveroot)
        if len(candsprop[0]) > 1:   # then long form of cand info
            impeakplot(d, candsprop, saveroot=saveroot)
        snrhist(snrs, saveroot=saveroot)
    if len(cands) == 2:
        return candsloc, candsprop
    elif len(cands) == 3:
        return n.array(dataloc)[ww], candsloc, candsprop

def parsepropinfo(prop):
    """ Function to find columns of interest
    """

    if len(prop[0]) == 6:
        snr0col = 0; x0col = 1; x0col = 2; snrcol = 3; xcol = 4; ycol = 5
    elif len(prop[0]) == 3:
        snrcol = 0; imcol = 1; speccol = 2
    return snrcol

def quickcands(pklname):
    pkl = open(pklname,'r')
    candloc = []; snr = []

    d = pickle.load(pkl)
    for cands in pickle_iter(pkl):
        if len(cands) == 2:
            loc, prop = cands
        elif len(cands) == 3:
            dataloc, loc, prop = cands
        candloc = candloc + loc
        snr = snr + [prop[i][0] for i in range(len(prop))]

    return n.array(candloc), n.array(snr)

def compresscands(pkllist, outname):
#    pkllist = glob.glob(pklroot)
    loclist = []; snrlist = n.array([]); fnlist = n.array([])
    for pklname in pkllist:
        print 'Processing %s (%d/%d)' % (pklname, pkllist.index(pklname)+1, len(pkllist))
# using quickcands doesn't give im peak pixel
#        loc,snr = quickcands(pklname)
#        for ll in loc:
#            loclist.append(ll)
# using filter_candidates allows full "loc" definition
        cands = filter_candidates(pklname)
        if len(cands) == 2:
            candsloc, candsprop = cands
        elif len(cands) == 3:
            dataloc, candsloc, candsprop = cands

        for i in range(len(loc)):
            peakx, peaky = n.where(prop[i][1] == prop[i][1].max())
            loclist.append(n.concatenate( (loc[i], peakx, peaky) ))

        snr = n.array([prop[i][0] for i in range(len(prop))])
        snrlist = n.append(snrlist, snr)

        fn = n.array([pklname] * len(loc))
        fnlist = n.append(fnlist, fn)

    pkl = open(outname, 'w')
    pickle.dump((fnlist, n.array(loclist), snrlist), pkl)
    pkl.close()

def reproducecand(pklname, threshold=-1, candnum=-1, cand = [], **kwargs):
    """ Takes merged pickle file and reproduces candidate candnum from list above thresh.
    Optionally takes candidate params as list [dt, int, dm] (dt/dm are values not indexes)
    kwargs used to modify parameters of search.
    """

    # get run metadata and basic candidate info
    pkl = open(pklname, 'r')
    d = pickle.load(pkl)
    pkl.close()
    d['savecands'] = False

    # overwrite d with kwargs. supports values of type int and string so far.
    for key in kwargs.keys():
        print 'Overwriting key %s to %s' % (key, kwargs[key])
        d[key] = kwargs[key]
    if 'excludeants' not in d.keys():
        d['excludeants'] = []

    window = 200
    maxiter = 400
    # define time and iteration parameters to get candidate
    if (len(cand) == 0) and (candnum != -1):  # if not forced, use cand from pkl file
        dataloc, locs, props = filter_candidates(pklname, threshold)
        snrcol = parsepropinfo(props)
        canddmind = locs[candnum][3]
        canddtind = locs[candnum][1]
        dtarr = [d['dtarr'][canddtind]]
        dmarr = [d['dmarr'][canddmind]]
        nints = d['datadelay'][canddmind] + window
        if 'cand' in d['filename']:
            filename = d['filename']
            nskip = 0
        else:
            scani = [int(name.rstrip('.pkl').split('_')[3][1:]) for name in d['pkllist']]
            filesplit = d['pkllist'][scani.index(dataloc[candnum])].rstrip('.pkl').split('_')
            filename = '_'.join(filesplit[1:4]) + '.ms'
            nskip = locs[candnum][2] - (d['datadelay'][canddmind] + window/2)
# method 1
#        iterint = nints/(((nints-n.mod(nints,maxiter))/maxiter)+2)    # keeps iterint under maxiter
# method 2
        iterint = nints-1
        print 'Reproducing candidate %d with SNR = %.1f at (%d,%d,%d) in file %s' % (candnum, props[candnum][snrcol], canddtind, nskip, canddmind, filename)
    elif len(cand) == 3:
        print 'Reproducing provided candidate with values %s' % (cand)
        filename = d['filename']
        candint = cand[1]
        canddm = cand[2]
        canddtind = d['dtarr'].index(cand[0])
        canddmind = min(range(len(d['dmarr'])), key = lambda i: abs(d['dmarr'][i] - canddm))
        dtarr = [d['dtarr'][canddtind]]
        dmarr = [d['dmarr'][canddmind]]
        nskip = candint - (d['datadelay'][canddmind] + window/2)
        nints = d['datadelay'][canddmind] + window
# method 1
#        iterint = nints/(((nints-n.mod(nints,maxiter))/maxiter)+2)    # keeps iterint under maxiter
# method 2
        iterint = nints-1
    elif (len(cand) == 0) and (('nskip' in kwargs.keys()) or ('nints' in kwargs.keys()) or ('iterint' in kwargs.keys())):
        filename = d['filename']
        print 'Reading data at given location for file %s' % d['filename']
        nskip = d['nskip']
        nints = d['nints']
        iterint = d['iterint']
        dtarr = d['dtarr']
        dmarr = d['dmarr']
    else:
        print 'Not given proper cand or data definition.'
        return

    print 'Running new search with %d ints, %d iterint, starting at %d.' % (nints, iterint, nskip)

    dd = leanpipedt.pipe_thread(filename=filename, nints=nints, iterint=iterint, nskip=nskip, spw=d['spw'], chans=d['chans'], scan=d['scan'], dmarr=dmarr, dtarr=dtarr, selectpol=d['selectpol'], datacol=d['datacol'], size=d['size'], res=d['res'], sigma_image=d['sigma_image'], searchtype=d['searchtype'], gainfile=d['gainfile'], bpfile=d['bpfile'], filtershape=d['filtershape'], secondaryfilter=d['secondaryfilter'], savecands=d['savecands'], flagmode=d['flagmode'], excludeants=d['excludeants'], nthreads=1)
    return dd

def reproducecanddata(d, candint=-1, dmind=0, dtind=0, twindow=30, **kwargs):
    """ Runs imaging as in initial detection. Needs to have dictionary and leanpipedt.data populated by run of reproducecand()
    Default images all with 'sizex/y'. Given integration, it will image once with 'full_sizex/y'.
    Assumes that entire dm track contained in data and datatrim is not needed.
    """

    import qimg_cython as qimg

    for key in kwargs.keys():
        print 'Overwriting key %s to %s' % (key, kwargs[key])
        d[key] = kwargs[key]
    if 'excludeants' not in d.keys():
        d['excludeants'] = []

    data0 = leanpipedt.dataprep(d, dmind, dtind, usetrim=False)[d['datadelay'][dmind]/d['dtarr'][dtind]:]    # returns masked array of dedispersed, ignores times with self-added data

    if candint == -1:
        print 'Imaging all integrations and thresholding.'
        ims,snr,candints = qimg.imgallfullfilterxy(n.outer(leanpipedt.u[d['iterint']/2], d['freq']/d['freq_orig'][0]), n.outer(leanpipedt.v[d['iterint']/2], d['freq']/d['freq_orig'][0]), data0.data, d['sizex'], d['sizey'], d['res'], d['sigma_image'])
        return ims,snr,candints
    else:
        print 'Imaging int %d.' % candint
        im = qimg.imgonefullxy(n.outer(leanpipedt.u[candint], d['freq']/d['freq_orig'][0]), n.outer(leanpipedt.v[candint], d['freq']/d['freq_orig'][0]), data0.data[candint], d['full_sizex'], d['full_sizey'], d['res'])
        # measure spectrum
        print 'Full image peak SNR=', im.max()/im.std()
        peakl, peakm = n.where(im == im.max())      # assumes new style u->u and v->v gridding
        l1 = (float((d['full_sizex'])/d['res'])/2. - peakl[0])/d['full_sizex']
        m1 = (float((d['full_sizey'])/d['res'])/2. - peakm[0])/d['full_sizey']
        # trim interesting int out
        minint = max(candint-twindow/2, 0)
        maxint = min(candint+twindow/2, len(data0))
        data0 = data0[minint:maxint].copy()
        lib.phaseshift_threaded(data0, d, l1, m1, leanpipedt.u[candint], leanpipedt.v[candint])

    if candint == -1:
        return ims,snr,candints
    else:
        return im, data0.mean(axis=1)

def extractms(pklname, threshold, candnum=-1):
    """ Takes summary pickle file and produces MS with snippet of data for each candidate with snr>threshold.
    Default is to extract for all cands above threshold. Otherwise take candnum (0-based count).
    """

    try:
        from casa import ms, split
        from casa import quanta as qa
        print 'Imported CASA'
        incasapy = True
    except ImportError:
        import tasklib as tl
        import casautil
        print 'Imported CASA sans casapy'
        ms = casautil.tools.ms()
        qa = casautil.tools.quanta()
        incasapy = False

    secondsinday = 24.*60*60
    daypersecond = 1.0/secondsinday

    # get run metadata and basic candidate info
    d = pickle.load(open(pklname, 'r'))
    dataloc, locs, props = filter_candidates(pklname, threshold)
    scani = [int(name.rstrip('.pkl').split('_')[3][1:]) for name in d['pkllist']]
    filesplit = d['pkllist'][scani.index(dataloc[candnum])].split('_')
    filename = '_'.join(filesplit[1:4]) + '.ms'
    snrcol = parsepropinfo(props)

    canddt = locs[candnum][1]
    candint = locs[candnum][2]
    canddm = locs[candnum][3]
    window = d['datadelay'][canddm]    # minimally clip data
    outms = filename[:-3] + '_cand' + str(canddt) + '-' + str(candint) + '-' + str(canddm) + '.ms'
    if not os.path.exists(outms):
        print 'Extracting candidate %d with SNR = %.1f at bins %s in file %s into %s' % (candnum, props[candnum][snrcol], locs[candnum], filename, outms)
        ms.open(filename)
        t0 = ms.range(items=['time'])['time'][0]/(24*3600)
        dt = d['inttime']/secondsinday   # int time in days
        ms.done()

        minint = max(0, candint - window - 100)
        maxint = candint + 100     # 100 is largest iterint in overall search so it represents high end of candint location error. include extra 15 for plotting window and fft time filter protection(?)
        candt0 = t0 + minint*dt
        candt1 = t0 + maxint*dt   
        candt0q = qa.time(qa.quantity(candt0, 'd'), form='ymd', prec=9)[0]
        candt1q = qa.time(qa.quantity(candt1, 'd'), form='ymd', prec=9)[0]

        print 'Splitting out int %d -- %d and times %s -- %s' % (minint, maxint, candt0q, candt1q)

        if incasapy:  # if in casapy
            split(vis=filename, outputvis=outms, datacolumn='data', timerange=candt0q+'~'+candt1q)
        else: # if in python
            cfg = tl.SplitConfig()  # configure split
            cfg.vis = filename
            cfg.out = outms
            cfg.col = 'data'
            cfg.timerange=candt0q+'~'+candt1q
            tl.split(cfg)  # run task
    else:
        print 'File %s exists. Not overwriting.' % outms
    return outms

def mergeplot(mergepkl, outroot=''):
    """ Take merge file and produces comprehensive candidate screening plots.
    Starts as dm-t plots, includes dt and peak pixel location.
    """

    if not outroot:
        outname = mergepkl[:-4] + '_dmt.png'
        outname2 = mergepkl[:-4] + '_dmcount.png'
        outname3 = mergepkl[:-4] + '_normprob.png'
        outname4 = mergepkl[:-4] + '_impeak.png'
        outname5 = mergepkl[:-4] + '_noise.png'
    else:
        outname = outroot + '_dmt.png'
        outname2 = outroot + '_dmcount.png'
        outname3 = outroot + '_normprob.png'
        outname4 = outroot + '_impeak.png'
        outname5 = outroot + '_noise.png'

    print 'Plotting to %s, %s, %s, %s, %s.' % (outname, outname2, outname3, outname4, outname5)

    pkl = open(mergepkl, 'r')
    d = pickle.load(pkl)
    stuff = pickle.load(pkl)
    if len(stuff) == 3:    # new format for merged pkl file
        dataloc = n.array(stuff[0])
        candsloc = n.array(stuff[1])
        candsprop = n.array(stuff[2])

    if 'remove' in d.keys():
        remove = d['remove']
    else:
        remove = {}

    # measure scan change points for labels
    changepoints = [dataloc[0]]
    for i in range(len(dataloc)-1):
        if (dataloc[i+1] - dataloc[i]):
            changepoints.append(dataloc[i+1])
    changepoints = n.array(changepoints)
    shifts = dataloc * d['nints']

    dts = n.unique(candsloc[:,1])
    tt = d['inttime']*(candsloc[:,2]+shifts)
    mint = tt.min(); maxt = tt.max()
    dd = n.array(d['dmarr'])[candsloc[:,3]]
    mindm = dd.min(); maxdm = dd.max()
    snrs = candsprop[:,3]
    snrmin = 0.9*snrs.min()

#    fig = plt.Figure(figsize=(8,8))
    fig = plt.Figure(figsize=(15,10))
    ax = {}
    for dtind in range(len(dts)):
#    for dtind in [0]:
        # dmt plot
        ax[dtind] = fig.add_subplot(str(len(dts)) + '1' + str(dtind+1))
#        ax[dtind] = fig.add_subplot('111')
        good = n.where(candsloc[:,1] == dtind)[0]
        times = d['inttime']*(candsloc[good,2]+shifts[good])
        dms = n.array(d['dmarr'])[candsloc[good,3]]
#        sizes = (snrs[good]-0.8*snrs[good].min())
        sizes = (snrs[good]-snrmin)**5   # set scaling to give nice visual sense of SNR
        ax[dtind].scatter(times, dms, s=sizes, facecolor='none', alpha=0.3, clip_on=False)
        ax[dtind].axis( (mint, maxt, mindm, maxdm) )
        ax[dtind].set_ylabel('DM (pc/cm3)')
        ax[dtind].text(0.9*maxt, 0.9*maxdm, 'dt='+str(d['dtarr'][dtind]))
        if dtind == dts[-1]:
            plt.setp(ax[dtind].get_xticklabels(), visible=True)
        elif dtind == dts[0]:
            ax[dtind].xaxis.set_label_position('top')
            ax[dtind].xaxis.set_ticks_position('top')
            ax[dtind].set_xticks(changepoints[::2]*d['inttime']*d['nints'])
            ax[dtind].set_xticklabels(changepoints[::2])
            plt.setp( ax[dtind].xaxis.get_majorticklabels(), rotation=90)
        else:
            plt.setp( ax[dtind].get_xticklabels(), visible=False)

    ax[dtind].set_xlabel('Time (s)', fontsize=20)
#    ax[dtind].set_xlabel('Scan number', fontsize=20)
    ax[dtind].set_ylabel('DM (pc/cm3)', fontsize=20) 
    canvas = FigureCanvasAgg(fig)
    canvas.print_figure(outname)

    # dmcount plot
    fig2 = plt.Figure(figsize=(15,10))
    ax2 = {}
    for dtind in range(len(dts)):
        good = n.where(candsloc[:,1] == dtind)[0]
        times = d['inttime']*(candsloc[good,2]+shifts[good])
        times = n.concatenate( (n.array([0]), times))
        bins = n.round(times).astype('int')
        counts = n.bincount(bins)

        ax2[dtind] = fig2.add_subplot(str(len(dts)) + '1' + str(dtind+1))
        ax2[dtind].scatter(range(len(counts)), counts, facecolor='none', alpha=0.5, clip_on=False)
        ax2[dtind].axis( (mint, maxt, 0, 1.1*counts.max()) )

        # label high points
        high = n.where(counts > n.median(counts) + 20*counts.std())[0]
        for ii in high:
            print '%d candidates for dt=%d at %d s' % (counts[ii], d['dtarr'][dtind], ii)
            ww = n.where(bins == ii)[0]
            print '\tFlag these:', zip(dataloc[good][ww], candsloc[good,2][ww])
            print

        if dtind == dts[-1]:
            plt.setp(ax2[dtind].get_xticklabels(), visible=True)
        elif (dtind == dts[0]) or (dtind == len(dts)/2):
            ax2[dtind].xaxis.set_label_position('top')
            ax2[dtind].xaxis.set_ticks_position('top')
            ax2[dtind].set_xticks(changepoints[::2]*d['nints']*d['inttime'])
            ax2[dtind].set_xticklabels(changepoints[::2])
            plt.setp( ax2[dtind].xaxis.get_majorticklabels(), rotation=90, size='small')
        else:
            plt.setp( ax2[dtind].get_xticklabels(), visible=False)

    ax2[dtind].set_xlabel('Time (s)')
    ax2[dtind].set_ylabel('Count') 
    canvas2 = FigureCanvasAgg(fig2)
    canvas2.print_figure(outname2)

    # DM SNR norm prob plot
    if 'goodintcount' in d.keys():
        Z = lambda quan: n.sqrt(2)*sp.erfinv( 2*quan - 1) 
        quan = lambda ntrials, i: (ntrials + 1/2. - i)/ntrials
        snrsort = n.array(sorted(snrs, reverse=True))     # high-res snr
#        snrsort = n.array(sorted(candsprop[:,0], reverse=True))    # low-res snr
        npix = (d['full_sizex']/d['res']) * (d['full_sizey']/d['res'])
        nints = d['goodintcount']
        ndms = len(d['dmarr'])
        ntrials = npix*nints*ndms*(1 + 1/2. + 1/4. + 1/8.)
        Zsort = n.array([Z(quan(ntrials, j+1)) for j in range(len(snrsort))])

        fig3 = plt.Figure(figsize=(10,10))
        ax3 = fig3.add_subplot(111)
        ax3.plot(snrsort, Zsort, 'r.')
        refl = n.linspace(min(snrsort.min(), Zsort.min()), max(snrsort.max(), Zsort.max()), 2)
        ax3.plot(refl, refl, 'k--')
        ax3.set_xlabel('SNR')
        ax3.set_ylabel('Normal quantile SNR')
        canvas = FigureCanvasAgg(fig3)
        canvas.print_figure(outname3)

    fig4 = plt.Figure(figsize=(10,10))
    ax4 = fig4.add_subplot(111)
    sizes = (snrs-0.9*snrs.min())**5   # set scaling to give nice visual sense of SNR
    xarr = 60*n.degrees(candsprop[:,4]); yarr = 60*n.degrees(candsprop[:,5])
    ax4.scatter(xarr, yarr, s=sizes, facecolor='none', alpha=0.5, clip_on=False)
    ax4.set_xlabel('Dec Offset (amin)')
    ax4.set_ylabel('RA Offset (amin)')
    fov = n.degrees(1./d['res'])*60.
    ax4.set_xlim(fov/2, -fov/2)
    ax4.set_ylim(-fov/2, fov/2)
    canvas4 = FigureCanvasAgg(fig4)
    canvas4.print_figure(outname4)

    noisehists('noise*pkl', outname5, remove=remove)

def calc_ntrials(candsdir):
    """ Function to read config dict from start of cand pickle file and calculate the total number of trials (nint * ndm/nresample * npix)
    Expects directory in which to look for cands*merge*pkl files.
    Returns total number of trials over all cand files.
    """

    pkllist = glob.glob(candsdir + 'cands*merge*pkl')
    ntot = 0
    print 'Counting over %d cands in directory %s...' % (len(pkllist), candsdir)
    for pkl in pkllist:
        d = pickle.load(open(pkl, 'r'))
        if 'resamplearr' in d.keys():     # some resamplearrs have only 1s, due to bug in merge step. if key exists, it should be filled with this.
            resamplearr = n.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2])
        else:     # for the few v6 processing files, resamplearr is effectively all 1s.
            resamplearr = n.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
        n_int_dm = (d['nints']*n.ones(len(d['dmarr']), dtype='int64')/resamplearr).sum()
        n_pix = d['sizex']/d['res'] * d['sizey']/d['res']
        n_trials = n_int_dm * n_pix
        print 'For file %s: %d trials' % (pkl, n_trials)
        ntot += n_trials
        
    print 'For sig=6.5, tail prob=4.016e-11 and false pos rate = ', 4.016e-11 * ntot
    return ntot

def impeakplot(d, candsprop, saveroot=''):
    """ Take images of all candidates and shows peak pixel location relative to each other on summed map.
    """

    imlist = [candsprop[i][1] for i in range(len(candsprop))]
    lenx = 0; leny = 0
    for im in imlist:
        sh = im.shape
        if sh[0] > leny:
            leny = sh[0]
        if sh[1] > lenx:
            lenx = sh[1]
    imsum = n.zeros((leny, lenx), dtype=n.dtype(candsprop[0][1][0,0]))
    impeak = []
    for im in imlist:
        sh = im.shape
        imsum[:sh[0],:sh[1]] = imsum[:sh[0],:sh[1]] + im
        peaky, peakx = n.where(im == im.max())
        impeak = impeak + [[peaky[0],peakx[0]]]
            
    impeak = n.array(impeak)
    if saveroot:
        savefile = saveroot + '_imsum.png'
        print 'Plotting to %s.' % savefile

        fov = n.degrees(1./d['res'])*60.
        ypix = len(candsprop[i][1])
        xpix = len(candsprop[i][1][0])
        if len(candsprop)> 0:
            fig = plt.Figure()
            ax = fig.add_subplot(111)
            im = ax.imshow(imsum, aspect='auto', origin='lower', interpolation='nearest', extent=[fov/2, -fov/2, -fov/2, fov/2])    
            ax.scatter(((ypix/2-impeak[:,1]))*fov/ypix, (impeak[:,0]-xpix/2)*fov/xpix, s=80, marker='o', facecolor='none')
            ax.set_xlabel('Dec Offset (amin)')
            ax.set_ylabel('RA Offset (amin)')
            ax.set_xlim(fov/2, -fov/2)
            ax.set_ylim(-fov/2, fov/2)
            canvas = FigureCanvasAgg(fig)
            canvas.print_figure(savefile)
        else:
            print 'No candidates to plot.'

    return impeak

def dmtplot(d, candsloc, snrs, saveroot=''):
    """ Take loc and snr lists to plot dm vs time with symbol size related to snr.
    """

    if saveroot:
        savefile = saveroot + '_dmt.png'
        print 'Plotting to %s.' % savefile

        fig = plt.Figure()
        ax = fig.add_subplot(111)
        if len(candsloc)> 0:
            ax.scatter(d['inttime']*candsloc[:,2], n.array(d['dmarr'])[candsloc[:,3]], s=2*snrs**2, facecolor='none', alpha=0.5, clip_on=False)
            ax.set_xlabel('Time (s)')
            ax.set_ylabel('DM (pc/cm3)')
            canvas = FigureCanvasAgg(fig)
            canvas.print_figure(savefile)
        else:
            print 'No candidates to plot.'

def locpluspeak(pklroot):
    """ Get all candidates from pklroot selection for compiling location with peak pixel (beam, dt, dm, time, x, y)
    """

    # need to get dictionary for imaging detail
    pkllist = glob.glob(pklroot)
    pkl = open(pkllist[0], 'rb')
    d = pickle.load(pkl)

    # get candidates and image peaks
    cands = filter_candidates(pklroot)
    if len(cands) == 2:
        candsloc, candsprop = cands
    elif len(cands) == 3:
        dataloc, candsloc, candsprop = cands

    impeaks = impeakplot(d, candsprop, saveroot='')
    
    candslocplus = n.concatenate((candsloc,impeaks), axis=1)

    return candslocplus

def snrhist(snr, saveroot=''):
    """ Plot SNR histogram of candidates. Expects snr list.
    """

    if saveroot:
        savefile = saveroot + '_snrhist.png'
        print 'Plotting to %s.' % savefile

        fig = plt.Figure()
        ax = fig.add_subplot(111)
        ax.hist(snr, bins=15, normed=True, histtype='step', alpha=0.5)
        ax.set_xlabel('SNR')
        ax.set_ylabel('Number')
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(savefile)
    return snr

def dmtpixplot(pklroot, mind, maxd):
    """ Get all candidates from pklroot (expected to include multiple scans and days) and look for dm-pix neighbors.
    mind and maxd define range of distances in dm,pix space to return.
    """

    locpl = locpluspeak(pklroot)

    # fill difference array
    diff = n.zeros( (len(locpl)**2,4), dtype='int')
    for i in range(len(locpl)):
        print 'Filling diff array for index ', i
        diff[len(locpl)*i:len(locpl)*(i+1)] = locpl[i,2:] - locpl[:,2:]    # include int,dm,x,y

    dist = n.sqrt(diff[:,1]**2 + diff[:,2]**2 + diff[:,3]**2)   # not really cartesian, but whatever
    ij = n.where((dist > mind) & (dist < maxd))[0]
    j = n.mod(ij, len(locpl))   # index of second
    i = (ij - j)/len(locpl)     # index of first

    print 'Candidates:'
    for ind in range(len(j)):
        print '*-%d-%d-*png and *-%d-%d-*png' % (locpl[i[ind],2], locpl[i[ind],3], locpl[j[ind],2], locpl[j[ind],3])

def singleplot(pklfile, candsloc, candsprop, i, save=1):
    """ Plot a single candidate out of a list.
    """

    candsloc = n.array(candsloc).astype(int)
    print 'Plotting for candidate at: ', candsloc[i]

    pkl = open(pklfile, 'rb')
    d = pickle.load(pkl)
    pkl.close()

    snr = n.array([cand[0] for cand in candsprop])

    fig = plt.Figure(figsize=(8.5,8))

    # text description of candidate
    ax = fig.add_subplot(221, axisbg='white')
#    ax.set_title('Candidate @ Tier 1')
    # first plot dm-t distribution beneath
    ax.scatter(d['inttime']*candsloc[:,2], n.array(d['dmarr'])[candsloc[:,3]], s=(snr-6)**2, facecolor='none', linewidth=0.05, clip_on=False)
    ax.scatter(d['inttime']*candsloc[i,2], n.array(d['dmarr'])[candsloc[i,3]], s=(snr[i]-6)**2, marker='x', facecolor='none', linewidth=2, clip_on=False)
    ax.set_ylim(0, 4000)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('DM (pc/cm3)')
    # then add annotating info
    beamra = n.round(d['delaycenters'][candsloc[i,0]][0], 1)
    beamdec = n.round(d['delaycenters'][candsloc[i,0]][1], 1)
    fov = n.degrees(1./d['res'])*60.
    xpix,ypix = candsprop[i][1].shape
    srcra,srcdec = n.where(candsprop[i][1] == candsprop[i][1].max())
    ax.text(0.1, 0.9, d['filename'], fontname='sans-serif', transform = ax.transAxes)
    ax.text(0.1, 0.8, 'dt %d, Int %d, DM %.1f' % (d['dtarr'][candsloc[i,1]], candsloc[i,2], d['dmarr'][candsloc[i,3]]), fontname='sans-serif', transform = ax.transAxes)
    ax.text(0.1, 0.7, 'Peak: (' + str(n.round((xpix/2-srcra[0])*fov/xpix, 1)) + '\' ,' + str(n.round((ypix/2-srcdec[0])*fov/ypix, 1)) + '\'), SNR: ' + str(n.round(snr[i], 1)), fontname='sans-serif', transform = ax.transAxes)

    # plot dynamic spectra
    left, width = 0.6, 0.2
    bottom, height = 0.2, 0.7
    rect_dynsp = [left, bottom, width, height]
    rect_lc = [left, bottom-0.1, width, 0.1]
    rect_sp = [left+width, bottom, 0.1, height]
    ax_dynsp = fig.add_axes(rect_dynsp)
    ax_lc = fig.add_axes(rect_lc)
    ax_sp = fig.add_axes(rect_sp)
    spectra = n.swapaxes(candsprop[i][2].data.real,0,1)      # seems that latest pickle actually contains complex values in spectra...
    dd = n.concatenate( (spectra[...,0], n.zeros_like(spectra[...,0]), spectra[...,1]), axis=1)    # make array for display with white space between two pols
    im = ax_dynsp.imshow(dd, origin='lower', interpolation='nearest', aspect='auto', cmap=plt.get_cmap('Greys'))
    ax_dynsp.text(0.5, 0.95, 'RR LL', horizontalalignment='center', verticalalignment='center', fontsize=16, color='w', transform = ax_dynsp.transAxes)
    ax_dynsp.set_yticks(range(0,len(d['freq']),30))
    ax_dynsp.set_yticklabels(d['freq'][::30])
    ax_dynsp.set_ylabel('Freq (GHz)')
    spectrum = spectra[:,len(spectra[0])/2].mean(axis=1)      # assume pulse in middle bin. get stokes I spectrum. **this is wrong in a minority of cases.**
    ax_sp.plot(spectrum, range(len(spectrum)), 'k.')
    ax_sp.plot(n.zeros(len(spectrum)), range(len(spectrum)), 'k:')
    ax_sp.set_ylim(0, len(spectrum))
    ax_sp.set_yticklabels([])
    xmin,xmax = ax_sp.get_xlim()
    ax_sp.set_xticks(n.linspace(xmin,xmax,3).round(2))
#    ax_sp.set_xticklabels()    # set labels to be rounded?
    ax_sp.set_xlabel('Flux (Jy)')
#    lc = spectra.mean(axis=2).mean(axis=0)    # to show stokes I
    lc = dd.mean(axis=0)
    lenlc = n.where(lc == 0)[0][0]
    ax_lc.plot(range(0,lenlc)+range(2*lenlc,3*lenlc), list(lc)[:lenlc] + list(lc)[-lenlc:], 'k.')
    ax_lc.plot(range(0,lenlc)+range(2*lenlc,3*lenlc), list(n.zeros(lenlc)) + list(n.zeros(lenlc)), 'k:')
    ax_lc.set_xlabel('Integration')
    ax_lc.set_ylabel('Flux (Jy)')
    ax_lc.set_xticks([0,0.5*lenlc,lenlc,1.5*lenlc,2*lenlc,2.5*lenlc,3*lenlc])
    ax_lc.set_xticklabels(['0',str(lenlc/2),str(lenlc),'','0',str(lenlc/2),str(lenlc)])
    ymin,ymax = ax_lc.get_ylim()
    ax_lc.set_yticks(n.linspace(ymin,ymax,3).round(2))

#    ax.plot(d['freq'], spectrum, '.')
#    ax.set_xlabel('Frequency (GHz)')
#    ax.set_ylabel('Flux density (Jy)')
#    sm = n.sqrt( ((spectrum**2).mean() - spectrum.mean()**2) / spectrum.mean()**2 )
#    ax.text(0.05, 0.05, 'specmod =' + str(n.round(sm, 1)), transform = ax.transAxes)
    
    # dm-time of candidates in provided list
    ax = fig.add_subplot(223)
    ax.scatter(((xpix/2-srcra[0])-0.05*xpix)*fov/xpix, (ypix/2-srcdec[0])*fov/ypix, s=80, marker='<', facecolor='none')
    ax.scatter(((xpix/2-srcra[0])+0.05*xpix)*fov/xpix, (ypix/2-srcdec[0])*fov/ypix, s=80, marker='>', facecolor='none')
    im = ax.imshow(candsprop[i][1].transpose(), aspect='equal', origin='upper', interpolation='nearest', extent=[fov/2, -fov/2, -fov/2, fov/2], cmap=plt.get_cmap('Greys'), vmin=0, vmax=0.5*candsprop[i][1].max())
#    ax.scatter(srcra[0]+20, srcdec[0], s=70, marker='<', facecolor='none')
#    ax.imshow(candsprop[i][1], aspect='auto', origin='lower', interpolation='nearest')
#    im.colorbar()
    ax.set_xlabel('RA Offset (arcmin)')
    ax.set_ylabel('Dec Offset (arcmin)')

    if save:
        candslocstr = '_tier1-%d-%d-%d-%d' % (candsloc[i,0],candsloc[i,1],candsloc[i,2],candsloc[i,3])
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(string.join(pklfile.split('.')[:-1], '.') + str(candslocstr) + '.png')

def singleplot2(pklfile, threshold, candnum, im, data, savefile=''):
    """ Plot a single candidate with image and data from reproducecanddata
    """

    pkl = open(pklfile, 'rb')
    d = pickle.load(pkl)
    dataloc0, candsloc0, candsprop0 = pickle.load(pkl)
    pkl.close()
    dataloc0 = n.array(dataloc0); candsloc0 = n.array(candsloc0); candsprop0 = n.array(candsprop0)

    dataloc, candsloc, candsprop = filter_candidates(pklfile, threshold)
    dataloc = n.array(dataloc); candsloc = n.array(candsloc); candsprop = n.array(candsprop)
    scani = [int(name.rstrip('.pkl').split('_')[3][1:]) for name in d['pkllist']]
    filesplit = d['pkllist'][scani.index(dataloc[candnum])].split('_')
    filename = '_'.join(filesplit[1:4]) + '.ms'
    snrcol = parsepropinfo(candsprop)
    snr0 = n.array([prop[snrcol] for prop in candsprop0])
    snr = n.array([prop[snrcol] for prop in candsprop])

    changepoints = []
    for i in range(len(dataloc0)-1):
        if (dataloc0[i+1] - dataloc0[i]):
            changepoints.append(dataloc0[i+1])
    changepoints = n.array(changepoints)

    fig = plt.Figure(figsize=(8.5,8))
    # text description of candidate
    ax = fig.add_subplot(221, axisbg='white')
#    ax.set_title('Candidate @ Tier 1')
    # first plot dm-t distribution beneath
    times0 = d['inttime']*(candsloc0[:,2] + d['nints']*dataloc0)
    times = d['inttime']*(candsloc[candnum,2] + d['nints']*dataloc[candnum])
    ax.scatter(times0, n.array(d['dmarr'])[candsloc0[:,3]], s=(snr0-0.9*snr0.min())**5, facecolor='none', linewidth=0.05, clip_on=False)
    ax.scatter(times, n.array(d['dmarr'])[candsloc[candnum,3]], s=(snr-0.9*snr0.min())**5, marker='x', facecolor='none', linewidth=2, clip_on=False)
    ax.set_ylim(0, 4500)
    ax.set_xlim(times0.min(), times0.max())
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('DM (pc/cm3)')
    ax.set_xticks(n.linspace(times0.min().astype('int'), times0.max().astype('int'), 2))
#    ax.set_xticklabels()

    # then add annotating info
    beamra = n.round(d['delaycenters'][candsloc[candnum,0]][0], 1)
    beamdec = n.round(d['delaycenters'][candsloc[candnum,0]][1], 1)
    fov = n.degrees(1./d['res'])*60.
    xpix,ypix = im.shape
    srcra,srcdec = n.where(im == im.max())
    ax.text(0.1, 0.9, filename, fontname='sans-serif', transform = ax.transAxes)
    ax.text(0.1, 0.8, 'dt %d, Int %d, DM %.1f' % (d['dtarr'][candsloc[candnum,1]], candsloc[candnum,2], d['dmarr'][candsloc[candnum,3]]), fontname='sans-serif', transform = ax.transAxes)
    ax.text(0.1, 0.7, 'Peak: (' + str(n.round((xpix/2-srcra[0])*fov/xpix, 1)) + '\' ,' + str(n.round((ypix/2-srcdec[0])*fov/ypix, 1)) + '\'), SNR: ' + str(n.round(snr[candnum], 1)), fontname='sans-serif', transform = ax.transAxes)

    # plot dynamic spectra
    left, width = 0.6, 0.2
    bottom, height = 0.2, 0.7
    rect_dynsp = [left, bottom, width, height]
    rect_lc = [left, bottom-0.1, width, 0.1]
    rect_sp = [left+width, bottom, 0.1, height]
    ax_dynsp = fig.add_axes(rect_dynsp)
    ax_lc = fig.add_axes(rect_lc)
    ax_sp = fig.add_axes(rect_sp)
    spectra = n.swapaxes(data.real,0,1)      # seems that latest pickle actually contains complex values in spectra...
    dd = n.concatenate( (spectra[...,0], n.zeros_like(spectra[...,0]), spectra[...,1]), axis=1)    # make array for display with white space between two pols
    impl = ax_dynsp.imshow(dd, origin='lower', interpolation='nearest', aspect='auto', cmap=plt.get_cmap('Greys'))
    ax_dynsp.text(0.5, 0.95, 'RR LL', horizontalalignment='center', verticalalignment='center', fontsize=16, color='w', transform = ax_dynsp.transAxes)
    ax_dynsp.set_yticks(range(0,len(d['freq']),30))
    ax_dynsp.set_yticklabels(d['freq'][::30])
    ax_dynsp.set_ylabel('Freq (GHz)')
    ax_dynsp.set_xlabel('Integration (rel)')
    spectrum = spectra[:,len(spectra[0])/2].mean(axis=1)      # assume pulse in middle bin. get stokes I spectrum. **this is wrong in a minority of cases.**
    ax_sp.plot(spectrum, range(len(spectrum)), 'k.')
    ax_sp.plot(n.zeros(len(spectrum)), range(len(spectrum)), 'k:')
    ax_sp.set_ylim(0, len(spectrum))
    ax_sp.set_yticklabels([])
    xmin,xmax = ax_sp.get_xlim()
    ax_sp.set_xticks(n.linspace(xmin,xmax,3).round(2))
#    ax_sp.set_xticklabels()    # set labels to be rounded?
    ax_sp.set_xlabel('Flux (Jy)')
#    lc = spectra.mean(axis=2).mean(axis=0)    # to show stokes I
    lc = dd.mean(axis=0)
    lenlc = n.where(lc == 0)[0][0]
    ax_lc.plot(range(0,lenlc)+range(2*lenlc,3*lenlc), list(lc)[:lenlc] + list(lc)[-lenlc:], 'k.')
    ax_lc.plot(range(0,lenlc)+range(2*lenlc,3*lenlc), list(n.zeros(lenlc)) + list(n.zeros(lenlc)), 'k:')
    ax_lc.set_xlabel('Integration')
    ax_lc.set_ylabel('Flux (Jy)')
    ax_lc.set_xticks([0,0.5*lenlc,lenlc,1.5*lenlc,2*lenlc,2.5*lenlc,3*lenlc])
    ax_lc.set_xticklabels(['0',str(lenlc/2),str(lenlc),'','0',str(lenlc/2),str(lenlc)])
    ymin,ymax = ax_lc.get_ylim()
    ax_lc.set_yticks(n.linspace(ymin,ymax,3).round(2))

#    ax.plot(d['freq'], spectrum, '.')
#    ax.set_xlabel('Frequency (GHz)')
#    ax.set_ylabel('Flux density (Jy)')
#    sm = n.sqrt( ((spectrum**2).mean() - spectrum.mean()**2) / spectrum.mean()**2 )
#    ax.text(0.05, 0.05, 'specmod =' + str(n.round(sm, 1)), transform = ax.transAxes)
    
    # dm-time of candidates in provided list
    ax = fig.add_subplot(223)
    ax.scatter(((xpix/2-srcra[0])-0.05*xpix)*fov/xpix, (ypix/2-srcdec[0])*fov/ypix, s=80, marker='<', facecolor='none')
    ax.scatter(((xpix/2-srcra[0])+0.05*xpix)*fov/xpix, (ypix/2-srcdec[0])*fov/ypix, s=80, marker='>', facecolor='none')
    impl = ax.imshow(im.transpose(), aspect='equal', origin='upper', interpolation='nearest', extent=[fov/2, -fov/2, -fov/2, fov/2], cmap=plt.get_cmap('Greys'), vmin=0, vmax=0.5*im.max())
    ax.set_xlabel('RA Offset (arcmin)')
    ax.set_ylabel('Dec Offset (arcmin)')

    if savefile:
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(savefile)

def readnoise(noisefile):
    f = open(noisefile,'r')
    itercount = []; noiseperbl = []; flagfrac = []; imnoise = []
    while True:         #for each file it loops through to get 'y' column data
        try:
            value=pickle.load(f)
            itercount.append(value[0])
            noiseperbl.append(value[1])
            flagfrac.append(value[2])
            imnoise.append(value[3]) # fixed to image pix
        except EOFError:           #when hits EOFError it breaks
            f.close()
            break
    return (n.array(itercount), n.array(noiseperbl), n.array(flagfrac), n.array(imnoise))

def noisehists(pklroot, savefile='', remove={}):
    """ Cumulative hist of image noise levels.
    """

    pkllist = glob.glob(pklroot)
    noises = []; noisec = []
    print 'Reading %d noise files' % len(pkllist)
    for pkl in pkllist:
        nn = readnoise(pkl)
        ii = nn[0]
        nn = list(nn[3])

        scani = int(pkl.split('_s')[1].split('.')[0])   # assumes scan name structure
        if scani in remove.keys():
            print 'Removing some noise measurements from ', pkl
            nranges = len(remove[scani])
            wwa = []
            for first in range(0,nranges,2):
                badrange0 = remove[scani][first]
                badrange1 = remove[scani][first+1]
                ww = list(n.where( (ii > badrange0) & (ii < badrange1) )[0])
                if len(ww):
                    wwa += ww
            for i in wwa[::-1]:
                junk = nn.pop(i)
        noises.append(nn)
        noisec = noisec + nn

    bins = n.linspace(min(noisec), max(noisec), 50)
    fig = plt.Figure(figsize=(10,10))
    ax = fig.add_subplot(211, axisbg='white')
    stuff = ax.hist(noises, bins=bins, histtype='bar', lw='none', ec='none')
    ax.set_title('Histograms of noise samples')
    ax.set_xlabel('Image RMS (Jy)')
    ax.set_ylabel('Number of noise measurements')
    ax2 = fig.add_subplot(212, axisbg='white')
    stuff = ax2.hist(noisec, bins=bins, cumulative=-1, normed=True, log=False, histtype='bar', lw='none', ec='none')
    ax2.set_xlabel('Image RMS (Jy)')
    ax2.set_ylabel('Number with noise > image RMS')

    if savefile:
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(savefile)

def candplot0(pklfile, threshold=0, plotnum=-1, save=1):
    """ Summary plot of candidate.
    This is a "tier 1" plot, since it only uses info saved in real time.
    Does not support island detection yet.
    save defines whether plots are also saved.
    """

    print 'Building Tier 1 candidate plot.'

    pkl = open(pklfile, 'rb')
    d = pickle.load(pkl)
    pkl.close()

    cands = filter_candidates(pklfile, threshold, saveroot='')
    if len(cands) == 2:
        candsloc, candsprop = cands
    elif len(cands) == 3:
        dataloc, candsloc, candsprop = cands
    if len(candsloc) == 0:
        print 'No candidates available...'
        return

    if plotnum == -1:
        canditer = range(len(candsloc))
        print 'Plotting all candidates.'
    else:
        canditer = [plotnum]
        print 'Plotting candidate number %d.' % plotnum

    # plot candidate info
    for i in canditer:
        singleplot(pklfile, candsloc, candsprop, i, save=1)

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
    (candsloc, candsprop) = pickle.load(pkl)
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
                ww = n.where( (n.array(candsloc[beam])[:,0]==island[i,0]) & (n.array(candsloc[beam])[:,1]==island[i,1]) & (n.array(candsloc[beam])[:,2]==island[i,2]) )
                islandind.append(ww[0][0])
            loc = n.array(candsloc[beam])[islandind]
            snr = n.array(candsprop[beam])[islandind, 0]

            # then plot snr as function of dmind and dtind for island
            fixed_dt = n.where(loc[:,0] == cand[0])
            fixed_dm = n.where(loc[:,2] == cand[2])
            dmdist = n.squeeze(n.array(d['dmarr'])[loc[fixed_dt, 2]])   # seems to have superfluous axis...?
            dmsnr = snr[fixed_dt]
            dtdist = n.squeeze(n.array(d['dtarr'])[loc[fixed_dm][:, 0]])   # seems to have superfluous axis...?
            dtsnr = snr[fixed_dm]

            # find cands index of candidate
            ww = n.where([n.all(candscand == cand) for candscand in candsloc[beam]])[0][0]

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
            p.text(0.1, 0.6, 'Beam: ' + str(beamra) + ', ' + str(beamdec) + ', dt: ' + str(d['dtarr'][cand[0]]), fontname='sans-serif')
            p.text(0.1, 0.4, 'Integration: ' + str(cand[1]) + ', DM: ' + str(d['dmarr'][cand[2]]), fontname='sans-serif')
            p.text(0.1, 0.2, 'SNR: ' + str(n.round(candsprop[beam][ww][0], 1)), fontname='sans-serif')

            p.subplot(222)
            p.plot(dmdist, dmsnr, 'b.', clip_on=False, label='DM at peak dt')
            p.xlabel('DM (pc/cm3)')
            p.ylabel('SNR')
            p.twiny()
            p.plot(dtdist, dtsnr, 'r+', clip_on=False, label='dt at peak DM')
            p.xlabel('dt (ints)')
            p.subplot(223)
            p.plot(candsdata[beam][ww][0], 'b.', label='Mean B')
            p.xlabel('Integration')
            p.ylabel('Mean, Std of Bispectra')
            p.twinx()
            p.plot(candsdata[beam][ww][1], 'r.', label='Std B')
            p.subplot(224)
            dataph = n.rot90(candsdata[beam][ww][2])
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

            # process!
            noiseperbl = estimate_noiseperbl(data)
            lib.dedisperse(data, d, d['dmarr'][dmind], verbose=1)        # dedisperses 'data' in place
            lib.phaseshift(data, d, n.radians(ra), n.radians(dec), u, v)            # shift phase center to candidate beam
#            time_filter(self.timescales[dtind], self.filtershape, bgwindow=self.bgwindow)    # tophat filter (bgwindow not used)
#            noiseperbl = self.obs.data.mean(axis=3).mean(axis=2).real.std()   # measure single noise for input to detect_bispectra

            # reproduce bispectrum results
            bispectra = lib.make_bispectra(data, triples)
            bispectra = n.ma.masked_array(bispectra, bispectra == 0j)   # set flags for zero data
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
            p.figure(figsize=(12,9))
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
            p.text(0.1, 0.6, 'Beam: ' + str(beamra) + ', ' + str(beamdec) + ', dt: ' + str(d['dtarr'][cand[0]]), fontname='sans-serif')
            p.text(0.1, 0.4, 'Integration: ' + str(cand[1]) + ', DM: ' + str(d['dmarr'][cand[2]]), fontname='sans-serif')
            p.text(0.1, 0.2, 'SNR: ' + str(n.round(candsnr, 1)), fontname='sans-serif')

            # image of dispersed visibilities for candidate
            p.subplot(222)
            fov = n.degrees(1./d['res'])*60.       # field of view in arcseconds
            p.imshow(im, aspect='auto', origin='lower', interpolation='nearest', extent=[fov/2, -fov/2, -fov/2, fov/2])
            p.colorbar()
            p.xlabel('RA/l Offset (arcmin)')
            p.ylabel('Dec/m Offset (arcmin)')

            p.subplot(223)
            dataph = spec(data, 0, len(data))
            p.plot(d['freq'], dataph[twindow/2], '.')
            p.text(0.05, 0.05, 'specmod =' + str(sm), transform = ax.transAxes)
            p.xlabel('Flux density (Jy)')
            p.ylabel('Frequency (GHz)')

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

def find_candidates(pklfile, d_neighbor=2, island_size_min=1, island_snr_min=5., show='', save=0):
    """ Visualize and define candidates to follow up from pickle file with search state and candidate dictionary.
    d_neighbor is the dm-time distance over which a neighbor is defined.
    island_size_min is the minimum number of detections requried to call an island a candidate
    island_snr_ min is the minimum snr requried to call an island a candidate
    returns the SNR peak of each island identified with d_neighbor.
    """

    # read in pickle file of candidates
    pkl = open(pklfile, 'rb')
    d = pickle.load(pkl)
    (candsloc, candsprop) = pickle.load(pkl)
    pkl.close()

    print d['filename'], d['spw'], d['chans'], d['dmarr'], d['sigma_image`']

    neighbors = {}
    for (ra,dec) in d['delaycenters']:
        neighbors[(ra,dec)] = []

    dmarr = n.array(d['dmarr'])

    # measure number of neighbors
    for beam in candsloc.keys():
        neighbors_lm = []
        cands_lm = candsloc[beam]
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
    for (ra,dec) in candsloc.keys():
        islands[(ra,dec)] = []
    for beam in candsloc.keys():
        cands_lm = candsloc[beam]
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

    print islands
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
                    islandind.append(n.where( (n.array(candsloc[beam])[:,0]==island[i,0]) & (n.array(candsloc[beam])[:,1]==island[i,1]) & (n.array(candsloc[beam])[:,2]==island[i,2]) )[0][0])
                maxsnrind = n.where(n.array(candsprop[beam])[islandind] == n.max(n.array(candsprop[beam])[islandind]))  # if snr is only thing in cands prop
#                maxsnrind = n.where(n.array(candsprop[beam])[islandind,0] == n.max(n.array(candsprop[beam])[islandind,0])) # if candsprop has more than one value
                if n.array(candsprop[beam])[islandind][maxsnrind][0] > island_snr_min:       # if snr is only prop
#                if n.array(candsprop[beam])[islandind][maxsnrind][0][0] > island_snr_min:   # more than just snr
                    islandmaxlocl.append(n.array(candsloc[beam])[islandind][maxsnrind][0])
                    islandmaxsnrl.append(n.array(candsprop[beam])[islandind][maxsnrind][0])    # if snr is only prop
#                    islandmaxsnrl.append(n.array(candsprop[beam])[islandind][maxsnrind][0][0])   # more than just snr
        islandmaxloc[beam] = n.array(islandmaxlocl).astype(int)
        islandmaxsnr[beam] = n.array(islandmaxsnrl)

    if show or save:
        cm = p.get_cmap('gist_rainbow')
        p.figure(1, figsize=(12,9))
        tl1 = p.subplot(221)
        tr1 = p.subplot(222)
        bl1 = p.subplot(223)
        br1 = p.subplot(224)

        beamind = 0
        for beam in candsloc.keys():  # iterate over beam candidate groups
            loc = n.array(candsloc[beam]).astype(int)
            prop = n.array(candsprop[beam])
            shiftind = float(beamind)/len(candsloc.keys())   # shift point location (and color) per beam for clarity
            if len(loc) == 0: break
            p.figure(1)
            p.subplot(tl1)
            p.scatter(d['inttime']*loc[:,1], dmarr[loc[:,2]] + 0.5*(dmarr[1]-dmarr[0])*shiftind, s=8*prop, facecolor='none', color=cm(shiftind), alpha=0.5, clip_on=False)  # if prop only has snr
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
            p.scatter(d['inttime']*loc[:,1], prop, color=cm(shiftind), alpha=0.3, clip_on=False)   # prop only has snr
#            p.scatter(d['inttime']*loc[:,1], prop[:,0], color=cm(shiftind), alpha=0.3, clip_on=False)
            p.scatter(d['inttime']*islandmaxloc[beam][:,1], islandmaxsnr[beam], color=cm(shiftind), marker='+', s=100, alpha=0.8, clip_on=False)
            axis_bl1 = p.axis()
            axis_bl1 = p.axis([axis_tl1[0], axis_tl1[1], axis_bl1[2], axis_bl1[3]])

            p.subplot(br1)
            p.hist(prop, orientation='horizontal', color=cm(shiftind), histtype='step', alpha=0.5, bins=len(prop)/2)  # prop only has snr
#            p.hist(prop[:,0], orientation='horizontal', color=cm(shiftind), histtype='step', alpha=0.5, bins=len(prop)/2)
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

def summarize_beams(pklname):
    """ Reads pickle file and summarizes where candidates are in beams.
    """

    with open(pklname, 'r') as pkl:
        shift = 0
        beams = []
        d = pickle.load(pkl)
        (candsloc, candsprop) = pickle.load(pkl)
        candsloc = n.array(candsloc)
        print 'Candidate beam, location, properties:'
        for beam in candsloc[:,0].unique():
            ww = n.where(candsloc[:,0] == beam)
            print beam, candsloc[ww]
            beams.append(beam)
            for cand in candsloc[ww]:
                p.text(beam[0]+shift, beam[1]+shift, str(cand[0]) + '-' + str(cand[1]) + '-' + str(cand[2]), horizontalalignment='center', verticalalignment='center', fontsize=9)
                shift += 0.00003

        for beam in beams:
            p.plot([beam[0], beam[0]+shift], [beam[1], beam[1]+shift], 'k-')

    p.show()

