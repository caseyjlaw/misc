import sdmreader, ephem
import numpy as n
import os, math, pickle, logging
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib import pyplot as plt
from functools import partial
from contextlib import closing
import multiprocessing as mp

# parameters
logging.basicConfig(level=logging.INFO)
nthread = 8
sdmname = 'TSKY0001.sb30889013.eb30992856.57214.70375685186'
roll = 512
outname = 'RFImap_57214.png'
intchunk = 8   # time size (in ints) to read. 8 is typical size 
spwfreq = ['2.00 GHz', '2.12 GHz', '2.24 GHz', '2.37 GHz', '2.50 GHz', '2.63 GHz', '2.76 GHz', '2.88 GHz', '3.00 GHz', '3.12 GHz', '3.24 GHz', '3.37 GHz', '3.50 GHz', '3.63 GHz', '3.76 GHz', '3.88 GHz']

def read_metadata(sdmname):
    """ Read scan and source metadata
    """

    sc,sr = sdmreader.read_metadata(sdmname)
    return sc, sr

def all(scanlist=[], pklname='rfiarray_57214.pkl'):
    """ Wrapper for functions to create or read data to visualize RFI
    pklname, if exists, will be read. Otherwise, values read/calculated and written to pklname.
    """

    if os.path.exists(pklname):
        logging.info('Found %s. Reading...' % pklname)
        with open(pklname, 'r') as pkl:
            rfiarray, altaz = pickle.load(pkl)
    else:
        logging.info('Reading/calculating RFI data...')
        rfiarray = read(scanlist)
        altaz = calcaltaz(scanlist)
        if pklname:   # if actually defined, then write there
            with open(pklname, 'w') as pkl:
                pickle.dump((rfiarray,altaz), pkl)

    plotrfi(rfiarray, altaz)

def read(scanlist):
    # iterate through scans (i.e., time) and build RFI summaries
    rfilist = []
    results = []
    dim0 = 0
    with closing(mp.Pool(nthread)) as readpool:
        for scan in scanlist:
            logging.info('Reading scan %d, %s' % (scan, sc[scan]))
            for nskip in range(0, sc[scan]['nints']-intchunk+1, intchunk):                   # iterate over intchunks (whole segments only, dropping end)
                results.append(readpool.apply_async(readreduce, [sdmname, scan, nskip]))

        logging.info('Enqueued %d jobs...' % len(results))
        for result in results:
            rfi = result.get()
            dim0 += len(rfi)
            rfilist.append(rfi)

    # aggregate summaries into array
    dim1 = rfi.shape[1]
    dtype = rfi.dtype
    rfiarray = n.zeros((dim0,dim1), dtype=dtype)
    loc0 = 0
    for i in range(len(rfilist)):
        rfi = rfilist[i]
        rfiarray[loc0:loc0+len(rfi)] = rfi
        loc0 += len(rfi)

    return rfiarray

def readreduce(sdmname, scan, nskip):
    def reducedata(data):
        return n.abs(data).max(axis=3).max(axis=1).max(axis=0)[None,:]    # 2d array, but collapsed by taking max over time
#        return n.abs(data).max(axis=3).max(axis=1)    # returns 2d array with all ints

    data = n.roll(sdmreader.read_bdf(sdmname, scan, nskip=nskip, readints=intchunk, writebdfpkl=True), roll, axis=2)   # roll needed to get increasing freq spw order
    return reducedata(data)

def calcaltaz(scanlist, sc, sr, format='str'):
    """ Calculates a single (alt,az) per scan in scanlist.
    """

    inttime = (sc[1]['endmjd'] - sc[1]['startmjd'])*24*3600/sc[1]['nints']

    vla = ephem.Observer()
    vla.lat = '34:04:43.497'
    vla.long = '-107:37:03.819'
    vla.elevation = 2124
    src = ephem.FixedBody()

    altaz = []
    for scan in scanlist:
        src._ra, src._dec = [(sr[srn]['ra'], sr[srn]['dec']) for srn in sr.keys() if sc[scan]['source'] == sr[srn]['source']][0]
        for nskip in range(0, sc[scan]['nints']-intchunk+1, intchunk):
            vla.date = ephem.date(jd_to_date(sc[scan]['startmjd'] + nskip*inttime/(24*3600) + 2400000.5))
            src.compute(vla)
            if format == 'str':
                altaz.append( '(%.1f, %.1f)' % (n.degrees(src.alt), n.degrees(src.az)) )
            elif format == 'float':
                altaz.append( (n.degrees(src.alt), n.degrees(src.az)) )

    return n.array(altaz)

def calcradec(scanlist, sc, sr):
    radec = []
    for scan in scanlist:
        for nskip in range(0, sc[scan]['nints']-intchunk+1, intchunk):
            radec.append([(n.degrees(sr[srn]['ra']), n.degrees(sr[srn]['dec'])) for srn in sr.keys() if sc[scan]['source'] == sr[srn]['source']][0])

    return n.array(radec)

def plotrfi(rfiarray, altaz):
    # plot
    fig,ax = plt.subplots(figsize=(12,10))
    im = ax.imshow(n.log(rfiarray), interpolation='nearest', origin='lower', cmap=plt.get_cmap('cubehelix'), aspect='auto')
    fig.colorbar(im)

    # nice ticks and labels
    ax.set_xticks(range(0,1024,64))
    ax.set_xticklabels(spwfreq)
    plt.setp( ax.xaxis.get_majorticklabels(), rotation=80)
    ax.set_xlabel('Frequency')
    ax.set_yticks(range(0,len(rfiarray),min(len(rfiarray), 150)))
    ax.set_yticklabels(altaz[::min(len(rfiarray), 150)])
    ax.set_ylabel('(Alt, Az)')

    # save to fig
    canvas2 = FigureCanvasAgg(fig)
    canvas2.print_figure(outname)


def jd_to_date(jd):
    """
    Convert Julian Day to date.
    
    Algorithm from 'Practical Astronomy with your Calculator or Spreadsheet', 
        4th ed., Duffet-Smith and Zwart, 2011.
    
    Parameters
    ----------
    jd : float
        Julian Day
        
    Returns
    -------
    year : int
        Year as integer. Years preceding 1 A.D. should be 0 or negative.
        The year before 1 A.D. is 0, 10 B.C. is year -9.
        
    month : int
        Month as integer, Jan = 1, Feb. = 2, etc.
    
    day : float
        Day, may contain fractional part.
        
    Examples
    --------
    Convert Julian Day 2446113.75 to year, month, and day.
    
    >>> jd_to_date(2446113.75)
    (1985, 2, 17.25)
    
    """
    jd = jd + 0.5
    
    F, I = math.modf(jd)
    I = int(I)
    
    A = math.trunc((I - 1867216.25)/36524.25)
    
    if I > 2299160:
        B = I + 1 + A - math.trunc(A / 4.)
    else:
        B = I
        
    C = B + 1524
    
    D = math.trunc((C - 122.1) / 365.25)
    
    E = math.trunc(365.25 * D)
    
    G = math.trunc((C - E) / 30.6001)
    
    day = C - E + F - math.trunc(30.6001 * G)
    
    if G < 13.5:
        month = G - 1
    else:
        month = G - 13
        
    if month > 2.5:
        year = D - 4716
    else:
        year = D - 4715
        
    return year, month, day
    
