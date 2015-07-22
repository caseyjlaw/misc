#
# script to do real-time transient detection on sdms found on cbe or aoc cluster
# claw, 13aug13

import numpy as n
from casa import importasdm, ms, split
import makeasdm, leanpipe
import datetime, os, sys, string, shutil
import xml.etree.ElementTree as et

# important paths on cbe
#workdir = '/lustre/evla/test/claw/'   # on cbe node
workdir = './'   # on postprocessing cluster
telcaldir = ''
#telcaldir = '/home/mchammer/evladata/telcal/'  # followed by yyyy/mm/sdmname.GN
datadir = '/home/mchammer/evla/mcaf/workspace/'  # where to find new sdms
bdfdir = '/lustre/evla/wcbe/data/bunker/'        # not present on aoc cluster, but would be ignored if binaries present

def prep(asdmfile, datadir=datadir, telcaldir=telcaldir, bdfdir=bdfdir, workdir=workdir, intentfilter="OBSERVE_TARGET"):
    """ Fills MS file and gets .GN telcal file.
    Takes asdm file name and derives all from that. 
    Copies to working area and makes it to proper asdm file, if needed.
    """

    telcalfile = asdmfile + '.GN'
    msfile = asdmfile + '.ms'

    os.chdir(workdir)
    if os.path.exists(msfile):
        print 'Found MS file %s' % (workdir+msfile)
    else:
        if os.path.exists(asdmfile):
            if os.path.exists(asdmfile+'/ASDMBinary/'):
                if len(os.listdir(asdmfile+'/ASDMBinary/')) > 0:
                    print 'ASDM file looks good.'
                else:
                    print 'Filling ASDM file.'
                    makeasdm.make(asdmfile, bdfdir)
            else:
                print 'Filling ASDM file.'
                makeasdm.make(asdmfile, bdfdir)
        else:
            print 'Copying ASDM file.'
            try:
                os.system('cp -r %s %s' % (datadir+asdmfile, workdir))   # copy new sdm over to working area
            except:
                raise StandardError('Data file %s not found.' % (datadir+asdmfile))
            makeasdm.make(asdmfile, bdfdir)

    if telcaldir:
# get telcal file
#        timeobj = time.gmtime(os.path.getctime(datadir+asdmfile))   # time of asdmfile creation, assuming we're in archive.
#        year = timeobj.tm_year
#        month = '%02d' % (timeobj.tm_mon)
        mjdref = datetime.datetime(1858, 11, 17, 0, 0, 0, 0, None)   # alternatively, parse file name for mjd
        mjdfile = datetime.timedelta(int((asdmfile).split('.')[-2]))   # assumes format for filenames
        timeobj = mjdref+mjdfile
        year = timeobj.year
        month = '%02d' % timeobj.month
        telcaldirfull = telcaldir+str(year)+'/'+str(month)+'/'     # full path to telcalfil

        if os.path.exists(telcaldirfull+telcalfile):
            if os.path.exists(telcalfile):
                print 'Telcal file already set.'
            else:
                print 'Copying Telcal file.'
                os.system('cp -r %s %s' % (telcaldirfull+telcalfile, workdir))   # copy new sdm over to working area
        else:
            raise StandardError('Telcal file %s not found.' % (telcaldirfull+telcalfile))

    # find scans
    scantree = et.parse(asdmfile + '/Scan.xml')
    scans = [ (int(row.find('scanNumber').text), row.find('sourceName').text, row.find('scanIntent').text) for row in scantree.getiterator('row')]
    # set total number of integrations
    scantree = et.parse(asdmfile + '/Main.xml')
    scanint = [(int(row.find('numIntegration').text),) for row in scantree.getiterator('row')]
    goodscans = n.hstack( (scans, scanint) )
    goodscans = [(int(num), name, int(nints), intent) for (num, name, intent, nints) in goodscans if intentfilter in intent]   # is this correct?
    print 'Found a total of %d scans and %d with given source name(s).' % (len(scans), len(goodscans))
    return goodscans

#def asdm2ms(asdmfile, msfile, scans, dropants='', inttime='0s'):
def asdm2ms(asdmfile, msfile, scans, inttime='0s'):
    """ Converts asdm to ms format for a single scan.
    msfile defines the name template for the ms. Should end in .ms, but "s<scans>" will be put in.
    scans should be a comma-delimited list of scans to include (fed directly to importasdm).
    inttime is string to feed to split command. gives option of integrated data down in time.
    """
# can't get split to work with antenna names...?
#    dropants is a list of antenna names (['ea01', 'ea17']) that will not be kept in the ms.
    # anttree = et.parse(asdmfile + '/Antenna.xml')
    # antlist = [row.find('name').text for row in anttree.getiterator('row')]
    # for drop in dropants.split(','):
    #     if drop in antlist:
    #         print 'Dropping ant %s.' % drop
    #         antlist.pop(antlist.index(drop))
    #     else:
    #         print 'Ant %s not in data.' % drop

    # # make string of ants comma delimited
    # ants = ''
    # for ant in antlist:
    #     ants = ants+ant+','
    # ants = ants[:-1]

    # fill ms file
    msfile2 = msfile[:-3] + '_s' + scans + '.ms'
    if os.path.exists(msfile2):
        print '%s already set.' % msfile2
    else:
        # if filtering ants, use tmp file to split. **NOT TESTED YET**
#        if (antlist != [row.find('name').text for row in anttree.getiterator('row')]) or (inttime != '0s'):
        if inttime != '0s':
            print 'Filtering by int time.'
            importasdm(asdm=asdmfile, vis='tmp_'+msfile2, scans=scans, ocorr_mode='co', savecmds=False, applyflags=True, process_flags=False)
#            split(vis='tmp_'+msfile2, outputvis=msfile2, datacolumn='data', timebin=inttime, antenna=ants)
            split(vis='tmp_'+msfile2, outputvis=msfile2, datacolumn='data', timebin=inttime)
            shutil.rmtree('tmp_'+msfile2)
            shutil.move('tmp_'+msfile2+'.flagversions', msfile2+'.flagversions')
        else:
            importasdm(asdm=asdmfile, vis=msfile2, scans=scans, ocorr_mode='co', corr_mode='co', savecmds=False, applyflags=True, process_flags=False, lazy=True)

    # optionally clean up asdm after checking ms file?
    return msfile2
    
def search(datadir):
    """ Function to ping data directory (as a daemon) and return names of new asdm files.
    """

    pass

def pipeline(asdmfile, workdir='', mode='prep', intentfilter='OBSERVE_TARGET'):
    """ Pipeline to simplify quasi-blind searching over many scans and/or files.
    mode of 'frb' does multi-node processing of frb data.
    mode of 'any' does search over all scans of arbitrary observation.
    """

    if workdir == '':
        workdir = os.chdir(workdir)

    # strip out path
    if '/' in asdmfile:
        asdmfile = asdmfile.split('/')[-1]

    # data filler steps
#   asdmfile = search(datadir)   # daemon process?
    msfile = asdmfile + '.ms'
    telcalfile = asdmfile + '.GN'

    try:
        goodscans = prep(asdmfile, workdir=workdir, intentfilter=intentfilter)    # set up asdm and telcal files, get scans of interest (tuple with scann
        msfile2 = asdm2ms(asdmfile, msfile, goodscans[scan][0])
    except:
        raise

    if mode == 'prep':
        # just filling data...
        return msfile2

    elif mode == 'frb':

        # set good channels
#        nch = 64
#        edgechan = n.round(0.05*nch).astype('int')
#        chans = n.arange(edgechan, nch-edgechan, dtype='int')
        dmarr = [0,19.2033,38.4033,57.6025,76.8036,96.0093,115.222,134.445,153.68,172.93,192.198,211.486,230.797,250.133,269.498,288.894,308.323,327.788,347.292,366.837,386.426,406.062,425.747,445.484,465.276,485.125,505.033,525.005,545.042,565.147,585.322,605.571,625.896,646.3,666.786,687.355,708.012,728.759,749.598,770.532,791.565,812.699,833.936,855.28,876.733,898.299,919.979,941.778,963.697,985.741,1007.91,1030.21,1052.64,1075.21,1097.92,1120.76,1143.76,1166.9,1190.19,1213.63,1237.23,1260.99,1284.92,1309.01,1333.27,1357.7,1382.31,1407.09,1432.06,1457.22,1482.56,1508.1,1533.83,1559.76,1585.89,1612.23,1638.77,1665.53,1692.51,1719.7,1747.12,1774.77,1802.64,1830.75,1859.1,1887.69,1916.53,1945.61,1974.95,2004.54,2034.39,2064.51,2094.9,2125.56,2156.49,2187.71,2219.21,2250.99,2283.07,2315.45,2348.13,2381.12,2414.41,2448.02,2481.94,2516.19,2550.77,2585.68,2620.92,2656.51,2692.44,2728.72,2765.36,2802.36,2839.72,2877.45,2915.55,2954.04,2992.91]
        chans = range(23,31)+range(32,50)+range(70,128)

        d = leanpipe.pipe_thread(filename=msfile, nints=nints, nskip=0, iterint=200, spw=[0,1], chans=chans, dmarr=dmarr, fwhmsurvey=0.5, fwhmfield=0.5, selectpol=['RR','LL'], scan=0, datacol='data', size=25600, res=50, sigma_image=6., searchtype='imageallstat', telcalfile=telcalfile, filtershape=None, savecands=True)

    elif mode == 'any':
        # read structure of data
        os.chdir(workdir)
        ms.open(msfile)
        scans = ms.getscansummary()
        spwinfo = ms.getspectralwindowinfo()
        summary = ms.summary()
        ms.close()

        print 'Searching over %d field(s) in %d scans and %d spws.' % (summary['nfields'], len(scans.keys()), len(spwinfo))

        # search
        scanlist = scans.keys()
        for scan in scanlist:
            nints = int(n.round((scans[scan]['0']['EndTime']-scans[scan]['0']['BeginTime'])*24*3600/scans[scan]['0']['IntegrationTime'])-1)
            spws = scans[scan]['0']['SpwIds']
            for spw in spws:
                nch = spwinfo[str(spw)]['NumChan']
                print
                print 'For scan %d and spw %d, %d integrations and %d channels' % (int(scan), spws[spw], nints, nch)
                edgechan = n.round(0.05*nch).astype('int')
                chans = n.arange(edgechan, nch-edgechan, dtype='int')
                d = leanpipe.pipe_thread(filename=msfile, nints=nints, nskip=0, iterint=min(nints,200), spw=[spws[spw]], chans=chans, dmarr=range(0,300,30), fwhmsurvey=0.5, fwhmfield=0.5, selectpol=['RR','LL'], scan=int(scanlist.index(scan)), datacol='data', size=50000, res=100, sigma_bisp=5., sigma_image=5., calibratedfilter=True, specmodfilter=1.5, searchtype='imageallstat', telcalfile=telcalfile, telcalcalibrator='3C48', filtershape='b')
