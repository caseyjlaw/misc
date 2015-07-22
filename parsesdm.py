#
# functions to read and convert sdm files
# claw, 14jun03
#

import numpy as n
import casautil, tasklib
import os, shutil, subprocess, glob
import xml.etree.ElementTree as et
import bdfparse

qa = casautil.tools.quanta()
me = casautil.tools.measures()

def get_metadata(filename, spw=[], chans=[]):
    """ Parses sdm file to define metadata for observation, including scan info, image grid parameters, pipeline memory usage, etc.
    Mirrors parsems.get_metadata().
    """

    # create primary state dictionary
    d = {}
    ss = filename.rstrip('/').split('/')
    if len(ss) > 1:
        d['filename'] = ss[-1]
        d['workdir'] = string.join(ss[:-1], '/') + '/'
        os.chdir(d['workdir'])
    else:
        d['filename'] = filename
        d['workdir'] = os.getcwd() + '/'


#    d['ants']
#    d['nants'] = len(d['ants'])
#    d['blarr']
#    d['nbl'] = len(d['blarr'])
#    d['scans']
#    d['urange'] = urange
#    d['vrange'] = vrange

    root = et.parse(d['workdir'] + d['filename'] + '/SpectralWindow.xml').getroot()
    reffreq = [float(row.find('chanFreqStart').text) for row in root.iter('row')]
    numchan = [float(row.find('numChan').text) for row in root.iter('row')]
    chansize = [float(row.find('chanFreqStep').text) for row in root.iter('row')]
    spectralwindow = [int(row.find('spectralWindowId').text.split('_')[1]) for row in root.iter('row')]

    if len(spw):
        d['spw'] = sorted(spectralwindow)[spw]
    else:
        d['spw'] = sorted(spectralwindow)
    spwch = []
    for freq in sorted(reffreq):
        ii = reffreq.index(freq)
        if spectralwindow[ii] in d['spw']:
            spwch.extend(list(n.linspace(reffreq[ii], reffreq[ii]+numchan[ii]*chansize[ii], numchan[ii])))
    d['freq_orig'] = n.array(spwch)/1e9
           
    d['nspw'] = len(d['spw'])
    if len(chans):
        d['freq'] = d['freq_orig'][chans]
        d['chans'] = chans
    else:
        d['freq'] = d['freq_orig']
        d['chans'] = range(len(d['freq']))
    d['nchan'] = len(d['freq'])

    return d

def read_bdf(d, scannum, readints=0, nskip=0):
    """ Uses Peter's bdfparse to get bdf data without CASA and load into numpy array.
    scannum is a scan number (1-based)
    """

    scans = filter_scans(d['workdir'] + d['filename'])   # find bdf from filtered list
    print 'Reading scan %d (index %d) of source %s.' % (scans[scannum][0], scannum, scans[scannum][1])

    bdfname = d['workdir'] + d['filename'] + '/ASDMBinary/*' + str(scans[scannum][3])  # bdf number is unique
    bdflist = glob.glob(bdfname)
    if len(bdflist) != 1:
        raise IndexError('Scan %s has no bdf.' % bdfname)
    fp = open(bdflist[0])
    bdf = bdfparse.BDFData(fp).parse()

    # read integrations
    if readints == 0:
        readints = bdf.n_integrations
    data_arr = n.empty( (readints, bdf.n_baselines, bdf.n_channels, bdf.n_basebands), dtype='complex64', order='C')
    for i in xrange(nskip, nskip+readints):
        data_arr[i] = bdf.get_data ('crossData.bin', i)

    fp.close()
    return data_arr

def get_uvw(d, scannum):
    """ Read antenna positions, calculate and return uvw in meters.
    scannum is scan number (1-based).
    """

    scans, sources = readscans(d['workdir'] + d['filename'])
    scan = scans[scannum]
    print 'Calculating uvw for scan %d of source %s' % (scannum, scan['source'])

    # define metadata for uvw calculation
    me.doframe(me.observatory('vla'))
    datetime = scan['start']
    me.doframe(me.epoch('utc', datetime))

    sourcenum = [kk for kk in sources.keys() if sources[kk]['source'] == scans[scannum]['source']][0]
    direction = me.direction('J2000', str(n.degrees(sources[sourcenum]['ra']))+'deg', str(n.degrees(sources[sourcenum]['dec']))+'deg')
    me.doframe(direction)

    print datetime, direction

    root = et.parse(d['workdir'] + d['filename'] + '/Station.xml').getroot()
    positions = [rr.find('position').text.split(' ') for rr in root.iter('row') if 'ANTENNA' in rr.find('type').text]
    x = [float(positions[i][2]) for i in range(len(positions))]
    y = [float(positions[i][3]) for i in range(len(positions))]
    z = [float(positions[i][4]) for i in range(len(positions))]
    ants = me.position('itrf', qa.quantity(x, 'm'), qa.quantity(y, 'm'), qa.quantity(z, 'm'))
    bls=me.asbaseline(ants)
    uvwlist = me.expand(me.touvw(bls)[0])[1]['value']
    u = uvwlist[0::3] * d['freq_orig'][0] * (1e9/3e8) * (-1)     # -1 to keep consistent with ms reading convention
    v = uvwlist[1::3] * d['freq_orig'][0] * (1e9/3e8) * (-1)
    w = uvwlist[2::3] * d['freq_orig'][0] * (1e9/3e8) * (-1)

    return u.astype('float32'), v.astype('float32'), w.astype('float32')

def filter_scans(sdmfile, namefilter='', intentfilter=''):
    """ Parses xml in sdmfile to get scan info for those containing 'namefilter' and 'intentfilter'
    """

    goodscans = {}
    # find scans
    scantree = et.parse(sdmfile + '/Scan.xml')
    scans = [ (int(row.find('scanNumber').text), row.find('sourceName').text, row.find('scanIntent').text) for row in scantree.getiterator('row')]
    # set total number of integrations
    root = et.parse(sdmfile + '/Main.xml').getroot()
    scanint = [int(row.find('numIntegration').text) for row in root.iter('row')]
    bdfnum = [int(rr.attrib['entityId'].split('/')[-1]) for rr in root.iter('EntityRef')]
    for i in range(len(scans)):
        if intentfilter in scans[i][2]:
            if namefilter in scans[i][1]:
                goodscans[scans[i][0]] = (scans[i][1], scanint[i], bdfnum[i])
#    goodscans = [(int(num), name, nints, bdfnum) for (num, name, intent, nints, bdfnum) in scans2 if intentfilter in intent if namefilter in name]
    print 'Found a total of %d scans and %d with name=%s and intent=%s.' % (len(scans), len(goodscans), namefilter, intentfilter)
    return goodscans

def sdm2ms(sdmfile, msfile, scan, inttime='0'):
    """ Converts sdm to ms format for a single scan.
    msfile defines the name template for the ms. Should end in .ms, but "s<scan>" will be put in.
    scan is string of (sdm counted) scan number.
    inttime is string to feed to split command. gives option of integrated data down in time.
    """

    # fill ms file
    msfile2 = msfile[:-3] + '_s' + scan + '.ms'
    if os.path.exists(msfile2):
        print '%s already set.' % msfile2
    else:
        # if filtering ants, use tmp file to split. **NOT TESTED YET**
#        if (antlist != [row.find('name').text for row in anttree.getiterator('row')]) or (inttime != '0s'):
        if inttime != '0':
            print 'Filtering by int time.'
            subprocess.call(['asdm2MS', '--ocm co --icm co --lazy --scans', scan, sdmfile, 'tmp_'+msfile2])
            cfg = tasklib.SplitConfig()  # configure split
            cfg.vis = 'tmp_'+msfile2
            cfg.out = msfile2
            cfg.timebin=inttime
            cfg.col = 'data'
            cfg.antenna='*&*'  # discard autos
            tasklib.split(cfg)  # run task
            # clean up
            shutil.rmtree('tmp_'+msfile2)
        else:
            subprocess.call(['asdm2MS', '--ocm co --icm co --lazy --scans', scan, sdmfile, msfile2])

    return msfile2

####################
# Functions below taken from "readscans.py" by Adam Deller, Steve Myers, and others...
####################
# Usage:
#         from readscans import *
#         sdmfile = '10B-148.sb2269046.eb2519548.55530.44014204861'
#         myscans, mysources = readscans(sdmfile)
#    then
#         listscans(myscans, mysources)
#    To show receiver bands used:
#         readrx(sdmfile)
#
#
# Notes: readflags
#    Returns a dictionary with the parameters
#    from the SDM Scan.xml tables indexed by scan_no: 
#       myscans[scan_no]['start']        start date/time (string)
#       myscans[scan_no]['end']          end date/time (string)
#       myscans[scan_no]['timerange']    start~end date/time (string)
#       myscans[scan_no]['source']       source name (string)
#       myscans[scan_no]['intent']       scan intent(s) (string)
#       myscans[scan_no]['nsubs']        number of subscans (int)

def call_qatime(arg, form='', prec=0):
    """
    This is a wrapper for qa.time(), which in casa 4.0 returns a list of 
    strings instead of just a scalar string.  In this case, return the first 
    value in the list.
    - Todd Hunter
    """

    result = qa.time(arg, form=form, prec=prec)
    if (type(result) == list or type(result)==n.ndarray):
        return(result[0])
    else:
        return(result)

def readscans(sdmfile):
    if (os.path.exists(sdmfile) == False):
        print "Could not find the SDM file = ", sdmfile
        return([],[])
    if (os.path.exists(sdmfile+'/Scan.xml') == False):
        print "Could not find the Scan.xml file.  Are you sure this is an SDM?"
        return([],[])
        
    try:
        from xml.dom import minidom
    except ImportError, e:
        print "failed to load xml.dom.minidom:\n", e
        exit(1)

    # read Scan.xml into dictionary also and make a list
    xmlscans = minidom.parse(sdmfile+'/Scan.xml')
    scandict = {}
    rowlist = xmlscans.getElementsByTagName("row")
    for rownode in rowlist:
        rowfid = rownode.getElementsByTagName("scanNumber")
        fid = int(rowfid[0].childNodes[0].nodeValue)
        # number of subscans
        try:
            # ALMA
            rowsubs = rownode.getElementsByTagName("numSubScan")
            nsubs = int(rowsubs[0].childNodes[0].nodeValue)
        except:
            # EVLA
            rowsubs = rownode.getElementsByTagName("numSubscan")
            nsubs = int(rowsubs[0].childNodes[0].nodeValue)
        # intents
        rownint = rownode.getElementsByTagName("numIntent")
        nint = int(rownint[0].childNodes[0].nodeValue)

        rowintents = rownode.getElementsByTagName("scanIntent")
        sint = str(rowintents[0].childNodes[0].nodeValue)
        sints = sint.split()
        rint = ''
        for r in range(nint):
            intent = sints[2+r]
            if rint=='':
                rint = intent
            else:
                rint += ' '+intent

        # start and end times in mjd ns
        rowstart = rownode.getElementsByTagName("startTime")
        start = int(rowstart[0].childNodes[0].nodeValue)
        startmjd = float(start)*1.0E-9/86400.0
        t = qa.quantity(startmjd,'d')
        starttime = call_qatime(t,form="ymd",prec=8)
        rowend = rownode.getElementsByTagName("endTime")
        end = int(rowend[0].childNodes[0].nodeValue)
        endmjd = float(end)*1.0E-9/86400.0
        t = qa.quantity(endmjd,'d')
        endtime = call_qatime(t,form="ymd",prec=8)
        # source name
        rowsrc = rownode.getElementsByTagName("sourceName")
        if (len(rowsrc) < 1):
            print "Scan %d appears to be corrupt." % (len(scandict)+1)
        else:
            src = str(rowsrc[0].childNodes[0].nodeValue)
            # to find out what all is available,
#            print rownode.getElementsByTagName("*")
            scandict[fid] = {}
            scandict[fid]['start'] = starttime
            scandict[fid]['startmjd'] = startmjd
            scandict[fid]['end'] = endtime
            scandict[fid]['endmjd'] = endmjd
#            print "starttime = ", starttime
#            print "endtime = ", endtime
            timestr = starttime+'~'+endtime
            scandict[fid]['timerange'] = timestr
            scandict[fid]['source'] = src
            scandict[fid]['intent'] = rint
            scandict[fid]['nsubs'] = nsubs
            scandict[fid]['duration'] = endmjd-startmjd
#    print '  Found ',rowlist.length,' scans in Scan.xml'

    # read Source.xml into dictionary also and make a list
    xmlsources = minidom.parse(sdmfile+'/Source.xml')
    sourcedict = {}
    sourcelist = []
    sourceId = []
    rowlist = xmlsources.getElementsByTagName("row")
    for rownode in rowlist:
        rowfid = rownode.getElementsByTagName("sourceId")
        fid = int(rowfid[0].childNodes[0].nodeValue)

        # source name
        rowsrc = rownode.getElementsByTagName("sourceName")
        src = str(rowsrc[0].childNodes[0].nodeValue)
        try:
            rowsrc = rownode.getElementsByTagName("directionCode")
            directionCode = str(rowsrc[0].childNodes[0].nodeValue)
        except:
            directionCode = ''
        rowsrc = rownode.getElementsByTagName("direction")
        (ra,dec) = rowsrc[0].childNodes[0].nodeValue.split()[2:4]
        ra = float(ra)
        dec = float(dec)
        if (src not in sourcelist):
            sourcelist.append(src)
            sourceId.append(fid)
            sourcedict[fid] = {}
#            sourcedict[fid]['sourceName'] = src
            sourcedict[fid]['source'] = src
            sourcedict[fid]['directionCode'] = directionCode
            sourcedict[fid]['ra'] = ra
            sourcedict[fid]['dec'] = dec
#            print "Loading source %s to index %d" % (src,fid)
        else:
            ai = sourceId[sourcelist.index(src)]
#            print "Source %s is already at index %d = ID:%d" % (src,sourcelist.index(src),ai)
            if (ra != sourcedict[ai]['ra'] or dec != sourcedict[ai]['dec']):
                print "WARNING: Multiple directions found for source %d = %s" % (fid,src)
                ras = (ra - sourcedict[ai]['ra'])*180*3600*math.cos(dec)/math.pi
                decs = (dec - sourcedict[ai]['dec'])*180*3600/math.pi
                print "The difference is (%f,%f) arcseconds." % (ras,decs)
#    for src in range(len(sourcedict)):
#        print "%s direction = %f, %f" % (sourcedict[src]['sourceName'],
#                                         sourcedict[src]['ra'],
#                                         sourcedict[src]['dec'])
        
    # return the dictionary for later use
    return [scandict, sourcedict]
# Done readscans

def listscans(dicts):
    myscans = dicts[0]
    mysources = dicts[1]
    if (myscans == []): return
    # Loop over scans
    for key in myscans.keys():
        mys = myscans[key]
        src = mys['source']
        tim = mys['timerange']
        sint= mys['intent']
        dur = mys['duration']*1440
        print '%8i %24s %48s  %.1f minutes  %s ' % (key, src, tim, dur, sint)
    durations = duration(myscans)
    print '  Found ', len(mysources),' sources in Source.xml'
    for key in durations:
        for mysrc in mysources.keys():
#            if (key[0] == mysources[mysrc]['sourceName']):
            if (key[0] == mysources[mysrc]['source']):
                ra = mysources[mysrc]['ra']
                dec = mysources[mysrc]['dec']
                directionCode = mysources[mysrc]['directionCode']
                break
        raString = qa.formxxx('%.12frad'%ra,format('hms'))
        decString = qa.formxxx('%.12frad'%dec,format('dms')).replace('.',':',2)
        print '   Total %24s (%d)  %5.1f minutes  (%.3f, %+.3f radian) %s: %s %s' % (key[0], int(mysrc), key[1], ra, dec, directionCode, raString, decString)
    durations = duration(myscans,nocal=True)
    for key in durations:
        print '   Total %24s      %5.1f minutes (neglecting pntg, atm & sideband cal. scans)' % (key[0],key[1])
    return
# Done

def duration(myscans, nocal=False):
    durations = []
    for key in myscans.keys():
        mys = myscans[key]
        src = mys['source']
        if (nocal and (mys['intent'].find('CALIBRATE_SIDEBAND')>=0 or
                       mys['intent'].find('CALIBRATE_POINTING')>=0 or
                       mys['intent'].find('CALIBRATE_ATMOSPHERE')>=0)):
            dur = 0
        else:
            dur = mys['duration']*1440
        new = 1
        for s in range(len(durations)):
            if (src == durations[s][0]):
                new = 0
                source = s
        if (new == 1):
            durations.append([src,dur])
        else:
            durations[source][1] = durations[source][1] + dur
    return(durations)
    
def readrx(sdmfile):

    # read Scan.xml into dictionary also and make a list
    xmlrx = minidom.parse(sdmfile+'/Receiver.xml')
    rxdict = {}
    rxlist = []
    rowlist = xmlrx.getElementsByTagName("row")
    for rownode in rowlist:
        a = rownode.getElementsByTagName("*")
        rowrxid = rownode.getElementsByTagName("receiverId")
        rxid = int(rowrxid[0].childNodes[0].nodeValue)
        rowfreqband = rownode.getElementsByTagName("frequencyBand")
        freqband = str(rowfreqband[0].childNodes[0].nodeValue)
        print "rxid = %d, freqband = %s" % (rxid,freqband)
    # return the dictionary for later use
    return rxdict
