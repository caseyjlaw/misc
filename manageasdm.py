#
# functions to read and convert asdm files
# claw, 14jun03
#

import numpy as n
import casautil, tasklib
import os, shutil, subprocess
import xml.etree.ElementTree as et

def filterscans(asdmfile, namefilter='', intentfilter=''):
    """ Parses xml in asdmfile to get scan info for those containing 'namefilter' and 'intentfilter'
    """

    # find scans
    scantree = et.parse(asdmfile + '/Scan.xml')
    scans = [ (int(row.find('scanNumber').text), row.find('sourceName').text, row.find('scanIntent').text) for row in scantree.getiterator('row')]
    # set total number of integrations
    scantree = et.parse(asdmfile + '/Main.xml')
    scanint = [(int(row.find('numIntegration').text),) for row in scantree.getiterator('row')]
    goodscans = n.hstack( (scans, scanint) )
    goodscans = [(int(num), name, int(nints), intent) for (num, name, intent, nints) in goodscans if intentfilter in intent if namefilter in name]
    print 'Found a total of %d scans and %d with name=%s and intent=%s.' % (len(scans), len(goodscans), namefilter, intentfilter)
    return goodscans

def asdm2ms(asdmfile, msfile, scan, inttime='0'):
    """ Converts asdm to ms format for a single scan.
    msfile defines the name template for the ms. Should end in .ms, but "s<scan>" will be put in.
    scan is string of (asdm counted) scan number.
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
            subprocess.call(['asdm2MS', '--ocm co --icm co --lazy --scans', scan, asdmfile, 'tmp_'+msfile2])
#           importasdm(asdm=asdmfile, vis='tmp_'+msfile2, scans=scans, ocorr_mode='co', savecmds=False, applyflags=True, process_flags=False) # replaced!
#            split(vis='tmp_'+msfile2, outputvis=msfile2, datacolumn='data', timebin=inttime, antenna=ants)
#            split(vis='tmp_'+msfile2, outputvis=msfile2, datacolumn='data', timebin=inttime)
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
            subprocess.call(['asdm2MS', '--ocm co --icm co --lazy --scans', scan, asdmfile, msfile2])
#           importasdm(asdm=asdmfile, vis=msfile2, scans=scans, ocorr_mode='co', corr_mode='co', savecmds=False, applyflags=True, process_flags=False, lazy=True) # replaced!

    return msfile2
