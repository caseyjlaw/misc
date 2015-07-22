#!/usr/bin/env python
# claw, 13aug19
# 
# Functions to take skeletal asdm file (as found on CBE) and fill with links to bdfs.
# Allows import of data via CASA MS filler.

import os, re, sys, string

def bdfnums(asdmpath):
    """ Takes path to an ASDM Main.xml file as input and gets all bdf numbers as a list.
    replaces 'grep bdf /home/mchammer/evla/mcaf/workspace/C_osro8_000.56516.855799363424/Main.xml'
    """

    file = open(asdmpath+'/Main.xml', "r")

    bdfnums = []
    for line in file:
        if re.search('bdf', line):
            bdfnum = line.split('///evla/bdf/')[1][:-4]   # hardwired based on one example!
            bdfnums.append(bdfnum)

    return bdfnums

def findbdfs(bdfnums, bdfdir='/lustre/evla/wcbe/data/bunker/'):
    """ Uses list of bdf numbers to make list of bdf files.
    """

    directory = os.listdir(bdfdir)

    bdflist = []
    for bdfnum in bdfnums:
        for line in directory:
            if re.search(bdfnum, line):
                bdflist.append(bdfdir + line)

    return bdflist

def linkbdfs(asdm, bdflist):
    """ Takes list of paths to bdfs and places them inside asdm as symlinks.
    asdm should be relative path to asdm directory that will be filled
    """

    if 'mchammer' in asdm:
        raise StandardError('You are trying to edit original data. Bad!')

    os.chdir(asdm)
    try:
        os.mkdir('ASDMBinary')
    except:
        pass
    os.chdir('ASDMBinary')
    for bdf in bdflist:
        os.system('ln -s ' + bdf + ' .')

def make(asdmfile, bdfdir='/lustre/evla/wcbe/data/bunker/'):
    """ Function to do all of the above.
    Argument is path to ASDM directory (with Main.xml file)
    """

    nums = bdfnums(asdmfile)  # get bdf numbers
    if len(nums) == 0:
        raise StandardError('No BDF numbers found')
    bdfs = findbdfs(nums, bdfdir=bdfdir)    # get path to bdfs
    if len(bdfs) == 0:
        raise StandardError('No BDF files found')
    linkbdfs(asdmfile, bdfs)   # symlink into asdm file

