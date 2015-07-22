#!/usr/bin/env python
import os, time, sys
import multiprocessing as mp
import xml.etree.ElementTree as et

# needs to be run by python like so: sbatch --wrap='python sbatch_proc.py asdmfile msfile'
arg0 = sys.argv.index('sbatch_proc.py')
asdmfile = sys.argv[arg0+1]
msfile = sys.argv[arg0+2]

# dm schedule needs to be optimized for memory footprint
#dmthreads = [ ('0,72 ', '17') ,('72,100 ', '14'), ('100,119 ', '10') ]
dmthreads = [ ('0,59 ', '30') ,('59,90 ', '16'), ('90,119 ', '15') ]   # 8-9 GB mem per run (empirically)
#dmthreads = [ ('0,90 ', '30') ,('90,119 ', '30') ]                    # 12 GB mem per run

# processing settings
nodestr = os.environ['SLURM_JOB_NODELIST'][2:]
if '-' in nodestr:  # if it is a range, take range string
    noderanges = os.environ['SLURM_JOB_NODELIST'][3:].split(']')[0].split(',')
else:   # if just one node, make a range string from it
    noderanges = [nodestr + '-' + nodestr]
nodelist = []
for noder in noderanges:
    if '-' in noder:
        for nodenum in range(int(noder.split('-')[0]),int(noder.split('-')[1])+1):
            nodelist.append('cn' + str(nodenum) + ' ')
    else:
        nodelist.append('cn' + str(noder) + ' ')
#nodelist = ['cn' + str(nodenum) +' ' for i in range(len(noderanges)) for nodenum in range(int(noderanges[i].split('-')[0]),int(noderanges[i].split('-')[1])+1)]
print 'Using nodelist', nodelist
proclist = []; lastpart = []; logname = []; nodestate = [0] * len(nodelist)
firstpart = 'srun -N1 -n1 --exclusive -D %s ' % os.getcwd()
middlepart = 'casapy --nogui --nologfile --nologger -c batch_leanpipet.py '
intentfilter = 'OBSERVE_TARGET'
#middlepart = 'sleep 5 '

# get number of scans from ASDM XML file
scantree = et.parse(asdmfile + '/Scan.xml')
totscans = [ (int(row.find('scanNumber').text), row.find('sourceName').text, row.find('scanIntent').text) for row in scantree.getiterator('row')]
goodscans = [(int(num), name, intent) for (num, name, intent) in totscans if intentfilter in intent]
print 'Found a total of %d scans and %d with given source name(s).' % (len(totscans), len(goodscans))
scans = [str(i)+ ' ' for i in range(len(goodscans))]

# build list of tail of processing description that functions as job list
for dmstr,nthreads in dmthreads:
    for scan in scans:
        lastpart.append(scan + dmstr + nthreads)
        logname.append('-o slurm-%j-%N-'+scan[:-1]+'-'+dmstr)

# process all jobs
while len(lastpart) > 0:

    # avoid collisions with ongoing CASA startup files
    if os.path.exists('casapy.log'):
        print 'CASA log present. Sleeping...'
        time.sleep(5)
    else:
        for i in range(len(nodelist)):
            if len(lastpart) == 0:   # escape if no jobs left
                break

            # process, if node is free
            if not nodestate[i]:
#                job = firstpart + '-w ' + str(nodelist[i]) + middlepart
                job = firstpart + logname.pop(0) + '-w ' + str(nodelist[i]) + middlepart + asdmfile + ' ' + msfile + ' ' + lastpart.pop(0)
                proc = mp.Process(target=os.system, args=(job,))
#                try:
                proc.start()
                # except (OSError, EOFError, AttributeError):    # try again, if ipython can't clean up casapy.log file
                #     print 'Caught exception. Waiting to start CASA again.'
                #     time.sleep(4)
                #     try:
                #         proc.start()
                #     except (OSError, EOFError, AttributeError):    # try again, if ipython can't clean up casapy.log file
                #         print 'Caught exception. Waiting to start CASA again.'
                #         time.sleep(4)
                #         proc.start()

                # assuming job started, modify tracking info
                proclist.append( (i, proc) )
                nodestate[i] = 1      # node is working
                print 'Submitted job %s' % job
                print '%d jobs remaining ' % len(lastpart)
                time.sleep(4)

    # look for completed processes and clean them up
    for i,proc in proclist:
        if not proc.is_alive():
            proc.join()
            proclist.remove( (i,proc) )
            nodestate[i] = 0      # node is free
            print 'Job complete.'

    # take a breather
    time.sleep(1)

# clean up processes after all jobs submitted
while len(proclist) > 0:
    for node,proc in proclist:
        proc.join()
        proclist.remove( (node,proc) )
        nodestate[node] = 0      # node is free
        print 'Job complete.'

print 'All Jobs complete.'
