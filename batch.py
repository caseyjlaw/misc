#!/usr/bin/env python
import os, time
import multiprocessing as mp

# processing settings
nodestr = os.environ['SLURM_JOB_NODELIST'][2:]
if '-' in nodestr:  # if it is a range, take range string
    noderanges = os.environ['SLURM_JOB_NODELIST'][2:].split('[')[1].split(']')[0].split(',')
else:   # if just one node, make a range string from it
    noderanges = [nodestr + '-' + nodestr]
nodelist = ['cn' + str(nodenum) +' ' for i in range(len(noderanges)) for nodenum in range(int(noderanges[i].split('-')[0]),int(noderanges[i].split('-')[1])+1)]
proclist = []; lastpart = []; nodestate = [0] * len(nodelist)
firstpart = 'srun -N1 -n1 --exclusive -D /home/caseyl/projects -o slurm-%j-%N '
middlepart = 'casapy --nogui --nologfile --nologger -c batch_leanpipet.py '
#middlepart = 'sleep 5 '
asdmfile = '13B-409_13sep17.asdm '
msfile = '13B-409_13sep17v1.ms '
#scans = ['0 ','1 ','2 ','3 ']
#dms = [ '0,31 ', '31,59 ', '59,83 ', '83,100 ', '100,111 ', '111,119' ]
scans = [str(i)+' ' for i in range(26)]
dmthreads = [ ('0,72 ', '17') ,('72,100 ', '14'), ('100,119 ', '10') ]

# build list of tail of processing description that functions as job list
for scan in scans:
    for dmstr,nthreads in dmthreads:
        lastpart.append(scan + dmstr + nthreads)

# process all jobs
while len(lastpart) > 0:

    # avoid collisions with ongoing CASA startup files
    if os.path.exists('casapy.log'):
        print 'CASA log present. Sleeping...'
        time.sleep(4)
    else:
        for i in range(len(nodelist)):
            if len(lastpart) == 0:   # escape if no jobs left
                break

            # process, if node is free
            if not nodestate[i]:
#                print lastpart.pop(0)
#                job = firstpart + '-w ' + str(nodelist[i]) + middlepart
                job = firstpart + '-w ' + str(nodelist[i]) + middlepart + asdmfile + msfile + lastpart.pop(0)
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
