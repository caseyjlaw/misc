#!/usr/bin/env python2.7
#
# split job into nsegment pieces and queue all up with rq
# each job is independent and shares memory. one worker per node.

from rq import Queue, Connection
import os, glob, time, argparse, pickle, string

parser = argparse.ArgumentParser()
parser.add_argument("filename", help="MS filename with full path")
parser.add_argument("--scans", help="scans to search. MS value, not index.", default=-1)
parser.add_argument("--mode", help="'run', 'failed', 'clear'", default='run')
args = parser.parse_args(); filename = args.filename; scans = args.scans

# parameters of search **should be from argv**
workdir = string.join(filename.rstrip('/').split('/')[:-1], '/') + '/'
if workdir == '/':
    workdir = os.getcwd() + '/'
filename = filename.rstrip('/').split('/')[-1]
fileroot = filename.split('_s')[0]
gainfile= workdir + fileroot + '.g2'
bpfile= workdir + fileroot + '.b1'

savenoise = True
savecands = True
nthread=16
nsegments=8
dmarr = [0,19.2033,38.4033,57.6025,76.8036,96.0093,115.222,134.445,153.68,172.93,192.198,211.486,230.797,250.133,269.498,288.894,308.323,327.788,347.292,366.837,386.426,406.062,425.747,445.484,465.276,485.125,505.033,525.005,545.042,565.147,585.322,605.571,625.896,646.3,666.786,687.355,708.012,728.759,749.598,770.532,791.565,812.699,833.936,855.28,876.733,898.299,919.979,941.778,963.697,985.741,1007.91,1030.21,1052.64,1075.21,1097.92,1120.76,1143.76,1166.9,1190.19,1213.63,1237.23,1260.99,1284.92,1309.01,1333.27,1357.7,1382.31,1407.09,1432.06,1457.22,1482.56,1508.1,1533.83,1559.76,1585.89,1612.23,1638.77,1665.53,1692.51,1719.7,1747.12,1774.77,1802.64,1830.75,1859.1,1887.69,1916.53,1945.61,1974.95,2004.54,2034.39,2064.51,2094.9,2125.56,2156.49,2187.71,2219.21,2250.99,2283.07,2315.45,2348.13,2381.12,2414.41,2448.02,2481.94,2516.19,2550.77,2585.68,2620.92,2656.51,2692.44,2728.72,2765.36,2802.36,2839.72,2877.45,2915.55,2954.04,2992.91]
dtarr = [1,2,4,8]   # integer to integrate in time for independent searches
searchtype = 'image1'    # search algorithm: 'image1' is single image snr threshold
sigma_image = 6.0
flagmode='standard'
flagantsol = []
spw = [0,1]
chans = range(6,122)+range(134,250)
uvres = 58   # imaging parameters. set res=size=0 to define from uv coords
npix = 512
scans = [int(i) for i in args.scans.split(',')]
if scans[0] == -1:
    scans = [int(ss) for ss in filename.rstrip('.ms').split('_s')[1].split(',')]  # attempt to extract scans from filename

def main():
    from realtime.RT import pipeline, set_pipeline
    from realtime.parsems import get_metadata

    # queue jobs
    for scan in scans:
        scanind = scans.index(scan)
        print 'Getting metadata for %s, scan %d' % (filename, scan)
        state = get_metadata(workdir+'/'+filename, chans=chans, spw=spw, scan=scanind)
        set_pipeline(state, nthread=nthread, nsegments=nsegments, gainfile=gainfile, bpfile=bpfile, dmarr=dmarr, dtarr=dtarr, savecands=savecands, savenoise=savenoise, sigma_image=sigma_image, flagmode=flagmode, flagantsol=flagantsol, searchtype=searchtype, uvres=uvres, npix=npix)
        print 'Sending %d segments to queue' % (state['nsegments'])
        for segment in range(state['nsegments']):
            q.enqueue_call(func=pipeline, args=(state, segment), timeout=24*3600, result_ttl=24*3600)


def cleanup(filename):
    # cleanup section
    if savecands:
        print 'Aggregating cands files...'
        # define input and outfile files (segmented to single)
        candsroot = 'cands_' + filename.rstrip('.ms') + '_seg*.pkl'
        seglist = glob.glob(candsroot)
        candsfile = 'cands_' + filename.rstrip('.ms') + '.pkl'

        # aggregate cands over segments
        cands = {}
        for cc in seglist:
            pkl = open(cc,'r')
            state = pickle.load(pkl)
            result = pickle.load(pkl)
            for kk in result.keys():
                cands[kk] = result[kk]
            pkl.close()

        # write cands to single file
        pkl = open(candsfile,'w')
        pickle.dump(state, pkl)
        pickle.dump(cands, pkl)
        pkl.close()

        # if successful, delete segment files
        if os.path.exists(candsfile):
            for cc in seglist:
                os.remove(cc)

    if savenoise:
        print 'Aggregating noise files...'
        # define input and outfile files (segmented to single)
        noiseroot = 'noise_' + filename.rstrip('.ms') + '_seg*.pkl'
        seglist = glob.glob(noiseroot)
        noisefile = 'noise_' + filename.rstrip('.ms') + '.pkl'

        # aggregate noise over segments
        noise = []
        for cc in seglist:
            pkl = open(cc,'r')
            result = pickle.load(pkl)
            noise.extend(result)
            pkl.close()

        # write noise to single file
        pkl = open(noisefile,'w')
        pickle.dump(noise, pkl)
        pkl.close()

        # if successful, delete segment files
        if os.path.exists(noisefile):
            for cc in seglist:
                os.remove(cc)


if __name__ == '__main__':
    # connect
    with Connection():

        if args.mode == 'run':
            q = Queue('low')
            main()

        elif args.mode == 'clear':
            q = Queue('high')
            q.empty()
            q = Queue('low')
            q.empty()
            q = Queue('failed')
            q.empty()

        elif args.mode == 'failed':
            q = Queue('failed')
            print 'Failed queue:'
            print q.jobs
            if len(q.jobs):
                print 'First in failed queue:'
                print q.jobs[0].exc_info

        elif args.mode == 'cleanup':
            cleanup(filename)
