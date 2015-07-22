import cPickle
import pylab as p

def parse(pklfile):
    pkl = open(pklfile,'r')
    d = cPickle.load(pkl)
    loc, prop = cPickle.load(pkl)
    ints = [loc[i][2] for i in range(len(loc))]
    dms = [loc[i][3] for i in range(len(loc))]
    snrs = [prop[i][0] for i in range(len(prop))]
    imgs = [prop[i][1] for i in range(len(prop))]
    specs = [prop[i][2].data.real for i in range(len(prop))]

    return ints, dms, snrs, imgs, specs

def plotdmt(ints, dms, snrs):
    p.figure(1)
    p.scatter(ints, dms, s=snrs, facecolor='none', alpha=0.5)
    p.show()
