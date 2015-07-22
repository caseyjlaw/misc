import pylab as plt
import pickle, glob, string

def noiseplot(filename,y):
    # give filename enclosed in single quotes
    # y is the column number in the noise data to pull and plot

    # assumes names like:
    #    noise_13B-409_13sep16v1_s61.pkl
    #    cands_13B-409_13sep16v1_s61_dm0-119.pkl

    imagedata = []
    filelist = glob.glob(filename)

    for filename in filelist:
        if y == 3:
            candfilelist = filename[:-4].split('_')[1:]
#            candfilelist.append('dm0-119')
            candfilelist.insert(0, 'cands')
            candfilename = string.join(candfilelist, '_') #+ '.pkl'
            cf = open(glob.glob(candfilename + '*')[0], 'r')
            d = pickle.load(cf)
            npix = d['sizex']/d['res'] * d['sizey']/d['res']
        
        f = open(filename, 'r')
        while True:
            try:
                value=pickle.load(f)
                if y == 3:
                    scale = npix
                else:
                    scale = 1.
                imagedata.append(scale*value[y])

            except EOFError:
                f.close()
                break
   
    if y == 0:
        b='Integration Number'
       
    if y ==1: 
        b = 'Estimate of noise in Raw Data'
        
    if y == 2:
        b = 'Fraction of data flagged'
    
    if y== 3:
        b = 'Estimate of noise in single 5ms image'

    plt.figure(1)
    plt.hist(imagedata, bins=len(imagedata)/10)
    plt.xlabel(b)
    plt.ylabel('Counts')
    plt.show()
