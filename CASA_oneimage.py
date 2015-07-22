"""
Pipeline to make CASA images to build lightcurve from EVLA data.
Copied from EVLA pipeline v3.3.
claw 18sep2012

Run like this:
python> execfile('CASA_makeimages.py')
"""

# prepare functions and parameters needed to image
import numpy as n
from EVLA_functions import *
doexe=True
fieldn = 6  # field number in main obs file
bbn = 'bb48'
#msname = '12A-339_sb9826039_1.56065.225594965275_f'+str(fieldn)+'c.ms'
msname = '12A-339_sb9826039_1.56065.225594965275_'+bbn+'.ms'
prefix = 'test'
#msfieldn = 0  # ms field number to image
msfieldn = fieldn  # ms field number to image
scann = '117~512' # scan to image
spwfreq = 18
spwn = '0~5:5~58' # ms spw number to image
#spwn = spwfreq # ms spw number to image
doclean = 1
antennas = ''
uvr = '>80klambda'
stokesparam = 'I'
psfmoden = 'clarkstokes'

c=2.997925e8
reference_frequencies = n.array([  8.33200000e+09,   8.46000000e+09,   3.74780000e+10,  3.76060000e+10,   3.77340000e+10,   3.78620000e+10, 3.79900000e+10,   3.81180000e+10,   3.82460000e+10,  3.83740000e+10,   2.69680000e+10,   2.70960000e+10,   2.72240000e+10,   2.73520000e+10,   2.74800000e+10,  2.76080000e+10,   2.77360000e+10,   2.78640000e+10,  4.79880000e+10,   4.81160000e+10,   4.82440000e+10,  4.83720000e+10,   4.85000000e+10,   4.86280000e+10,  4.87560000e+10,   4.88840000e+10,   3.90380000e+10, 3.91660000e+10,   3.92940000e+10,   3.94220000e+10,  3.95500000e+10,   3.96780000e+10,   3.98060000e+10,  3.99340000e+10,   2.54880000e+10,   2.56160000e+10,  2.57440000e+10,   2.58720000e+10,   2.60000000e+10,  2.61280000e+10,   2.62560000e+10,   2.63840000e+10,  1.84880000e+10,   1.86160000e+10,   1.87440000e+10,  1.88720000e+10,   1.90000000e+10,   1.91280000e+10,  1.92560000e+10,   1.93840000e+10])

#ms.open(msname)
#uv_range = ms.range(["uvdist"])
#uv_max = uv_range['uvdist'][1]
#ms.close()

uv_max = 6944.8193653207109   # for 12A-339

def getOptimumSize(size):
    '''
    This method takes as input the a size parameter.  The return is the smallest
    integer Y which satisfies the following conditions:
    * Y > size
    * Y = 2^a*3^b*5^c where a,b, and c are non-negative integers and at least one
    of a or b is 0 and c is nonzero
    '''
    def evaluate(pow2, pow3, pow5):
        # Convience method to calculate the value given multiples
        return int(math.pow(2,pow2) *math.pow(3,pow3)*math.pow(5,pow5))
    
    max5 = int(math.ceil(math.log(size,5)))
    returnValue = evaluate(0, 0, max5)
    for pow5 in range(max5,0,-1):
        pow2 = math.ceil(math.log(size/math.pow(5,pow5),2))
        if not pow2 < 0:
            returnValue = min(returnValue, evaluate(pow2,0,pow5))

        pow3 = math.ceil(math.log(size/math.pow(5,pow5),3))
        if not pow3 < 0:
            returnValue = min(returnValue, evaluate(0,pow3,pow5))
    return returnValue

print 'field ', fieldn
print 'scan ', scann
print 'spw ', spwn
print 'bb ', bbn
wave=c/reference_frequencies[spwfreq]
cellsize=206265.*wave/uv_max/3.
mycell=str(cellsize)+'arcsec'
fwhm=206265.*wave/25.0
myimsize=getOptimumSize(fwhm/cellsize)
imname=prefix+str(fieldn)+'_scan'+str(scann)+"_spw"+str(spwn)+"_"+bbn+'_'+str(stokesparam)

if os.path.exists(imname+'.image'):
    print 'Image exists. Continuing clean...'

default('clean')
vis=msname
imagename=imname
outlierfile=''
field=str(msfieldn)
scan=str(scann)
spw=str(spwn)
antenna=antennas
selectdata=True
mode='mfs'
nterms=1
reffreq=''
gridmode=''
niter=0
gain=0.1
threshold='0.0mJy'
psfmode=psfmoden
imagermode=''
multiscale=[]
mask='cleanbox'+str(fieldn)+'.txt'
interactive=False
imsize=[myimsize,myimsize]
cell=[mycell,mycell]
phasecenter=''
restfreq=''
stokes=stokesparam
weighting='natural'
uvrange=uvr
uvtaper=False
modelimage=''
restoringbeam=['']
pbcor=False
minpb=0.2
usescratch=False
calready=False
allowchunk=False
async=False
if doexe: clean()

if (doclean & os.path.exists(imname+'.residual') & (imstat(imagename=imname+'.residual')['max'][0] != 0)):

    default('imstat')
    imagename=imname+'.residual'
    axes=-1
    region=''
    box=''
    chans=''
    stokes=stokesparam
    listit=False
    verbose=False
    tmpout=imstat()

    mythresh=str(3.*tmpout['rms'][0])+'Jy'

    default('clean')
    vis=msname
    imagename=imname
    outlierfile=''
    field=str(msfieldn)
    scan=str(scann)
    spw=str(spwn)
    antenna=antennas
    selectdata=True
    mode='mfs'
    nterms=1
    reffreq=''
    gridmode=''
    niter=10000
    gain=0.1
    threshold=mythresh
    psfmode=psfmoden
    imagermode=''
    multiscale=[]
    mask='cleanbox'+str(fieldn)+'.txt'
    interactive=False
    imsize=[myimsize,myimsize]
    cell=[mycell,mycell]
    phasecenter=''
    restfreq=''
    stokes=stokesparam
    weighting='natural'
    uvrange=uvr
    uvtaper=False
    modelimage=''
    restoringbeam=['']
    pbcor=False
    minpb=0.2
    usescratch=False
    calready=False
    allowchunk=False
    async=False
    if doexe: clean()

    default('imstat')
    imagename=imname+'.residual'
    axes=-1
    region=''
    box=''
    chans=''
    stokes=stokesparam
    listit=False
    verbose=False
    tmpout=imstat()
             
    finalrms=tmpout['rms'][0]

    # Measure max of image
    
    default('imstat')
    imagename=imname+'.image'
    axes=-1
    region=''
    box=''
    chans=''
    stokes=stokesparam
    listit=False
    verbose=False
    tmpout=imstat()

    peak=tmpout['max'][0]
    loc = str(tmpout['maxpos'][0]) + ',' + str(tmpout['maxpos'][1])

    print(imname + ', peak=' + str(peak) + ', err=' + str(finalrms) + ', at pixel ' + loc)

             
