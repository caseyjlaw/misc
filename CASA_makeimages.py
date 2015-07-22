"""
Pipeline to make CASA images to build lightcurve from EVLA data.
Copied from EVLA pipeline v3.3.
claw 18sep2012

Run like this:
python> msname='...'
python> execfile('CASA_makeimages.py')
"""

# prepare functions and parameters needed to image
import numpy as np
from EVLA_functions import *
doexe=True
#fieldn = 6
#msname = '12A-339_sb9826039_1.56065.225594965275_f'+str(fieldn)+'c.ms'
bbn = 'bb48'
msname = '12A-339_sb9826039_1.56065.225594965275_'+bbn+'.ms'
nspw = 7      # number of spw in file
fields = [6,7]  # fields to image
start_scan=0 # first scan number to image
stokespar='I'
fluxlogname = 'flux_'+string.join(msname.split('.')[:-1],'.')+'.txt'     # peak flux saved in this file
uvr = '>80klambda'     # uvrange used in imaging

execfile('EVLA_pipe_msinfo_min.py')

def logprint(msg):
    print (msg)
    outlog.write(msg+"\n")
    outlog.flush()

ms.open(msname)
uv_range = ms.range(["uvdist"])
uv_max = uv_range['uvdist'][1]
ms.close()

c=2.997925e8

# set cell size to 1/(3.*Bmax)

#for jj in range(numFields):
for jj in fields:
  for target_scan in field_scans[jj]:
    if target_scan >= start_scan:
      for i in spws_per_scan[str(target_scan)]:
       if i in [2,10,18,28,34,45]:
        if i == 2: ii = '2~9'
        if i == 10: ii = '10~17'
        if i == 18: ii = '18~25'
        if i == 28: ii = '28~33'
        if i == 34: ii = '34~41'
        if i == 45: ii ='45~49'
        ii = '0~'+str(nspw-1)
        print 'field ', jj
        print 'scan ', target_scan
        print 'spw ', ii
        wave=c/reference_frequencies[nspw-1]
        cellsize=206265.*wave/uv_max/3.
        mycell=str(cellsize)+'arcsec'
        fwhm=206265.*wave/25.0
        myimsize=getOptimumSize(fwhm/cellsize)
        imname="field"+str(jj)+'_scan'+str(target_scan)+"_"+bbn
        suffix='.tt0'

        if os.path.exists(imname+'.image'+suffix):
         continue
        else:
         default('clean')
         vis=msname
         imagename=imname
         outlierfile=''
         field=str(jj)
         scan=str(target_scan)
         spw=str(ii)
         selectdata=True
         mode='mfs'
         nterms=2
         reffreq=''
         gridmode=''
         niter=0
         gain=0.1
         threshold='0.0mJy'
         psfmode='clark'
         imagermode=''
         multiscale=[]
         mask='cleanbox'+str(jj)+'.txt'
         interactive=False
         imsize=[myimsize,myimsize]
         cell=[mycell,mycell]
         phasecenter=''
         restfreq=''
         stokes=stokespar
         weighting='natural'
         uvrange=uvr
         uvtaper=False
         modelimage=''
         restoringbeam=['']
         pbcor=False
         minpb=0.2
         calready=False
         allowchunk=False
         async=False
         usescratch=False
         if doexe: clean()

         if (os.path.exists(imname+'.residual'+suffix) & (imstat(imagename=imname+'.residual'+suffix)['max'][0] != 0)):

             default('imstat')
             imagename=imname+'.residual'+suffix
             axes=-1
             region=''
             box=''
             chans=''
             stokes=stokespar
             listit=False
             verbose=False
             tmpout=imstat()

             mythresh=str(3.*tmpout['rms'][0])+'Jy'

             default('clean')
             vis=msname
             imagename=imname
             outlierfile=''
             field=str(jj)
             scan=str(target_scan)
             spw=str(ii)
             selectdata=True
             mode='mfs'
             nterms=2
             reffreq=''
             gridmode=''
             niter=10000
             gain=0.1
             threshold=mythresh
             psfmode='clark'
             imagermode=''
             multiscale=[]
             mask='cleanbox'+str(jj)+'.txt'
             interactive=False
             imsize=[myimsize,myimsize]
             cell=[mycell,mycell]
             phasecenter=''
             restfreq=''
             stokes=stokespar
             weighting='natural'
             uvrange=uvr
             uvtaper=False
             modelimage=''
             restoringbeam=['']
             pbcor=False
             minpb=0.2
             calready=False
             allowchunk=False
             async=False
             usescratch=False
             if doexe: clean()
         
             default('imstat')
             imagename=imname+'.residual'+suffix
             axes=-1
             region=''
             box=''
             chans=''
             stokes=stokespar
             listit=False
             verbose=False
             tmpout=imstat()

             mythresh=str(3.*tmpout['rms'][0])+'Jy'

             default('clean')
             vis=msname
             imagename=imname
             outlierfile=''
             field=str(jj)
             scan=str(target_scan)
             spw=str(ii)
             selectdata=True
             mode='mfs'
             nterms=2
             reffreq=''
             gridmode=''
             niter=10000
             gain=0.1
             threshold=mythresh
             psfmode='clark'
             imagermode=''
             multiscale=[]
             mask='cleanbox'+str(jj)+'.txt'
             interactive=False
             imsize=[myimsize,myimsize]
             cell=[mycell,mycell]
             phasecenter=''
             restfreq=''
             stokes=stokespar
             weighting='natural'
             uvrange=uvr
             uvtaper=False
             modelimage=''
             restoringbeam=['']
             pbcor=False
             minpb=0.2
             calready=False
             allowchunk=False
             async=False
             usescratch=False
             if doexe: clean()

             default('imstat')
             imagename=imname+'.residual'+suffix
             axes=-1
             region=''
             box=''
             chans=''
             stokes=stokespar
             listit=False
             verbose=False
             tmpout=imstat()
             
             finalrms=tmpout['rms'][0]

             # Measure max of image

             default('imstat')
             imagename=imname+'.image'+suffix
             axes=-1
             region=''
             box=''
             chans=''
             stokes=stokespar
             listit=False
             verbose=False
             tmpout=imstat()

             peak=tmpout['max'][0]
             loc = str(tmpout['maxpos'][0]) + ',' + str(tmpout['maxpos'][1])
             outlog=open(fluxlogname,'a')
             logprint(imname + ', peak=' + str(peak) + ', err=' + str(finalrms) + ', at pixel ' + loc)
             outlog.close()
             
         else:
             print("Clean error for Field/spw "+str(jj)+"/"+str(ii))
             continue

  print("Field "+str(jj)+" complete")
