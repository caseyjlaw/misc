# claw, 15 May 2012
#
# CASA script to do common calibration steps
# run from within casapy with 'execfile()'.

# 23may2012: oh shit! cal coordinate is off by 6m, which is > 1 primary beam... use fluxcal to do phase cal, too.

import string

# data file names
masterms = '12A-336_sb9667618_1b.56040.87127945602.ms'  # second good try.
#'12A-336_sb9653146_1.56034.879022291665.ms'  # first good try. no flux cal
targetms = '12A-336_1b_j0628.ms'
calms = '12A-336_1b_j0628cal.ms'
fluxcalms = '12A-336_1b_3c147.ms'     # this observation got no flux cal
refant = 'ea21'

fluxcal = '0542+498=3C147'
cal = 'J0632+103'
target = 'J0628+09'

# calibration product file names
bpms = string.join(fluxcalms.split('.')[:-1], '.') + '_bp.ms'
gain0ms = string.join(fluxcalms.split('.')[:-1], '.') + '_gain0.ms'
gain1ms = string.join(fluxcalms.split('.')[:-1], '.') + '_gain1.ms'
gain2ms = string.join(calms.split('.')[:-1], '.') + '_gain2.ms'

# prep for rerun
clearstat()
clearcal(vis=masterms)   # if repeating this script, this step is probably a good idea

# manually check on field numbers. comment this out afterwards.
#ms.open(masterms)
#ms.summary()

# flag data
#flagautocorr(vis=masterms)
#flagdata(vis=masterms, mode='quack', flagbackup=True, quackinterval=15, quackmode='beg')
#flagdata(vis=masterms, mode='manualflag', clipminmax=[-1e-10,1e-10], clipoutside=False, clipexpr='ABS RR')  # flag zeros
#flagdata(vis=masterms, mode='manualflag', clipminmax=[-1e-10,1e-10], clipoutside=False, clipexpr='ABS LL')  # flag zeros
#flagdata(vis=masterms, mode='manualflag', scan='10')   # flag first 3c147 scan that has wonky calibration result
#flagdata(vis=masterms, mode='manualflag', antenna='ea19,ea23', correlation='LL')
#flagdata(vis=masterms, mode='manualflag', antenna='ea13', correlation='RR')
# could also flag highest of rfi...
#flagdata(vis=masterms, mode='manualflag', clipminmax=[1e-10,1e-1?], clipoutside=True, clipexpr='ABS RR')  # flag high points
#flagdata(vis=masterms, mode='manualflag', clipminmax=[1e-10,1e-1?], clipoutside=True, clipexpr='ABS LL')  # flag high points

# flux cal
setjy(vis=masterms, modimage='3C147_L.im', standard='Perley-Butler 2010', fluxdensity=-1, field=fluxcal)

# phase cal before bp cal (prevent wrapping)
gaincal(vis=masterms, field=fluxcal, caltable=gain0ms, refant=refant, spw='0:25~36', calmode='p', solint='inf', combine='', minsnr=5)  # one solution per scan if combine=''

# visualize calibration
plotcal(caltable=gain0ms, field=fluxcal, xaxis='time', yaxis='phase', poln='R', figfile='plotcal-G0-phase-R.png')
plotcal(caltable=gain0ms, field=fluxcal, xaxis='time', yaxis='phase', poln='L', figfile='plotcal-G0-phase-L.png')

# bandpass calibrate
bandpass(vis=masterms, field=fluxcal, caltable=bpms, spw='0:5~58', refant=refant, solnorm=True, combine='scan', solint='inf', bandtype='BPOLY', degamp=8, gaintable=[gain0ms])

# visualize calibration
plotcal(caltable=bpms, field=fluxcal, xaxis='chan', yaxis='amp', poln='R', figfile='plotcal-B0-amp-R.png')
plotcal(caltable=bpms, field=fluxcal, xaxis='chan', yaxis='amp', poln='L', figfile='plotcal-B0-amp-L.png')
plotcal(caltable=bpms, field=fluxcal, xaxis='chan', yaxis='phase', poln='R', figfile='plotcal-B0-phase-R.png')
plotcal(caltable=bpms, field=fluxcal, xaxis='chan', yaxis='phase', poln='L', figfile='plotcal-B0-phase-L.png')

# gain calibrate
gaincal(vis=masterms, field=fluxcal, caltable=gain1ms, spw='0:5~58', solint='inf', refant=refant, gaintype='G', calmode='ap', solnorm=False, gaintable=[bpms])   # calibrate flux calibrator first
#gaincal(vis=masterms, field=fluxcal, caltable=gain1ms, spw='0~1:5~58', solint='inf', refant=refant, gaintype='G', calmode='ap', solnorm=False, gaintable=[bpms], combine='spw')   # calibrate flux calibrator first

# transfer flux scale to phase calibrator for new cal table
#fluxscale(vis=masterms, caltable=gain1ms, fluxtable=gain2ms, reference=[fluxcal], transfer=[cal])

# visualize calibration
plotcal(caltable=gain1ms, xaxis='time', yaxis='phase', poln='R', figfile='plotcal-G2-phase-R.png')
plotcal(caltable=gain1ms, xaxis='time', yaxis='phase', poln='L', figfile='plotcal-G2-phase-L.png')
plotcal(caltable=gain1ms, xaxis='time', yaxis='amp', poln='R', figfile='plotcal-G2-amp-R.png')
plotcal(caltable=gain1ms, xaxis='time', yaxis='amp', poln='L', figfile='plotcal-G2-amp-L.png')

# apply calibration.
applycal(vis=masterms, gaintable=[gain1ms, bpms], field=fluxcal, spw='0', interp=['nearest',''], gainfield=[fluxcal, ''], parang=False, calwt=False)
applycal(vis=masterms, gaintable=[gain1ms, bpms], field=target, spw='0', interp=['nearest', ''], gainfield=[fluxcal, ''], parang=False, calwt=False)  # using nearest, since fluxcal does not span targets scans. probably fine, since fluxcal calibration is stable in time.
applycal(vis=masterms, gaintable=[gain1ms, bpms], field=cal, spw='0', interp=['nearest', ''], gainfield=[fluxcal, ''], parang=False, calwt=False)  # using nearest, since fluxcal does not span targets scans. probably fine, since fluxcal calibration is stable in time.
#applycal(vis=masterms, gaintable=[gain2ms, bpms], field=cal, spw='0', interp=['nearest', ''], gainfield=[fluxcal, ''], parang=False, calwt=False, spwmap=[[0,0],[0,1]])  # using nearest, since fluxcal does not span targets scans. probably fine, since fluxcal calibration is stable in time.

split(vis=masterms, outputvis=calms, field=0, spw='0', uvrange='>0 lambda', datacolumn='corrected')
split(vis=masterms, outputvis=targetms, field=1, spw='0', uvrange='>0 lambda', datacolumn='corrected')
split(vis=masterms, outputvis=fluxcalms, field=2, spw='0', uvrange='>0 lambda', datacolumn='corrected')
#split(vis=masterms, outputvis=targetms, field=1, spw='0~1', uvrange='>0 lambda', datacolumn='corrected')
#split(vis=masterms, outputvis=fluxcalms, field=2, spw='0~1', uvrange='>0 lambda', datacolumn='corrected')

# test images
#clean(vis=fluxcalms, imagename='img_3c147_sb0_I', mode='mfs', niter=1000, threshold='3mJy', imsize=600, cell='2.5arcsec', stokes='I', weighting='natural', interactive=False, imagermode='')
#clean(vis=calms, imagename='img_j0628cal_sb0_I', mode='mfs', niter=1000, threshold='3mJy', imsize=600, cell='2.5arcsec', stokes='I', weighting='natural', interactive=False, imagermode='')
#clean(vis=targetms, imagename='img_j0628_sb0_I', mode='mfs', niter=1000, threshold='3mJy', imsize=600, cell='2.5arcsec', stokes='I', weighting='natural', interactive=False, imagermode='')
