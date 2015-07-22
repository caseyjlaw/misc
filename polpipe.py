# quick and dirty script (log) for redoing sgra* calibration with polarimetry
import numpy as n

bbn = 'bb39'
msname = '12A-339_sb9826039_1.56065.225594965275_' + bbn + '.ms'
calname = 'cal7'

if bbn == 'bb18':
    bandname = 'K'
    nspw = 5
    i0=2.9324 # from setjy for 3c286 (field=2)
    p0=i0*0.1257 # from polcal page
    reffreqname = '18872MHz'  # at spw=45
    refantn = 'ea19'
elif bbn == 'bb25':
    bandname = 'K'
    nspw = 7
    i0=2.3379
    p0=i0*0.1257
    reffreqname = '25488MHz'   # at spw=34
#    refantn = 'ea21'
    refantn = 'ea19'
elif bbn == 'bb27':
    bandname = 'K'
    nspw = 8
    i0=2.2392
    p0=i0*0.1257
    reffreqname = '26968MHz'    # at spw=10
#    refantn = 'ea21'
    refantn = 'ea19'
elif bbn == 'bb37':
    bandname = 'K'
    nspw = 8
    i0=1.7359
    p0=i0*0.1257
    reffreqname = '37478MHz'     # at spw=2
#    refantn = 'ea21'
    refantn = 'ea19'
elif bbn == 'bb39':
    bandname = 'K'
    nspw = 6
    i0=1.6729
    p0=i0*0.1257
    reffreqname = '39294MHz'     # at spw=28
    refantn = 'ea19'
elif bbn == 'bb48':
    bandname = 'Q'
    nspw = 7
    i0=1.4295
    p0=i0*0.1257
    reffreqname = '47988MHz'     # at spw=18
#    refantn = 'ea05'
    refantn = 'ea19'
        
# old CASA
#import Script_plotWX
#tau = Script_plotWX.plotWX(vis=msname)
# middle CASA
#tau=plotweather(vis=msname,seasonal_weight=0.5,doPlot=True)
# >=4.0 CASA
# use opacity.ms file:
#gencal(vis='12A-339_sb9826039_1.56065.225594965275.ms', caltable='cal7_opacity.ms', caltype='opac', spw='0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49', parameter=tau)
# and gaincurve
#gencal(vis='12A-339_sb9826039_1.56065.225594965275.ms', caltable='cal7_gaincurve.ms', caltype='gceff')

def basicblock():
    # set flux scale for 3c286. note that we are not writing to the model_data column, which is only possible after casa 3.4.
    setjy(vis=msname, field='2', modimage='3C286_'+bandname+'.im', standard='Perley-Butler 2010', scalebychan=True, usescratch=False)

    # delay cal
    gaincal(vis=msname, caltable=calname+'_'+bbn+'_K0.ms', field='2', spw='0~'+str(nspw-1)+':5~58', refant=refantn, gaintype='K', solint='inf', combine='scan', gaintable=[calname+'_antpos.ms', calname+'_gaincurve.ms', calname+'_opacity.ms'])

    # first phase cal before delay and bp cal
    gaincal(vis=msname, caltable=calname+'_'+bbn+'_G0.ms', field='2', refant=refantn, spw='0~'+str(nspw-1)+':27~36', gaintype='G', calmode='p', solint='int', minsnr=2, gaintable=[calname+'_antpos.ms', calname+'_gaincurve.ms', calname+'_opacity.ms',calname+'_'+bbn+'_K0.ms'])

    # bp cal
    bandpass(vis=msname, caltable=calname+'_'+bbn+'_B0.ms', field='2', spw='0~'+str(nspw-1)+':5~58', refant=refantn, bandtype='B', solint='inf', combine='scan', solnorm=True, gaintable=[calname+'_antpos.ms', calname+'_gaincurve.ms', calname+'_opacity.ms',calname+'_'+bbn+'_K0.ms',calname+'_'+bbn+'_G0.ms'], interp=['','','','nearest','nearest'])
    
    # second phase cal
    gaincal(vis=msname, caltable=calname+'_'+bbn+'_G1scanph.ms', field='2,5,9', refant=refantn, spw='0~'+str(nspw-1)+':5~58', gaintype='G', calmode='p', solint='inf', minsnr=2, solnorm=False, gaintable=[calname+'_antpos.ms', calname+'_gaincurve.ms', calname+'_opacity.ms',calname+'_'+bbn+'_K0.ms',calname+'_'+bbn+'_B0.ms'], interp=['','','','nearest','nearest'])
    gaincal(vis=msname, caltable=calname+'_'+bbn+'_G1intph.ms', field='2,5,9', refant=refantn, spw='0~'+str(nspw-1)+':5~58', gaintype='G', calmode='p', solint='int', minsnr=2, solnorm=False, gaintable=[calname+'_antpos.ms', calname+'_gaincurve.ms', calname+'_opacity.ms',calname+'_'+bbn+'_K0.ms',calname+'_'+bbn+'_B0.ms'], interp=['','','','nearest','nearest'])
    gaincal(vis=msname, caltable=calname+'_'+bbn+'_G1scanamp.ms', field='2,5,9', refant=refantn, spw='0~'+str(nspw-1)+':5~58', gaintype='G', calmode='ap', solint='inf', minsnr=2, solnorm=False, gaintable=[calname+'_antpos.ms', calname+'_gaincurve.ms', calname+'_opacity.ms',calname+'_'+bbn+'_K0.ms',calname+'_'+bbn+'_B0.ms',calname+'_'+bbn+'_G1intph.ms'], interp=['','','','nearest','nearest','nearest'])

def polblock():
    # polcal. most params here come from either listobs or first use of setjy
    # first set the pol flux scale

    q0 = p0*cos(66*pi/180)
    u0 = p0*sin(66*pi/180)
    alpha0 = log(2.93237/2.87442)/log(18872/19384.)
        
    # field 5 (phase cal) parameterization of I, PA, and pol from first pass cal. freq (x) in MHz.
    line = lambda p, x: p[0] + p[1]*x
    # from pass 1
#    p_I = n.array([  8.07702338e-01,  -9.13761211e-06])
#    p_PA = n.array([  3.19703007e-01,   3.95346277e-05])
#    p_pol = n.array([  9.69603378e-03,   3.20534138e-08])
    # from pass 3
    p_I = n.array([  7.49212635e-01,  -7.59641505e-12 ])
    p_PA = n.array([ -8.19787881e-03,   4.92769912e-11 ])
    p_pol = n.array([  1.08576714e-02,  -2.32546199e-14 ])

    setjy(vis=msname, field='2', modimage='', standard='Perley-Butler 2010', scalebychan=True, fluxdensity=[i0,q0,u0,0], spix=alpha0, reffreq=reffreqname, usescratch=False)

#    for spwn in range(nspw):
#        freq = int(reffreqname[:-3]) + 128 * spwn
#        i1 = line(p_I, freq)
#        q1 = line(p_pol, freq) * cos(line(p_PA, freq))
#        u1 = line(p_pol, freq) * sin(line(p_PA, freq))
#        setjy(vis=msname, field='5', modimage='', standard='Perley-Butler 2010', spw=str(spwn), scalebychan=False, usescratch=False, fluxdensity=[i1,q1,u1,0], reffreq=reffreqname)

    # next calibrate cross hand delay
    gaincal(vis=msname, caltable=calname+'_'+bbn+'_Kcross.ms', field='2', spw='0~'+str(nspw-1)+':5~58', gaintype='KCROSS', solint='inf', combine='scan', refant=refantn, gaintable=[calname+'_antpos.ms', calname+'_gaincurve.ms', calname+'_opacity.ms',calname+'_'+bbn+'_K0.ms', calname+'_'+bbn+'_B0.ms', calname+'_'+bbn+'_G1intph.ms', calname+'_'+bbn+'_G1scanamp.ms'], gainfield=['','','','','','2','2'], interp=['','','','nearest','nearest', 'nearest', 'nearest'], parang=True)
    # delays around 1 ns. seems ok

    # next solve for leakages
    # option 1 is to use field 5 (gain cal) with full pol solution, since 3c286 and 3c48 not transferring well to target fields.
    polcal(vis=msname, caltable=calname+'_'+bbn+'_D1.ms', field='5', spw='0~'+str(nspw-1)+':5~58', solint='inf', combine='scan', refant=refantn, poltype='D+QU', gaintable=[calname+'_antpos.ms', calname+'_gaincurve.ms', calname+'_opacity.ms',calname+'_'+bbn+'_K0.ms', calname+'_'+bbn+'_B0.ms', calname+'_'+bbn+'_G1intph.ms', calname+'_'+bbn+'_G1scanamp.ms', calname+'_'+bbn+'_Kcross.ms'], gainfield=['','','','','','5','5',''], interp=['','','','nearest','nearest','nearest', 'nearest', 'nearest'])
    #option 2 is to bootstrap by setting QU from model of field 5 (set with setjy)
#    polcal(vis=msname, caltable=calname+'_'+bbn+'_D1.ms', field='5', spw='0~'+str(nspw-1)+':5~58', solint='inf', combine='scan', refant=refantn, poltype='D', gaintable=[calname+'_antpos.ms', calname+'_gaincurve.ms', calname+'_opacity.ms',calname+'_'+bbn+'_K0.ms', calname+'_'+bbn+'_B0.ms', calname+'_'+bbn+'_G1intph.ms', calname+'_'+bbn+'_G1scanamp.ms', calname+'_'+bbn+'_Kcross.ms'], gainfield=['','','','','','5','5',''], interp=['','','','nearest','nearest','nearest','nearest','nearest'])
    
# solve for R-L pol angle
    polcal(vis=msname, caltable=calname+'_'+bbn+'_X1.ms', field='2', spw='0~'+str(nspw-1)+':5~58', solint='inf', combine='scan', refant=refantn, poltype='X', gaintable=[calname+'_antpos.ms', calname+'_gaincurve.ms', calname+'_opacity.ms',calname+'_'+bbn+'_K0.ms', calname+'_'+bbn+'_B0.ms', calname+'_'+bbn+'_G1intph.ms', calname+'_'+bbn+'_G1scanamp.ms', calname+'_'+bbn+'_Kcross.ms', calname+'_'+bbn+'_D1.ms'], gainfield=['','','','','','2','2','',''], interp=['','','','nearest','nearest','nearest', 'nearest', 'nearest', 'nearest'])

    # scale other calibrator flux scales
    myscale = fluxscale(vis=msname, caltable=calname+'_'+bbn+'_G1scanamp.ms', fluxtable=calname+'_'+bbn+'_fluxscale.ms', reference=['2'], transfer=['5,9'])

    # selfcal
    gaincal(vis=msname, caltable=calname+'_'+bbn+'_G2.ms', field='6,7', uvrange='>80klambda', gaintype='G', calmode='p', solint='inf', refant=refantn, spw='0~'+str(nspw-1)+':5~58', solnorm=False, gaintable=[calname+'_antpos.ms', calname+'_gaincurve.ms', calname+'_opacity.ms',calname+'_'+bbn+'_G1scanph.ms',calname+'_'+bbn+'_fluxscale.ms',calname+'_'+bbn+'_K0.ms',calname+'_'+bbn+'_B0.ms',calname+'_'+bbn+'_Kcross.ms',calname+'_'+bbn+'_D1.ms',calname+'_'+bbn+'_X1.ms'], gainfield=['','','','5','5','','','','',''], interp=['','','','linear','linear','nearest','nearest','nearest','nearest','nearest'])


def applyblock():
    # apply to each field. first 'nearest' for each cal field, then 'linear' interp for target fields
    applycal(vis=msname, field='2', gaintable=[calname+'_antpos.ms', calname+'_gaincurve.ms', calname+'_opacity.ms',calname+'_'+bbn+'_G1intph.ms',calname+'_'+bbn+'_fluxscale.ms',calname+'_'+bbn+'_K0.ms',calname+'_'+bbn+'_B0.ms',calname+'_'+bbn+'_Kcross.ms',calname+'_'+bbn+'_D1.ms',calname+'_'+bbn+'_X1.ms'], gainfield=['','','','2','2','','','','',''], interp=['','','','nearest','nearest','nearest','nearest','nearest','nearest','nearest'], parang=True, calwt=False)
    applycal(vis=msname, field='5', gaintable=[calname+'_antpos.ms', calname+'_gaincurve.ms', calname+'_opacity.ms',calname+'_'+bbn+'_G1intph.ms',calname+'_'+bbn+'_fluxscale.ms',calname+'_'+bbn+'_K0.ms',calname+'_'+bbn+'_B0.ms',calname+'_'+bbn+'_Kcross.ms',calname+'_'+bbn+'_D1.ms',calname+'_'+bbn+'_X1.ms'], gainfield=['','','','5','5','','','','',''], interp=['','','','nearest','nearest','nearest','nearest','nearest','nearest','nearest'], parang=True, calwt=False)
    applycal(vis=msname, field='9', gaintable=[calname+'_antpos.ms', calname+'_gaincurve.ms', calname+'_opacity.ms',calname+'_'+bbn+'_G1intph.ms',calname+'_'+bbn+'_fluxscale.ms',calname+'_'+bbn+'_K0.ms',calname+'_'+bbn+'_B0.ms',calname+'_'+bbn+'_Kcross.ms',calname+'_'+bbn+'_D1.ms',calname+'_'+bbn+'_X1.ms'], gainfield=['','','','9','9','','','','',''], interp=['','','','nearest','nearest','nearest','nearest','nearest','nearest','nearest'], parang=True, calwt=False)
    applycal(vis=msname, field='6,7', gaintable=[calname+'_antpos.ms', calname+'_gaincurve.ms', calname+'_opacity.ms',calname+'_'+bbn+'_G1scanph.ms',calname+'_'+bbn+'_fluxscale.ms',calname+'_'+bbn+'_G2.ms',calname+'_'+bbn+'_K0.ms',calname+'_'+bbn+'_B0.ms',calname+'_'+bbn+'_Kcross.ms',calname+'_'+bbn+'_D1.ms',calname+'_'+bbn+'_X1.ms'], gainfield=['','','','5','5','6,7','','','','',''], interp=['','','','linear','linear','nearest','nearest','nearest','nearest','nearest','nearest'], parang=True, calwt=False)

# run blocks sequentially
basicblock()
polblock()
applyblock()
