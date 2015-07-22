#!/usr/bin/env python
# claw, 9jul09
# 
# script to estimate fast transient rates

import pylab as p
import numpy as n

def parkes(dt=10.):
    """
    Estimate of rate of transients detection based on Parkes MB pulsar survey.
    Use results for basic PMPS data and RRATs published in McLaughlin+ 2006, Keane+ 2011. 
    Normalize by area coverage and sensitivity for given dt, then scaled to VLA L-band.
    """

    # pmps RRATs
    # **may need to add from keane+ 2010.**
    # **also, this excludes the roughly 30% of normal pulsars seen in single-pulse searches.**
    # **also, no zero-dm events (flare stars, dwarfs, etc.) are expected to be seen.**
    rratrates = [27/19., 108/24., 32/41., 18/30., 229/13., 18/17., 8/13., 11/10., 10/8., 4/13., 66/14.]  # mclaughlin+ 2006
    rratrates = rratrates + [3.8, 3.2, 18.7, 52.0, 0.7, 1.7, 8.9, 3.2, 0.6, 8.2, 63.4, 8.9, 1.8, 0.4, 0.6, 0.3, 1.1, 0.3, 0.2, 18.9, 19.2, 57.9] # keane+ 2011
    rratw50 = [30.,10.,20.,16.,3.,2.,15.,16.,2.,5.,2.] # mclaughlin+ 2006
    rratw50 = rratw50 + [4.,5.,3.,1.,64.,9.,12.,6.,7.,4.,3.,16.,2.,16.,2.,16.,4.,16.,7.,2.,12.,1.] # keane+ 2011

    # rate of events for flux greater than fluxarr. plaw gives (assumed) powerlaw of pulse energy distribution.
    slimarr = n.arange(35,110.,2)  # s_lim in mJy
    plaw = -2 # assumed powerlaw slope of cumulative pulse distribution. based on powerlaw of -1 for pulsar luminosity distribution.

    # two functions to calculate base sensitivity of parkes and cumulative rate for new s
#    s_parkes = lambda dt: 100*(10./dt)**(1/2.) # 5sigma is 100 mJy in 10 ms
#    c_rate = lambda fluxarr, dt, s_parkes, plaw, pmpsrate, w50: pmpsrate * (fluxarr*(dt/w50)/s_parkes(w50))**plaw
    # two other functions to do scale from limits at fixed dt
    s_parkes = 100 # 5sigma is 100 mJy in 10 ms
    c_rate = lambda fluxarr, dt, s_parkes, plaw, pmpsrate, w50: pmpsrate * (fluxarr*(dt/w50)/s_parkes)**plaw

    # build total rate
    p.figure(1)
    p.clf()
    total_dt = c_rate(slimarr, max(dt,rratw50[0]), s_parkes, plaw, rratrates[0], rratw50[0]) # line for inferring rate at given dt
    total_orig = c_rate(slimarr, max(0.1,rratw50[0]), s_parkes, plaw, rratrates[0], rratw50[0]) # line for inferring original rate for Parkes
    for i in range(1,len(rratrates)):
        total_dt += c_rate(slimarr, max(dt,rratw50[i]), s_parkes, plaw, rratrates[i], rratw50[i])
        total_orig += c_rate(slimarr, max(0.1,rratw50[i]), s_parkes, plaw, rratrates[i], rratw50[i])
    p.plot(slimarr, total_dt, label='Rate @ %s ms' % (str(dt)))   # rate for apparent flux and dt
    p.plot(slimarr, total_orig, label='Rate for Parkes dt (0.1 ms or something)')   # rate for apparent flux and dt
    p.xlabel('Flux limit at 10 ms (mJy)')
    p.ylabel('Total Parkes pulse rate (hr$^{-1}$)')
    p.text(100, total_orig[(n.abs(slimarr-100)).argmin()], 'Parkes', rotation=90, horizontalalignment='center', verticalalignment='center')
    p.text(38, total_dt[(n.abs(slimarr-38)).argmin()], 'VLA L bisp (1 GHz)', rotation=90, horizontalalignment='center', verticalalignment='center')
    p.legend()

    # add second y axis scaled by total PMPS area to single VLA l-band field of view (0.2/1500 sq deg)
    p.twinx()
    p.plot(slimarr, total_orig*(0.2/1500), ls='None') # make it white to blank it out
    p.ylabel('VLA L-band pulse rate (hr$^{-1}$ f.o.v.$^{-1}$)')
    p.title('Rates with %s ms integrations' % (str(dt)))



def crabgp():
    """
    Plot the rate of crab giant pulses visible to ata-256 as a function of distance
    See McLaughlin & Cordes 2003
    """

    flux = n.arange(50.)/10. + 2
    alpha = -3.4

#    mode = 'f'  # plot flux vs. frequency
    mode = 'd'  # plot flux vs. distance

    if mode == 'd':
        # the distance estimate also assumes dedispersion and/or matched filter detection (i.e., idealized)
        rate = (24/70.)*10**11.1*(10**flux)**(-2.5)  # cumulative number with S > "flux" per day for rate in McLaughlin & Cordes 2003
        dist812 = 0.85 * (26/5.)**(-1/2.) * (286/10.)**(1/4.) * ((10.**flux)/(10.**5))**(1/2.)  # BW=60 MHz, Ssys=26Jy, dt=0.1ms, SNR=5, 2pols, bottom at 800 MHz, weighted mean at 812 MHz (as in Lundgren)
        dist632 = 0.85 * (26/5.)**(-1/2.) * (500/10.)**(1/4.) * ((((632/812.)**alpha)*10.**flux)/(10.**5))**(1/2.)  # BW=30 MHz, Ssys=26Jy, dt=0.1ms, SNR=5, 2pols, bottom at 500 MHz, weighted mean at 514 MHz
        dist527 = 0.85 * (26/5.)**(-1/2.) * (60/10.)**(1/4.) * ((((527.5/812.)**alpha)*10.**flux)/(10.**5))**(1/2.)  # BW=60 MHz, Ssys=26Jy, dt=0.1ms, SNR=5, 2pols, bottom at 500 MHz, weighted mean at 527 MHz
        dist588 = 0.85 * (26/5.)**(-1/2.) * (60/10.)**(1/4.) * ((((588/812.)**alpha)*10.**flux)/(10.**5))**(1/2.)  # BW=60 MHz, Ssys=26Jy, dt=0.1ms, SNR=5, 2pols, bottom at 560 MHz, weighted mean at 588 MHz
        dist1614 = 0.85 * (26/5.)**(-1/2.) * (500/10.)**(1/4.) * ((((1614/812.)**alpha)*10.**flux)/(10.**5))**(1/2.)  # BW=500 MHz, Ssys=26Jy, dt=0.1ms, SNR=5, 2pols, bottom at 1420 MHz, weighted mean at 1614 MHz
#dist501_qb = 0.85 * (26/5.)**(-1/2.) * (15/10.)**(1/4.) * ((((501/812.)**alpha)*10.**flux)/(10.**5))**(1/2.)  # BW=60 MHz, Ssys=26Jy, dt=0.1ms, SNR=5, 2pols, bottom at 500 MHz, weighted mean at 502 MHz

    # alt numbers for 10 ms time resolution.  assumes GPs are underresolved by factor of three (1.5 ms integration/10 ms burst), but gains in noise and bandwidth
        flux = n.arange(50,100000,5.)
        ratealt = (24/70.)*10**11.1*(flux*3*(10/1.5))**(-2.5)  # cumulative number with S > "flux" per day for rate in Lundgren et al. 1995 averaged from 1.5 ms down to 10 ms
        dist812alt = 0.85 * (156/5.)**(-1/2.) * (600/10. * 10/0.1)**(1/4.) * (flux/(10.**5))**(1/2.)  # BW=60 MHz, Ssys=26Jy, dt=0.1ms, SNR=5, 2pols, bottom at 800 MHz, weighted mean at 812 MHz (as in Lundgren)
        dist674alt = 0.85 * (156/5.)**(-1/2.) * (600/10. * 10/0.1)**(1/4.) * ((((674/812.)**alpha)*flux)/(10.**5))**(1/2.)  # BW=30 MHz, Ssys=26Jy, dt=0.1ms, SNR=5, 2pols, bottom at 500 MHz, weighted mean at 514 MHz
        dist1624alt = 0.85 * (156/5.)**(-1/2.) * (600/10. * 10/0.1)**(1/4.) * ((((1624/812.)**alpha)*flux)/(10.**5))**(1/2.)  # BW=500 MHz, Ssys=26Jy, dt=0.1ms, SNR=5, 2pols, bottom at 1420 MHz, weighted mean at 1624 MHz

#l1 = p.loglog(dist812alt,ratealt, label='1 Crab @ 812 MHz')  # plot for 1 Crab at 812 MHz
        l1 = p.loglog(dist674alt,ratealt, label='1 Crab @ 674 MHz')
        l1 = p.loglog(dist674alt,10*ratealt, label='10 Crab @ 674 MHz')
        l1 = p.loglog(dist1624alt,10*ratealt, label='10 Crab @ 1624 MHz')
#l1 = p.loglog(dist812,rate, label='1 Crab @ 812 MHz')  # plot for 1 Crab at 812 MHz
#l2 = p.loglog(dist812,10*rate, label='10 Crabs @ 812 MHz')
#l3 = p.loglog(dist527,10*rate, label='10 Crabs @ 527 MHz')
#l4 = p.loglog(dist632,10*rate, label='10 Crabs @ 632 MHz')
#l4 = p.loglog(dist588,10*rate, label='10 Crabs @ 588 MHz')
#l5 = p.loglog(dist1614,10*rate, label='10 Crabs @ 1614 MHz')
#l6 = p.loglog(dist527,1000*rate, label='1000 Crabs @ 527 MHz')
#l6 = p.loglog(dist501_qb,10*rate, label='10 Crabs @ 501 MHz (quarter band)')

        top = 6e3
#p.loglog([0.78,0.78],[0.1*top,top],'b--',label='M31')
        p.text(0.78,top,'M31 ',rotation='vertical',verticalalignment='top',horizontalalignment='center')
#p.loglog([0.86,0.86],[0.1*top,top],'b--',label='M33')
        p.text(0.86,top,'M33 ',rotation='vertical',verticalalignment='top',horizontalalignment='center')
#p.loglog([0.1,0.1],[0.1*top,top],'b--',label='UMa I')
        p.text(0.1,top,'UMa I ',rotation='vertical',verticalalignment='top',horizontalalignment='center')
#p.loglog([0.06,0.06],[0.1*top,top],'b--',label='UMi, Bootes')
        p.text(0.06,top,'UMi, Bootes ',rotation='vertical',verticalalignment='top',horizontalalignment='center')
#p.loglog([0.09,0.09],[0.1*top,top],'b--',label='Sextans I')
        p.text(0.09,top,'Sextans I ',rotation='vertical',verticalalignment='top',horizontalalignment='center')
#p.loglog([0.08,0.08],[0.1*top,top],'b--',label='Draco')
        p.text(0.08,top,'Draco ',rotation='vertical',verticalalignment='top',horizontalalignment='center')
#p.loglog([0.25,0.25],[0.1*top,top],'b--',label='Leo I')
        p.text(0.25,top,'Leo I ',rotation='vertical',verticalalignment='top',horizontalalignment='center')
#p.loglog([0.21,0.21],[0.1*top,top],'b--',label='Leo II')
        p.text(0.21,top,'Leo II ',rotation='vertical',verticalalignment='top',horizontalalignment='center')
#p.loglog([0.5,0.5],[0.1*top,top],'b--',label='NGC 6822')
        p.text(0.5,top,'NGC 6822 ',rotation='vertical',verticalalignment='top',horizontalalignment='center')
#p.loglog([18,18],[0.1*top,top],'b--',label='Virgo Cl.')
#p.text(18,0.1*top,'Virgo Cl.',rotation='vertical')

# pretty up

        p.legend(('1 Crab @ 674 MHz','10 Crab @ 674 MHz'),loc=3)
#p.legend(('10 Crabs @ 527 MHz','10 Crabs @ 588 MHz','10 Crabs @ 1614 MHz','1000 Crabs @ 527 MHz'),loc=4)
#p.legend(('1 Crab @ 812 MHz','10 Crabs @ 812 MHz','10 Crabs @ 527 MHz','10 Crabs @ 588 MHz'),loc=3)
#p.legend(('10 Crabs @ 527 MHz','10 Crabs @ 632 MHz'),loc=3)
        p.axis([0.055,0.91,1e-2,top])
        p.title('Rate of Detectable Crab-like GPs (for BW=600 MHz, dt=10 ms)')
#p.title('Rate of Detectable Crab-like GPs (for DM=50 pc cm-3, dt=0.1 ms)')
        p.ylabel('Rate (>D; per day)')
        p.xlabel('Distance visible to ATA-42 (Mpc)')
        p.show()

    elif mode == 'f':
        freqs = n.arange(5,100)/10.
        rate = lambda flux,alpha,freq,fsc,tmin,diam: (206265*3e-1/freq/diam/60/5.0)**2 * 5e10 * ((freq/0.812)**(alpha) * (8/2.)**2 * 10*n.sqrt(tmin**2 + (fsc/freq)**(4.4*2)) * flux)**(-2.5)

    # note big dependence on relative values of spectral index and pulse scattering index.
        fsc = 4.0

        norm100 = 1/(rate(0.0065, 3.4, 2.0, 4.0, 100., 100.) * 29/24.)
        norm10 = 1/(rate(0.0065*n.sqrt(100/10.), 3.4, 2.0, 4.0, 10., 100.) * 29/24.)
    # limits and rates of number of pulses from a single gc crab pulsar
#    p.plot(freqs, 10*rate(0.0065, freqs, 17., 10 ))  # macquart et al. 2010 pulsars in central 50" for 10 hours (single pulse search)
        p.plot(freqs, norm10*rate(0.035, 3.4, freqs, fsc, 10., 25.), 'b', label='10 ms @ EVLA')  # evla proposal including survey speed 6x better than d09
        p.plot(freqs, norm100*rate(0.035/n.sqrt(10), 3.4, freqs, fsc, 100., 25.), 'b--', label='100 ms @ EVLA')  # evla proposal including survey speed 6x better than d09
        p.plot(freqs, norm10*rate(0.0065*n.sqrt(100/10.), 3.4, freqs, fsc, 10., 100.), 'g', label='10 ms @ D09') # deneva pulsars in 5.8' beam at 2 ghz for 29x1 hours
        p.plot(freqs, norm100*rate(0.0065, 3.4, freqs, fsc, 100., 100.), 'g--', label='100 ms @ D09') # deneva pulsars in 5.8' beam at 2 ghz for 29x1 hours
        p.plot(2., 24/29., 'g*')
#    p.errorbar(2., 24/29., yerr=0.3, capsize=4, uplims=0, lolims=1, fmt='g')
        p.loglog()
        p.xlabel('Frequency (GHz)')
        p.ylabel('GC Crab pulse rate (day$^{-1}$)')
        p.legend()
        p.show()
