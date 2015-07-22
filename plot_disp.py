#! /usr/bin/env python
# claw, 4aug11
#
# Script to visualize effects of dispersion and scattering on pulse detection

import numpy as n
import pylab as p

# define functional forms for dispersion drift and scattering
def drift(dm,f,bw=0.128):   # f,bw in GHz
    return 4.15e-3 * dm * (1/f**2 - 1/(f+bw)**2)

def scat1(dm,f): 
    return 0.01 * (f/3.)**(-4) * (dm/1000.)**2     # needs to be refined

def scat(dm,f):   # bhat et al. 2004. gives t_d in s
    return 10**(-9.46 + 0.154*n.log10(dm) + 1.07*(n.log10(dm))**2 - 3.86*n.log10(f))

bw = 0.1   # bandwidth in GHz
tint = 1.   # integration time in seconds
dmarr = n.arange(10,1000,10)   # dm grid in pc/cm3
farr = n.arange(10,30)/100.    # frequency grid in GHz

# initialize arrays
driftarr = [drift(dmarr, farr[0], bw=bw)]
scatarr = [scat(dmarr, farr[0])]

for f in farr[1:]:
    driftarr = n.append(driftarr, [drift(dmarr, f, bw=bw)], axis=0)
    scatarr = n.append(scatarr, [scat(dmarr, f)], axis=0)

# ratio of drift to scatter
raarr = driftarr/scatarr
zz = n.zeros(n.shape(raarr))
mask = n.where(driftarr > tint, raarr, zz)   # area where drift is larger than integration time

# plot
p.imshow(n.log10(mask), aspect='auto', interpolation='nearest', origin='lower', extent= (dmarr.min(), dmarr.max(), farr.min(), farr.max()))

p.title('Log of Drift time over scattering time, %d integrations' % (tint))
p.xlabel('DM (pc cm$^{-3}$)')
p.ylabel('Frequency (GHz)')
p.colorbar()
p.show()

# can add plot of SNR of pulse peak, assuming spectral index
# can add plot of drift SNR
