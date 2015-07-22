import time

def flux(field, scan, spw, bb, chanst=':5~58', prefix='test'):
    Qm = imstat(prefix+field+'_scan'+scan+'_spw'+spw+chanst+'_bb'+bb+'_IQUV.image', stokes='Q', box='500,500,500,500')['max'][0]
    Qe = imstat(prefix+field+'_scan'+scan+'_spw'+spw+chanst+'_bb'+bb+'_IQUV.residual', stokes='Q', box='600,600,900,900')['rms'][0]
    Um = imstat(prefix+field+'_scan'+scan+'_spw'+spw+chanst+'_bb'+bb+'_IQUV.image', stokes='U', box='500,500,500,500')['max'][0]
    Ue = imstat(prefix+field+'_scan'+scan+'_spw'+spw+chanst+'_bb'+bb+'_IQUV.residual', stokes='U', box='600,600,900,900')['rms'][0]
    I = imstat(prefix+field+'_scan'+scan+'_spw'+spw+chanst+'_bb'+bb+'_IQUV.image', stokes='I', box='400,400,600,600')['max'][0]
    pol = n.sqrt(Qm**2 + Um**2)
    polfrac = pol/I
    theta = n.angle(n.complex(Qm, Um))

    print 'field, bb, spw, Qm, Qe, Um, Ue, I'
    print str(field)+','+str(bb)+','+str(spw)+','+str(Qm)+','+str(Qe)+','+str(Um)+','+str(Ue)+','+str(I)
    print 'pol, polfrac, theta'
    print pol, polfrac, theta
    return [float(field), float(bb), float(spw), Qm, Qe, Um, Ue, I]

def plot(field='6', scan='117~512', chanst=':5~58', prefix='test'):
    freqdict = {'48': [4.79880000e+10,   4.81160000e+10,   4.82440000e+10,  4.83720000e+10,   4.85000000e+10,   4.86280000e+10,  4.87560000e+10,   4.88840000e+10],
                '39': [3.92940000e+10,   3.94220000e+10,  3.95500000e+10,   3.96780000e+10,   3.98060000e+10,  3.99340000e+10],
                '37': [3.74780000e+10,  3.76060000e+10,   3.77340000e+10,   3.78620000e+10, 3.79900000e+10,   3.81180000e+10, 3.82460000e+10,  3.83740000e+10],
                '27': [2.69680000e+10,   2.70960000e+10,   2.72240000e+10,   2.73520000e+10,   2.74800000e+10,  2.76080000e+10, 2.77360000e+10,   2.78640000e+10],
                '25': [2.54880000e+10,   2.56160000e+10,  2.57440000e+10,   2.58720000e+10,   2.60000000e+10,  2.61280000e+10, 2.62560000e+10,   2.63840000e+10],
                '18': [1.88720000e+10,   1.90000000e+10,   1.91280000e+10,  1.92560000e+10,   1.93840000e+10]
                }

    fluxes = []
    for bb in ['18','25','27','37','39','48']:
        for spw in ['0','1','2','3','4','5','6','7']:
            try:
                fluxvals = flux(field, scan, spw, bb, chanst=chanst, prefix=prefix)
                fluxes.append(fluxvals)
            except:
                print 'Skipping (spw, bb):', spw, bb
    fluxes = n.array(fluxes)

    freqs = []; q = []; qe = []; u = []; ue = []; angle = []; pol = []; I = []
    for ff in fluxes:
        freqs.append(freqdict[str(int(ff[1]))][int(ff[2])])
        q.append(float(ff[3]))
        qe.append(float(ff[4]))
        u.append(float(ff[5]))
        ue.append(float(ff[6]))
        pol.append(n.sqrt(float(ff[3])**2 + float(ff[5])**2))
        I.append(float(ff[7]))
        angle.append(n.angle(n.complex(float(ff[3]),float(ff[5]))))

    freqs = n.array(freqs)
    q = n.array(q); qe = n.array(qe); u = n.array(u); ue = n.array(ue); I = n.array(I); pol = n.array(pol)
    perr = n.sqrt((q*qe)**2+(u*ue)**2)/n.sqrt(q**2+u**2)
    polfrac = pol/I

    print 'resting...'
    time.sleep(3)

    p.figure(1)
    p.subplot(211)
    p.errorbar(freqs/1e9, q, yerr=qe, fmt='.', label='field '+str(field))
    p.ylabel('Q (Jy)')
    p.subplot(212)
    p.errorbar(freqs/1e9, u, yerr=ue, fmt='.', label='field '+str(field))
    p.ylabel('U (Jy)')
    p.xlabel('Freq (GHz)')
    p.legend()

    p.figure(2)
    p.subplot(211)
    p.plot(freqs/1e9, angle, '.', label='field '+str(field))
    p.ylabel('Pol. angle (rad)')
    p.subplot(212)
    p.errorbar(freqs/1e9, polfrac, yerr=perr/I, fmt='.', label='field '+str(field))
    p.ylabel('Pol frac')
    p.xlabel('Freq (GHz)')
    p.legend()
    
    p.figure(3)
    p.errorbar(q, u, xerr=qe, yerr=ue, label='field '+str(field))
    p.xlabel('Q (Jy)')
    p.ylabel('U (Jy)')
    p.legend()

#    return n.array(q), n.array(qe), n.array(u), n.array(ue)
    return freqs, angle, polfrac, I
