import numpy as n
import leanpipe_cython as lib
try:
    import pyfftw
except:
    pass

try:
    import fftw3
except:
    pass

def imgone(u, v, data, size, res, method='', kernel=''):
    """ Hack of aipy imager
    u,v,data should be 1d
    size is extent of uv grid
    res is uv resolution of grid
    (so size/res is number of uv bins)
    kernel can be '' (only binning) or 'half'
    """

    # initial definitions
    ndim = n.round(size/res)
    grid = n.zeros((ndim,ndim), dtype='complex128')

    # append hermition to get real image. not really needed?
#    u = n.concatenate([u, -u], axis=0)
#    v = n.concatenate([v, -v], axis=0)
#    data = n.concatenate([data, n.conj(data)], axis=0)

    # put uv data on grid
# testing simpler(?) ways
    u = n.mod(n.round(u/res).astype(n.int), ndim)
    v = n.mod(n.round(v/res).astype(n.int), ndim)
    inds = n.vstack((v,u))
# old way
#    u = n.round(u/res).astype(n.int)
#    v = n.round(v/res).astype(n.int)#    inds = n.concatenate([-v,u],).transpose()
#    ok = n.logical_and(n.abs(inds[:,0]) < ndim, n.abs(inds[:,1]) < ndim)
#    data = data.compress(ok)
#    inds = inds.compress(ok, axis=0)

    # add uv data to grid
    if kernel == '':
        for i in xrange(len(u)):
            cellu, cellv = inds[:,i]
            if ( (n.abs(cellu) < ndim) & (n.abs(cellv) < ndim) ):
                grid[cellu,cellv] += data[i]
    elif kernel == 'half':
        for i in xrange(len(u)):
            cellu, cellv = inds[:,i]
            if ( (n.abs(cellu) > 1) & (n.abs(cellv) > 1) & (n.abs(cellu) < ndim-1) & (n.abs(cellv) < ndim-1) ):
                data0 = grid[cellu-1:cellu+2,cellv-1:cellv+2]
                grid[cellu-1:cellu+2,cellv-1:cellv+2] = data0 + convvis(data[i])

    if method=='pyfftw':
        im = pyfftw.interfaces.numpy_fft.ifft2(grid).real.astype(n.float32)
    elif method=='fftw3':
        im = n.zeros((ndim,ndim), dtype='complex128')
        ifft = fftw3.Plan(grid, im, direction='backward')
        ifft.execute()
    else:
        im = n.fft.ifft2(grid).real.astype(n.float32)

    print 'Pixel size %.1f\", Field size %.1f\"' % (3600*n.degrees(1./size), 3600*n.degrees(1./res))

    return n.roll(n.roll(im, ndim/2, axis=0), ndim/2, axis=1)   # recenter then return image

def imgall(u, v, data, size, res):
    """ Hack of aipy imager
    u,v should be 1d. data should have extra time dimension
    size is extent of uv grid
    res is uv resolution of grid
    (so size/res is number of uv bins)
    kernel can be '' (only binning) or 'half'
    """

    # initial definitions
    ndim = n.round(size/res)
    grid = n.zeros((len(data),ndim,ndim), dtype='complex128')

    # put uv data on grid
    u = n.mod(n.round(u/res).astype(n.int), ndim)
    v = n.mod(n.round(v/res).astype(n.int), ndim)
    inds = n.vstack((v,u))

    # add uv data to grid
    ims = []
    for i in xrange(len(u)):
        cellu, cellv = inds[:,i]
#        if ( (n.abs(cellu) < ndim) & (n.abs(cellv) < ndim) ):
        grid[:,cellu,cellv] += data[:,i]

    for t in xrange(len(data)):
        im = n.fft.ifft2(grid[t]).real.astype(n.float32)
        ims.append(n.roll(n.roll(im, ndim/2, axis=0), ndim/2, axis=1))
    print 'Pixel size %.1f\", Field size %.1f\"' % (3600*n.degrees(1./size), 3600*n.degrees(1./res))

    return ims

def convvis(vis):
    """ Defines normalized convolution kernel that is multiplied by visibility.
    """
    conv = n.array([ [0., 0.125, 0.], [0.125, 0.5, 0.125], [0., 0.125, 0.] ])

    return vis*conv

def gendata(nint, nbl, nch, npol, bytenum, pyfftw=False):
    """ Generate data of shape (nbl, nch, npol)
    """

#    if pyfftw:
#        data = pyfftw.n_byte_align_empty((nint, nbl, nch, npol), bytenum, 'complex128')
#    else:
    data = n.zeros((nint, nbl, nch, npol), dtype='complex128')
#    data[:] = n.random.normal(-1,1,(nint,nbl,nch,npol)) + 1j*n.random.normal(-1,1,(nint,nbl,nch,npol))
    data = n.random.randn(nint,nbl,nch,npol) + 1j*n.random.randn(nint,nbl,nch,npol)
    return data

def genuv(nint, nbl, nch, maxuv, freqs):
    """ Generate u, v
    """

    u = n.zeros((nint, nbl)); v = n.zeros((nint, nbl))
    u[:] = n.random.uniform(-maxuv, maxuv, nbl)
    v[:] = n.random.uniform(-maxuv, maxuv, nbl)
    
#    u = n.outer(u, freqs/freqs[0])
#    v = n.outer(v, freqs/freqs[0])

    return u, v

def addsource(data, i=1.):
    data = data + (i+0j)
    return data

def phaseshift(data, u, v, dl, dm):
    """ Shift phase center
    """

    ang = lambda dl,dm,u,v: dl*u + dm*v
    for i in xrange(data.shape[0]):
        for j in xrange(data.shape[1]):
            data[i,j] = data[i,j] * n.exp(-2j*n.pi*ang(dl, dm, u[i,j], v[i,j]))
    return data

def genall(nint, nbl, nch, npol, maxuv, bwfrac):
    """ freq0 is first freq in MHz, chsize is channel size in MHz
    """

    d={}
    d['freq'] = n.arange(1200,1200+128)/1000.; d['freq_orig'] = n.arange(1200,1200+128)/1000.; d['l0'] = 0.; d['m0'] = 0.
    u, v = genuv(nint, nbl, nch, maxuv, n.arange(1000,1000*(1+bwfrac), 1000*bwfrac/nch)/1000.)
    data = gendata(nint, nbl, nch, npol, 16)
    data = addsource(data)
    lib.phaseshift(data, d, n.radians(0.125), 0., u, v)
    data = addsource(data)
    lib.phaseshift(data, d, n.radians(0.25), 0., u, v)
    data = addsource(data)

#    t1 = addsource(t1)
#    t2 = lib.phaseshift(data, d, n.radians(0.2), 0., u[:,0], v[:,0])
#    t2 = addsource(t2)
#    t3 = lib.phaseshift(data, d, n.radians(0.3), 0., u[:,0], v[:,0])
#    t3 = addsource(t3)
#    t4 = lib.phaseshift(data, d, n.radians(0.4), 0., u[:,0], v[:,0])
#    t4 = addsource(t4)
#    data = n.concatenate( ([data], [data], [data], [t1], [data], [data], [t2], [data], [data], [t3], [data], [data], [t4], [data], [data], [data]))

    return data, u, v
