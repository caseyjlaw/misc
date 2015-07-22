import numpy as n

def imgone(u, v, w, data, size, res):

    # initial definitions
    ndim = n.round(size/res)
    uv = n.zeros((ndim,ndim), dtype='complex128')

    # append hermition to get real image
    u = n.concatenate([u, -u], axis=0)
    v = n.concatenate([v, -v], axis=0)
    w = n.concatenate([w, -w], axis=0)
    data = n.concatenate([data, n.conj(data)], axis=0)

    # put uv data on grid
    u = n.round(u/res).astype(n.int)
    v = n.round(v/res).astype(n.int)
    inds = n.array([-v,u],).transpose()
    ok = n.logical_and(n.abs(inds[:,0]) < ndim, n.abs(inds[:,1]) < ndim)
    data = data.compress(ok)
    inds = inds.compress(ok, axis=0)

    # add uv data to grid
    for i in range(len(inds)):
        cellu, cellv = inds[i]
        uv[cellu,cellv] += data[i]

    im = n.fft.ifft2(uv).real.astype(n.float32)
    return n.roll(n.roll(im, ndim/2, axis=0), ndim/2, axis=1)   # recenter then return image

def imgall(u0, v0, w0, data0, size, res):

    # initial definitions
    ndim = n.round(size/res)
    sh = data0.shape
    uv = n.zeros((sh[0], ndim, ndim), dtype='complex128')

    # put uv data on grid
    u0 = n.round(u0/res).astype(n.int)
    v0 = n.round(v0/res).astype(n.int)
    inds = n.array([-v0,u0],).transpose()
#    ok = n.logical_and(n.abs(inds[:,0]) < ndim, n.abs(inds[:,1]) < ndim)
#    data = data.compress(ok)
#    inds = inds.compress(ok, axis=0)

    # add uv data to grid
    for j in xrange(len(inds)):
        cellu, cellv = inds[j]
        if ((cellu > -ndim) & (cellu < ndim) & (cellv > -ndim) & (cellv < ndim)):
            uv[:,cellu,cellv] += data0[:,j].mean(axis=1)  # mean over pols

    ims = []
    for i in xrange(len(uv)):
        im = n.fft.ifft2(uv[i]).real.astype(n.float32)
        im = n.roll(n.roll(im, ndim/2, axis=0), ndim/2, axis=1)   # recenter then return image
        ims.append(im)
    return ims
