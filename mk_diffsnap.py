#! /usr/bin/env python
"""
This is a general-purpose script for making images from MIRIAD UV files.  Data
(optionally selected for baseline, channel) are read from the file, phased
to a provided position, normalized for passband/primary beam effects, gridded
to a UV matrix, and imaged
"""

import aipy as a, numpy as n, sys, optparse, ephem, os

o = optparse.OptionParser()
o.set_usage('mk_img.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, ant=True, pol=True, chan=True, cal=True,
    src=True, dec=True)
o.add_option('-o', '--output', dest='output', default='dim,dbm',
    help='Comma delimited list of data to generate FITS files for.  Can be: dim (dirty image), dbm (dirty beam), uvs (uv sampling), or bms (beam sampling).  Default is dim,dbm.')
o.add_option('--list_facets', dest='list_facets', action='store_true',
    help='List the coordinates of all the pointings that will be used.')
o.add_option('--facets', dest='facets', type='int', default=200,
    help='If no src is provided, facet the sphere into this many pointings for making a map.  Default 200.')
o.add_option('--fmt', dest='fmt', default='im%04d',
    help='A format string for counting successive images written to files.  Default is im%04d (i.e. im0001).')
o.add_option('--skip_phs', dest='skip_phs', action='store_true',
    help='Do not phase visibilities before gridding.')
o.add_option('--skip_amp', dest='skip_amp', action='store_true',
    help='Do not use amplitude information to normalize visibilities.')
o.add_option('--skip_bm', dest='skip_bm', action='store_true',
    help='Do not weight visibilities by the strength of the primary beam.')
o.add_option('--bgsub', dest='bgsub', type='int', default=1,
    help='Define background subtraction window of +-1 integration. Default is no bg subtraction.')
o.add_option('--size', dest='size', type='int', default=500,
    help='Size of maximum UV baseline.')
o.add_option('--res', dest='res', type='float', default=0.5,
    help='Resolution of UV matrix.')
o.add_option('--no_w', dest='no_w', action='store_true',
    help="Don't use W projection.")
o.add_option('--wres', dest='wres', type='float', default=0.5,
    help="W-Plane projection resolution.  Default 0.5")
o.add_option('--altmin', dest='altmin', type='float', default=0,
    help="Minimum allowed altitude for pointing, in degrees.  When phase center is lower than this altitude, data is omitted.  Default is 0.")
o.add_option('--minuv', dest='minuv', type='float', default=0,
    help="Minimum distance from the origin in the UV plane (in wavelengths) for a baseline to be included.  Default is 0.")
o.add_option('--buf_thresh', dest='buf_thresh', default=2e6, type='float',
    help='Maximum amount of data to buffer before gridding.  Excessive gridding takes performance hit, but if buffer exceeds memory available... ouch.')
opts, args = o.parse_args(sys.argv[1:])

# Parse command-line options
uv = a.miriad.UV(args[0])
(j,t,j),j = uv.read()
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
a.scripting.uv_selector(uv, opts.ant, opts.pol)
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
aa.select_chans(chans)
aa.set_active_pol(opts.pol)
afreqs = aa[0].beam.afreqs
cfreq = n.average(afreqs)
aa.set_jultime(t)
del(uv)
outputs = opts.output.split(',')

# Get all sources that will be used as phase centers.  If no sources are
# specified, define phase centers for faceting a sphere.
if opts.src == 'zen':
    srcs = [a.phs.RadioFixedBody(aa.sidereal_time(), aa.lat, name='zen')]
    cat = a.phs.SrcCatalog(srcs)
elif not opts.src is None: 
    srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
    cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
else:
    ras,decs = a.map.facet_centers(opts.facets, ncrd=2)
    srcs = [a.phs.RadioFixedBody(ra,dec,name=str(i)) 
        for i,(ra,dec) in enumerate(zip(ras,decs))]
    cat = a.phs.SrcCatalog(srcs)

if opts.list_facets:
    cat.compute(aa)
    srcs = cat.keys(); srcs.sort()
    for cnt, src in enumerate(cat.values()):
        cen = ephem.Equatorial(src.ra, src.dec, epoch=aa.epoch)
        cen = ephem.Equatorial(cen, epoch=ephem.J2000)
        print '# %3d >  RA=%s  DEC=%s  (%f, %f in deg)' % \
            (cnt, cen.ra, cen.dec, 
            a.img.rad2deg*cen.ra, a.img.rad2deg*cen.dec)

# Generate the image object that will be used.
us,vs,ws,ds,wgts = [],[],[],[],[]
if opts.no_w:
    im = a.img.Img(opts.size, opts.res, mf_order=0)
else:
    im = a.img.ImgW(opts.size, opts.res, mf_order=0, wres=opts.wres)
L,M = im.get_LM()
DIM = int(opts.size/opts.res)
n_ints = 0

#print 'Calculating image of primary beam'
#top = im.get_eq(0, aa.lat)
#mask = top[0].mask
#m = a.coord.eq2top_m(0, aa.lat)
#top = top.transpose([1,0,2])
#x,y,z = n.dot(m, top)
#aa.select_chans([120])
#d = aa.ants[0].bm_response((x.flatten(),y.flatten(),z.flatten()), pol='y')[0]**2
#aa.select_chans(chans)
#d.shape = (DIM,DIM)
#bm_im = n.where(mask, 0, d)
#print 'done'

# Define a quick function writing an image to a FITS file
def fname(ftag, cnt): return '%s.%s.fits' % (opts.fmt % cnt, ftag)
def to_fits(ftag,i,src,cnt,history=''):
    filename = fname(ftag,cnt)
    print 'Saving data to', filename
    while len(i.shape) < 4: i.shape = i.shape + (1,)
    cen = ephem.Equatorial(src.ra, src.dec, epoch=aa.epoch)
    # We precess the coordinates of the center of the image here to
    # J2000, just to have a well-defined epoch for them.  For image coords to
    # be accurately reconstructed, precession needs to be applied per pixel
    # and not just per phase-center because ra/dec axes aren't necessarily
    # aligned between epochs.  When reading these images, to be 100% accurate,
    # one should precess the ra/dec coordinates back to the date of the
    # observation, infer the coordinates of all the pixels, and then
    # precess the coordinates for each pixel independently.
    cen = ephem.Equatorial(cen, epoch=ephem.J2000)
    a.img.to_fits(filename, i, clobber=True,
        object=src.src_name, obs_date=str(aa.date),
        ra=cen.ra*a.img.rad2deg, dec=cen.dec*a.img.rad2deg, epoch=2000.,
        d_ra=L[-1,-1]*a.img.rad2deg, d_dec=M[1,1]*a.img.rad2deg,
        freq=n.average(aa[0].beam.afreqs),history=history)

def grid_it(im,us,vs,ws,ds,wgts):
    #print 'Gridding %d integrations' % n_ints
    sys.stdout.write('|'); sys.stdout.flush()
    if len(ds) == 0: raise ValueError('No data to use.')
    ds,wgts = n.concatenate(ds), n.concatenate(wgts).flatten()
    us,vs,ws = n.concatenate(us), n.concatenate(vs), n.concatenate(ws)
    # Grid data into UV matrix
    (us,vs,ws),ds,wgts = im.append_hermitian((us,vs,ws),ds,wgts)
    im.put((us,vs,ws), ds, wgts)
    #im.put((us,vs,ws), ds, wgts, invker2=bm_im)

def img_it(im):
    global n_ints
    #print 'Imaging with %d integrations' % n_ints
    n_ints = 0
    # Form dirty images/beams
    uvs = a.img.recenter(n.abs(im.uv).astype(n.float), (DIM/2,DIM/2))
    bms = a.img.recenter(n.abs(im.bm[0]).astype(n.float), (DIM/2,DIM/2))
    dim = im.image((DIM/2, DIM/2))
    dbm = im.bm_image(term=0, center=(DIM/2,DIM/2))
    return uvs,bms, dim,dbm

# Loop through all specified sources, generating images
imgcnt = 0
for srccnt, s in enumerate(cat.values()):
    s.compute(aa)
    print '%d / %d' % (srccnt + 1, len(cat.values()))
    print 'Pointing (ra, dec):', s.ra, s.dec
    src = a.fit.SrcCatalog([s])
    # Gather data
    # Read each file
    for filename in args:
      sys.stdout.write('.'); sys.stdout.flush()

      if opts.bgsub:
          # grab background visibilities before integration
          ds_bg1,wgts_bg1 = [],[]
          uv = a.miriad.UV(filename)
          (j,t,j),j = uv.read()
          aa_bg1 = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
          aa_bg1.select_chans(chans)
          aa_bg1.set_active_pol(opts.pol)
          aa_bg1.set_jultime(t)
          src.compute(aa_bg1)
          if s.alt < opts.altmin * a.img.deg2rad: continue
          s_eq = src.get_crds('eq', ncrd=3)
          aa_bg1.sim_cache(s_eq)
          a.scripting.uv_selector(uv, opts.ant, opts.pol)
          uv.select('decimate', opts.decimate, opts.decphs-1)    # sets data to iterate over by selecting subset of integrations
          for (crd,t,(i,j)),d,f in uv.all(raw=True):
              pol = a.miriad.pol2str[uv['pol']]
              d,f = d.take(chans), f.take(chans)
              if not opts.skip_amp: d /= aa_bg1.passband(i,j)
          # Throws PointingError if not up:
              if not opts.skip_phs: d = aa_bg1.phs2src(d, s, i, j)
              u,v,w = aa_bg1.gen_uvw(i,j,src=s)
              longenough = n.where(n.sqrt(u**2+v**2) < opts.minuv, 0, 1).squeeze()
              if not opts.skip_bm:
                  # Calculate beam strength for weighting purposes
                  wgt = aa_bg1.bm_response(i,j).squeeze()
                  # Optimal SNR: down-weight beam-attenuated data 
                  # by another factor of the beam response.
                  d *= wgt; wgt *= wgt
              else: wgt = n.ones(d.shape, dtype=n.float)
              valid = n.logical_and(n.logical_not(f), longenough)
              d = d.compress(valid)
              if len(d) == 0: continue
              ds_bg1.append(d)
              wgts_bg1.append(wgt.compress(valid))

          # grab background visibilities after integration
          ds_bg2,wgts_bg2 = [],[]
          uv = a.miriad.UV(filename)
          (j,t,j),j = uv.read()
          aa_bg2 = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
          aa_bg2.select_chans(chans)
          aa_bg2.set_active_pol(opts.pol)
          aa_bg2.set_jultime(t)
          src.compute(aa_bg2)
          if s.alt < opts.altmin * a.img.deg2rad: continue
          s_eq = src.get_crds('eq', ncrd=3)
          aa_bg2.sim_cache(s_eq)
          a.scripting.uv_selector(uv, opts.ant, opts.pol)
          uv.select('decimate', opts.decimate, opts.decphs-1)    # sets data to iterate over by selecting subset of integrations
          for (crd,t,(i,j)),d,f in uv.all(raw=True):
              pol = a.miriad.pol2str[uv['pol']]
              d,f = d.take(chans), f.take(chans)
              if not opts.skip_amp: d /= aa_bg2.passband(i,j)
              # Throws PointingError if not up:
              if not opts.skip_phs: d = aa_bg2.phs2src(d, s, i, j)
              u,v,w = aa_bg2.gen_uvw(i,j,src=s)
              longenough = n.where(n.sqrt(u**2+v**2) < opts.minuv, 0, 1).squeeze()
              if not opts.skip_bm:
                  # Calculate beam strength for weighting purposes
                  wgt = aa_bg2.bm_response(i,j).squeeze()
                  # Optimal SNR: down-weight beam-attenuated data 
                  # by another factor of the beam response.
                  d *= wgt; wgt *= wgt
              else: wgt = n.ones(d.shape, dtype=n.float)
              valid = n.logical_and(n.logical_not(f), longenough)
              d = d.compress(valid)
              if len(d) == 0: continue
              ds_bg2.append(d)
              wgts_bg2.append(wgt.compress(valid))

          ds_bg = (n.array(ds_bg1) + n.array(ds_bg2))/2.

      uv = a.miriad.UV(filename)
      a.scripting.uv_selector(uv, opts.ant, opts.pol)
      uv.select('decimate', opts.decimate, opts.decphs)    # sets data to iterate over by selecting subset of integrations
      # Read all data from each file
      for (crd,t,(i,j)),d,f in uv.all(raw=True):
          pol = a.miriad.pol2str[uv['pol']]
          history = uv['history']
          history = history +  sys.argv[0].split('/')[-1].strip()+' ' + ' '.join(sys.argv[1:])
          aa.set_jultime(t)
          src.compute(aa)
          if s.alt < opts.altmin * a.img.deg2rad: continue
          s_eq = src.get_crds('eq', ncrd=3)
          aa.sim_cache(s_eq)
          aa.set_active_pol(pol)
          d,f = d.take(chans), f.take(chans)
          if not opts.skip_amp: d /= aa.passband(i,j)
          # Throws PointingError if not up:
          if not opts.skip_phs: d = aa.phs2src(d, s, i, j)
          u,v,w = aa.gen_uvw(i,j,src=s)
          longenough = n.where(n.sqrt(u**2+v**2) < opts.minuv, 0, 1).squeeze()
          if not opts.skip_bm:
              # Calculate beam strength for weighting purposes
              wgt = aa.bm_response(i,j).squeeze()
              # Optimal SNR: down-weight beam-attenuated data 
              # by another factor of the beam response.
              d *= wgt; wgt *= wgt
          else: wgt = n.ones(d.shape, dtype=n.float)
          valid = n.logical_and(n.logical_not(f), longenough)
          d = d.compress(valid)
          if len(d) == 0: continue
          n_ints += 1
          ds.append(d)
          us.append(u.compress(valid))
          vs.append(v.compress(valid))
          ws.append(w.compress(valid))
          wgts.append(wgt.compress(valid))
          # If data buffer is full, grid data
          if len(ds) * len(chans) > opts.buf_thresh:
              if opts.bgsub:
                  grid_it(im,us,vs,ws,ds-ds_bg,wgts)
              else:
                  grid_it(im,us,vs,ws,ds,wgts)
              us,vs,ws,ds,wgts = [],[],[],[],[]

    # Grid remaining data into UV matrix
    try:
        if opts.bgsub:
            grid_it(im,us,vs,ws,ds-ds_bg,wgts)
        else:
            grid_it(im,us,vs,ws,ds,wgts)
        uvs,bms,dim,dbm = img_it(im)
    except(ValueError):
        print 'No data: skipping output file.'
        continue
    for k in ['uvs','bms','dim','dbm']:
        if k in outputs: to_fits(k, eval(k), s, imgcnt,history=history)
    imgcnt += 1
    us,vs,ws,ds,wgts = [],[],[],[],[]
    if opts.no_w:
        im = a.img.Img(opts.size, opts.res, mf_order=0)
    else:
        im = a.img.ImgW(opts.size, opts.res, mf_order=0, wres=opts.wres)
    

