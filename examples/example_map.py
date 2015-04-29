import qpoint as qp
import numpy as np
import healpy as hp

# initialize, maybe change a few options from their defaults
Q = qp.QMap(nside=256, pol=True, accuracy='low',
            fast_math=True, mean_aber=True)

# dumb simulation
n = 100000
ctime = 1418662800. + np.arange(n)/100.
az = 100. + 40.*np.sin(2*np.pi*np.arange(n)/4000.)
el = 32. + 10.*np.mod(np.arange(n,dtype=float),500000.)/500000.
el = np.floor(el/0.1)*0.1
pitch = None
roll = None
lat = -77.6*np.ones_like(ctime)
lon = 165.7 - np.arange(n)*3/850000.

# step waveplate twice a day...
lmst = Q.lmst(ctime, lon)
hwp = np.ones_like(lmst)
hwp[lmst<=12] = 22.5
hwp[lmst>12] = 45.0
q_hwp = Q.hwp_quat(hwp)

# calculate boresight quaternions
q_bore = Q.azel2bore(az, el, pitch, roll, lon, lat, ctime)

# store pointing for mapmaking
Q.init_point(q_bore=q_bore, q_hwp=q_hwp)

cls = np.loadtxt('wmap7_r0p03_lensed_uK_ext.txt', unpack=True)
ell, cls = cls[0], cls[1:]
map_in = np.vstack(hp.synfast(cls, nside=1024, pol=True, new=True))
print map_in.mean(), map_in.sum()

# initialize source map
Q.init_source(map_in, pol=True)

# several detector offsets in degrees
delta_az_list = [-2,-1.0,-1.0,0,0,1.0,1.0,2];
delta_el_list = [2,-1.0,1.0,0,0,-1.0,1.0,-2];
delta_psi_list = [22.5,22.5,22.5,22.5,-22.5,-22.5,-22.5,-22.5];
q_off_list = Q.det_offset(delta_az_list, delta_el_list, delta_psi_list)

# run
print 'generate tod'
tod = Q.to_tod(q_off_list)
print tod.mean(), tod.sum()

# initialize and calculate hits and data maps
print 'bin to map'
smap, pmap = Q.from_tod(q_off_list, tod=tod)
print 'map stats'
print smap.mean(), smap.sum()
print pmap.mean(), pmap.sum()

print 'add to map'

# add noise
tod2 = 10 * np.random.randn(*tod.shape)

# update map, but don't double-count the hits map
Q.from_tod(q_off_list, tod=tod2, count_hits=False)
print 'map stats'
print smap.mean(), smap.sum()
print pmap.mean(), pmap.sum()

print 'plotting'

# extract columns
hits, p01, p02, p11, p12, p22 = pmap
m1, m2, m3 = smap

# better cartview
from spider_analysis.map.tools import cartview

opt = dict(lonra=[-170, -70], latra=[-50, -20], coord='C', cbar='h', cbar_size='8%')
lopt = dict(min=-300, max=300, **opt)

# plot stuff

cartview(map_in[0], title='input map', unit='Temperature [uK]', **lopt)
cartview(hits, title='hits', **opt)

mm = m1.copy()
mm[hits.astype(bool)] /= hits[hits.astype(bool)]
mm[~hits.astype(bool)] = np.nan
cartview(mm, title='normalized output map', unit='Temperature [uK]', **lopt)

import pylab
pylab.show()
