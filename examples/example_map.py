import qpoint as qp
import numpy as np
import healpy as hp
from spider_analysis.map import synfast, cartview
import pylab

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
map_in = np.vstack(synfast(cls, nside=512, pol=True, new=True, deriv=0))
print map_in.mean(), map_in.sum()

# initialize source map
Q.init_source(map_in, pol=True)

# detector offsets
d_psi0 = [0, 45.2, 90.5, 135]
n = len(d_psi0)
d_az = [0]*n + [1]*n + [2]*n + [-1]*n + [-2]*n
d_el = [0]*n + [-1]*n + [-2]*n + [1]*n + [2]*n
d_psi = d_psi0*5
d_az, d_el = (d_az + d_az[n:],
              d_el + d_az[n:])
d_psi += d_psi[n:]

q_off_list = Q.det_offset(d_az, d_el, d_psi)

# run
print 'generate tod'
tod = Q.to_tod(q_off_list)
print tod.mean(), tod.sum()

# initialize and calculate hits and data maps
print 'bin to map'
vec, proj = Q.from_tod(q_off_list, tod=tod)
print 'map stats'
print vec.mean(), vec.sum()
print proj.mean(), proj.sum()

print 'add to map'

# add noise
tod2 = 10 * np.random.randn(*tod.shape)

# update map, but don't double-count the hits map
Q.from_tod(q_off_list, tod=tod2, count_hits=False)
print 'map stats'
print vec.mean(), vec.sum()
print proj.mean(), proj.sum()

print 'solving'

# solve
cond = Q.proj_cond()
print cond.min(), cond.max()
map_out = Q.solve_map(fill=np.nan)
print map_out.min(), map_out.max()

print 'plotting'

# plot stuff
opt = dict(lonra=[-170, -70], latra=[-50, -20], coord='C', cbar='h', cbar_size='8%')
lopt = dict(min=-300, max=300, **opt)

cartview(map_in[0], title='input map', unit='Temperature [uK]', **lopt)
cartview(proj[0], title='hits', **opt)
cartview(cond, title='condition number', **opt)
cartview(map_out[0], title='solved map', unit='Temperature[uK]', **lopt)

# difference map
md = map_out[0] - hp.ud_grade(map_in[0], 256)
md[map_out[0]==0] = np.inf
cartview(md, title='difference', unit='Temperature [uK]', min=-20, max=20, **opt)

pylab.show()
