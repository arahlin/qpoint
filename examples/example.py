import qpoint as qp
import numpy as np

# initialize, maybe change a few options from their defaults
Q = qp.QPoint(accuracy='low', fast_math=True, pair_dets=True)

# update a bunch of parameters
# Q.set(rate_daber='always',
#       rate_lonlat='always',
#       rate_wobble='never',
#       rate_dut1='never',
#       rate_erot='always',
#       rate_npb=10,
#       rate_aaber=100,
#       rate_ref='never',
#       accuracy='low',
#       mean_aber=False,
#       fast_math=True,
#       polconv='healpix')

# dumb simulation
n = 100000
ctime = 1418662800. + np.arange(n)/100.
az = 80. + 100.*np.sin(2*np.pi*np.arange(n)/4000.)
el = 32. + 10.*np.mod(np.arange(n,dtype=float),500000.)/500000.
el = np.floor(el/0.1)*0.1
pitch = None # np.zeros_like(ctime)
roll = None # np.zeros_like(ctime)
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

# detector offset in degrees
delta_az = 1.0
delta_el = -1.0
delta_psi = 22.5

q_off = Q.det_offset(delta_az, delta_el, delta_psi)

# calculate detector pointing
ra, dec, sin2psi, cos2psi = Q.bore2radec(q_off, ctime, q_bore, q_hwp=q_hwp)

# several detector offsets in degrees
delta_az_list = [-1.0,0,0,1.0];
delta_el_list = [1.0,0,0,-1.0];
delta_psi_list = [22.5,22.5,-22.5,-22.5];

q_off_list = Q.det_offset(delta_az_list, delta_el_list, delta_psi_list)

# noise
tod = 300 * np.random.randn(4,n)

# initialize and calculate hits and data maps
smap, pmap = Q.bore2map(q_off_list, ctime, q_bore, nside=256, q_hwp=q_hwp, tod=tod)

# several other detector offsets in degrees
delta_az_list2 = [-3.0,-2.0,2.0,3.0];
delta_el_list2 = [3.0,2.0,-2.0,-3.0];
delta_psi_list2 = [22.5,22.5,-22.5,-22.5];

q_off_list2 = Q.det_offset(delta_az_list2, delta_el_list2, delta_psi_list2)

tod2 = 300 * np.random.randn(4,n)

# update pmap
Q.bore2map(q_off_list2, ctime, q_bore, pmap=pmap, smap=smap, q_hwp=q_hwp, tod=tod2)

# extract columns
hits, p01, p02, p11, p12, p22 = pmap.T
m1, m2, m3 = smap.T

import healpy as hp
hp.mollview(hits)
hp.mollview(m1/hits,min=-300,max=300)
import pylab
pylab.show()
