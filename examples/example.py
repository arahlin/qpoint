import qpoint as qp
import numpy as np

# initialize, maybe set a few options
Q = qp.QPoint(accuracy='low', fast_math=True)

# update a bunch of parameters
# Q.set_params(daber='always',
#              lonlat='always',
#              wobble='never',
#              dut1='never',
#              erot='always',
#              npb=10,
#              aaber=100,
#              refro='never',
#              accuracy='low',
#              mean_aber=False,
#              fast_math=True,
#              polconv='healpix')

# dumb simulation
n = 10000
ctime = 1418662800. + np.arange(n)/100.
az = 80. + 100.*np.sin(2*np.pi*np.arange(n)/4000.)
el = 32. + 10.*np.mod(np.arange(n,dtype=float),500000.)/500000.
el = np.floor(el/0.1)*0.1
pitch = np.zeros_like(ctime)
roll = np.zeros_like(ctime)
lat = -77.6*np.ones_like(ctime)
lon = 165.7 - np.arange(n)*3/850000.

# calculate boresight quaternions
q = Q.azel2bore(az, el, pitch, roll, lon, lat, ctime)

# detector offset in degrees
delta_az = 1.0
delta_el = -1.0
delta_psi = 22.5

# calculate detector pointing
ra,dec,sin2psi,cos2psi = Q.bore2radec(delta_az, delta_el, delta_psi, ctime, q)
