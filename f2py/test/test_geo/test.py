from pygacode import geo
import numpy as np

np.set_printoptions(precision=16)

theta=np.array([0.1,0.2])

geo.geo_interp(theta,True)

print(geo.geo_b)
