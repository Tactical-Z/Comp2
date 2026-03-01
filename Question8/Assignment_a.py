
import numpy as np

# for colored charts
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d import art3d
from matplotlib import cm
import matplotlib.collections as mcoll



a = 6378137.0 # semi major axis - m
f = 1/298.257223563 # flatteninng of earth 
e2 = 2*f - f**2 # eccentricity of earth

r_eci = np.array([5.6891, 1.9453, 3.4283]) * 1e6  # meters
theta = np.deg2rad(123.211)

R_ECI2ECEF = np.array([
    [np.cos(theta),  np.sin(theta), 0],
    [-np.sin(theta), np.cos(theta), 0],
    [0,              0,             1]
])

# execution ----

# Answere 1. (ECEF) 
r_ecef = R_ECI2ECEF @ r_eci
print("ECEF")
print(r_ecef)

# Answere 2. (Geocentric) 
x, y, z = r_ecef
r_mag = np.sqrt(x**2 + y**2 + z**2)
lat_gc = np.arctan2(z, np.sqrt(x**2 + y**2))
lon_gc = np.arctan2(y, x)
h_gc = r_mag - a   # spherical Earth assumption
print("Geocentric")
print("Latitude (deg):", np.rad2deg(lat_gc))
print("Longitude (deg):", np.rad2deg(lon_gc))
print("Altitude (m):", h_gc)


# Answere 3. (Geodetic) 
x, y, z = r_ecef
rho = np.sqrt(x**2 + y**2)

# Initial guess
phi = np.arctan2(z, rho)

for _ in range(5):
    N = a / np.sqrt(1 - e2 * np.sin(phi)**2)
    h = rho/np.cos(phi) - N
    phi = np.arctan2(z, rho*(1 - e2*N/(N+h)))

lon = np.arctan2(y, x)
print("Geodetic")
print("Latitude (deg):", np.rad2deg(phi))
print("Longitude (deg):", np.rad2deg(lon))
print("Altitude (m):", h)