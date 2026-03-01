
import numpy as np

a = 6378137.0 # semi major axis - m
f = 1/298.257223563 # flatteninng of earth 
e2 = 2*f - f**2 # eccentricity of earth
b = a*(1-f)

r_eci = np.array([5.6891, 1.9453, 3.4283]) * 1e6  # sat position ECI - meters
theta = np.deg2rad(123.211)
v_eci = np.array([-0.90839, 0.21127, -0.36081]) # line of sight vector

def RangeToIntesection(_r_eci, _v_eci):
    r_ecef = R_ECI2ECEF @ _r_eci
    v_ecef = R_ECI2ECEF @ _v_eci

    x0, y0, z0 = r_ecef
    vx, vy, vz = v_ecef

    # Quadratic coefficients
    A = (vx**2 + vy**2)/a**2 + (vz**2)/b**2
    B = 2*((x0*vx + y0*vy)/a**2 + (z0*vz)/b**2)
    C = (x0**2 + y0**2)/a**2 + (z0**2)/b**2 - 1

    # discriminant
    disc = B**2 - 4*A*C

    if disc < 0:
        print("No intersection")
        return -1
    else:
        s1 = (-B + np.sqrt(disc)) / (2*A)
        s2 = (-B - np.sqrt(disc)) / (2*A)
        s = min(val for val in [s1, s2] if val > 0)
        return s

R_ECI2ECEF = np.array([
    [np.cos(theta),  np.sin(theta), 0],
    [-np.sin(theta), np.cos(theta), 0],
    [0,              0,             1]
])



# execution ----
r_ecef = R_ECI2ECEF @ r_eci
v_ecef = R_ECI2ECEF @ v_eci

s = RangeToIntesection(r_eci, v_eci)

# Answere 2. ECEF intersection
r_int_ecef = r_ecef + s * v_ecef
print("\nIntersection ECEF (m):", r_int_ecef)

# Answere 1. ECI intersection
R_ECEF2ECI = R_ECI2ECEF.T
r_int_eci = R_ECEF2ECI @ r_int_ecef
print("\nIntersection ECI (m):", r_int_eci)

# Answere 3. Geocentric
x, y, z = r_int_ecef
lon_gc = np.arctan2(y, x)
lat_gc = np.arctan2(z, np.sqrt(x**2 + y**2))
print("\nGeocentric:")
print("Latitude (deg):", np.rad2deg(lat_gc))
print("Longitude (deg):", np.rad2deg(lon_gc))

# Answere 4. Geodetic
rho = np.sqrt(x**2 + y**2)
phi = np.arctan2(z, rho)

for _ in range(6):
    N = a / np.sqrt(1 - e2*np.sin(phi)**2)
    h = rho/np.cos(phi) - N
    phi = np.arctan2(z, rho*(1 - e2*N/(N+h)))

lon_gd = np.arctan2(y, x)

print("\nGeodetic:")
print("Latitude (deg):", np.rad2deg(phi))
print("Longitude (deg):", np.rad2deg(lon_gd))

# Answere 4. Altitude
print("\nAltitude at intersection (m):", h, "= 0")
