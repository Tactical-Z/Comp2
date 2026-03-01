
import numpy as np
import matplotlib.pyplot as plt

# init starting constants ----
G = 6.674e-11 # gravitational constant
M = 6.972e23 # Earth mass
R = 6378# Earth radius
mu = 398600.4418
a = 6885 # Semi major axis - km
e = 0.01 # Eccentricity
anomoly = np.deg2rad(40) # true anomoly
n = np.sqrt(mu/a**3) # mean motion in rad/s
i = np.deg2rad(97.4) # Inclination
𝛺 = np.deg2rad(280) # RAAN 
𝜔 = np.deg2rad(40) # Argument of periapsis

# Times --
T = 2 * np.pi / n # orbital period in s
t_final = 10 * T # since we are simulating 10 orbits for b as well
dt = 60 # seconds
time = np.arange(0,t_final,dt)

# Anomalies ---
E0 = 2 * np.arctan(np.sqrt((1-e) / (1+e)) * np.tan(anomoly / 2)) # Eccentric anomaly
M0 = E0 -e * np.sin(E0) # mean anomaly # represents the fraction of an orbital period that has elapsed since pasing the periapsis # (comes from keplers equation) M = E − e sin E

# Functions ---

# solve the keplers equation using newton raphson, returns the eccentric anomaly at mean anomaly input
# converges to the tolerance value for new E using f(E) and f'(E)
def solve_kepler(_meanAnomaly, _eccentricity, tolerance = 1e-10):
    E = _meanAnomaly # inital guess
    while True:
        f = E - _eccentricity * np.sin(E) - _meanAnomaly
        f_prime = 1 - _eccentricity*np.cos(E)
        E_new = E - f / f_prime
        if abs(E_new - E) < tolerance:
            return E_new
        E = E_new

def rotationZ(_angle):
    return np.array([
        [np.cos(_angle), -np.sin(_angle), 0],
        [np.sin(_angle), np.cos(_angle), 0],
        [0, 0, 1],
    ])

def rotationX(_angle):
    return np.array([
        [1, 0, 0],
        [0, np.cos(_angle), -np.sin(_angle)],
        [0, np.sin(_angle), np.cos(_angle)],
    ])

# init list for orbit values
Q = rotationZ(𝛺) @ rotationX(i) @ rotationZ(𝜔)
x_vals = []
y_vals = []
z_vals = []

# itterate over time where t is each min (60 sek)
for t in time:
    M = M0 + n * t # Mean anomaly 
    E = solve_kepler(M,e) # Find local ecentric anomaly
    nu = 2 * np.arctan2(np.sqrt(1+e) * np.sin(E/2), np.sqrt(1-e)*np.cos(E/2)) # use ecentric anomaly to find true anomaly
    r = a * (1 - e*np.cos(E)) # find radius 
    r_pqw = np.array([r*np.cos(nu), r*np.sin(nu), 0])
    r_eci = Q @ r_pqw

    x_vals.append(r_eci[0])
    y_vals.append(r_eci[1])
    z_vals.append(r_eci[2])

# Create spherical angles
theta = np.linspace(0, 2*np.pi, 100)   # longitude
phi = np.linspace(0, np.pi, 100)      # latitude

theta, phi = np.meshgrid(theta, phi)

# Sphere parametric equations
x_earth = R * np.sin(phi) * np.cos(theta)
y_earth = R * np.sin(phi) * np.sin(theta)
z_earth = R * np.cos(phi)

# Plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x_vals, y_vals, z_vals, 'r', linewidth=2)
ax.plot_surface(x_earth, y_earth, z_earth, color='b', alpha=0.2)
ax.plot([0, 0], [0, 0], [-R-100, R+100], color='k', linewidth=2)
ax.set_xlabel("x (km)")
ax.set_ylabel("y (km)")
ax.set_zlabel("z (km)")
ax.set_title("Satellite Trajectory in ECI Frame (10 Orbits)")
ax.set_box_aspect([1,1,1])
plt.grid(True)
plt.show()