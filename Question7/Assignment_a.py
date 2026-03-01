
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

# Times --
T = 2 * np.pi / n # orbital period in s
t_final = 10 * T # since we are simulating 10 orbits for a 
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

# init list for orbit values
x_vals = []
y_vals = []

# itterate over time where t is each min (60 selk)
for t in time:
    M = M0 + n * t # Mean anomaly 
    E = solve_kepler(M,e) # Find local ecentric anomaly
    nu = 2 * np.arctan2(np.sqrt(1+e) * np.sin(E/2), np.sqrt(1-e)*np.cos(E/2)) # use ecentric anomaly to find true anomaly
    r = a * (1 - e*np.cos(E)) # find radius 
    x = r * np.cos(nu) # get x location from radius and true anomaly
    y = r * np.sin(nu) # get y location from radius and true anomaly
    x_vals.append(x)
    y_vals.append(y)

theta = np.linspace(0, 2*np.pi, 500)
x_earth = R * np.cos(theta)
y_earth = R * np.sin(theta)

# Plot
plt.figure()
plt.plot(x_earth, y_earth, 'b')
plt.plot(x_vals, y_vals, 'r')
plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title("Satellite Trajectory in Perifocal Frame (10 Orbits)")
plt.axis("equal")
plt.grid(True)
plt.show()