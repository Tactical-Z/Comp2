
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# for colored charts
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d import art3d
from matplotlib import cm
import matplotlib.collections as mcoll

# init starting constants ----
G = 6.674e-11 # gravitational constant
M = 6.972e23 # Earth mass
R = 6378 # Earth radius - km
J2 = 1.08262668e-3 # J2 element 
mu = 398600.4418
a = 6885 # Semi major axis - km
e = 0.01 # Eccentricity
anomoly = np.deg2rad(40) # true anomoly
n = np.sqrt(mu/a**3) # mean motion in rad/s
i = np.deg2rad(97.4) # Inclination
𝛺 = np.deg2rad(280) # RAAN 
𝜔 = np.deg2rad(40) # Argument of periapsis
theta_GMST = np.deg2rad(123.211)

cd = 2.2 # drag coefficient
A = 0.6*0.6 # drag area - m^2
m = 38 # mass - kg
rho0 = 3.614e-13 # atmospheric density - kg/m^3
h = a-R # orbital altitude - km
H = h / 8.5 # km 

# Times --
T = 2 * np.pi / n # orbital period in s
t_final = 100 * T # since we are simulating 100 orbits for c
dt = 10 # seconds
time = np.arange(0, t_final, dt)

# Anomalies ---
E0 = 2 * np.arctan(np.sqrt((1-e) / (1+e)) * np.tan(anomoly / 2)) # Eccentric anomaly
M0 = E0 -e * np.sin(E0) # mean anomaly # represents the fraction of an orbital period that has elapsed since pasing the periapsis # (comes from keplers equation) M = E − e sin E

# Functions ---

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

def eom(t, state):
    x,y,z,vx,vy,vz = state
    r_vec = np.array([x,y,z])
    r = np.linalg.norm(r_vec)

    # two body acceleration
    a_2b = -mu * r_vec / (r**3)

    # J2 pertubation
    factor = 3/2 * J2 * mu * (R**2) / (r**5)
    zx_ratio = (z**2)/(r**2)
    a_j2 = factor * np.array([
        x*(5*zx_ratio -1),
        y*(5*zx_ratio -1),
        z*(5*zx_ratio -3)
    ])

    # Atmospheric Drag
    rho = rho0 * np.exp(-h/H)
    v_vec = np.array([vx,vy,vz])
    v_mag = np.linalg.norm(v_vec)
    a_drag = -0.5 * cd * A/m * rho* v_mag* v_vec

    a_total = a_2b + a_j2 + a_drag
    return [vx,vy,vz, a_total[0], a_total[1], a_total[2]]


# init list for orbit values
Q = rotationZ(𝛺) @ rotationX(i) @ rotationZ(𝜔)
x_vals = []
y_vals = []
z_vals = []

# solve ----------
p = a*(1-e**2)
v_r = np.sqrt(mu/p) * e * np.sin(anomoly)          # radial velocity
v_theta = np.sqrt(mu/p) * (1 + e*np.cos(anomoly))  # transverse velocity
v_pqw = np.array([v_r, v_theta, 0])

r = a*(1-e**2)/(1 + e*np.cos(anomoly))  # radius at true anomaly
r0 = Q @ np.array([r, 0, 0])            # PQW vector along x-axis in PQW frame
v0 = Q @ v_pqw

state_0 = np.concatenate((r0,v0))

sol_J2 = solve_ivp(eom, (0,t_final), state_0, t_eval=time, rtol=1e-9)


# Convert to both solutions from ECI to ECEF
R_ECI2ECEF = np.array([
    [np.cos(theta_GMST),  np.sin(theta_GMST), 0],
    [-np.sin(theta_GMST), np.cos(theta_GMST), 0],
    [0,                   0,                  1]
])

eci_coordsJ2 = sol_J2.y[:3, :]  # Get first 3 rows (position)
ecef_coordsJ2 = R_ECI2ECEF @ eci_coordsJ2 # apply rotation to all points
x_ecefJ2, y_ecefJ2, z_ecefJ2 = ecef_coordsJ2 # make individual components for chart

# Earth Plots -------
# Creating 3D Earth sphere for chart
theta = np.linspace(0, 2*np.pi, 100)   # longitude
phi = np.linspace(0, np.pi, 100)      # latitude
theta, phi = np.meshgrid(theta, phi)
x_earth = R * np.sin(phi) * np.cos(theta)
y_earth = R * np.sin(phi) * np.sin(theta)
z_earth = R * np.cos(phi)

# Createing 2D earth circle for chart
earth_circle = plt.Circle((0, 0), R, color='b', alpha=0.2) # earth radius

# J2 pertubation plot cordinator
pointsJ2_2d = np.array([x_ecefJ2, y_ecefJ2]).T.reshape(-1, 1, 2)
segmentsJ2_2d = np.concatenate([pointsJ2_2d[:-1], pointsJ2_2d[1:]], axis=1)
normJ2 = plt.Normalize(0, len(time) - 1)
colorsJ2 = cm.viridis(normJ2(range(len(time)-1)))
lc_2dJ2 = mcoll.LineCollection(segmentsJ2_2d, colors=colorsJ2, linewidths=2)

# ------------------------------ J2 3D Plot ------------------------------
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# Chart orbit over color map rather than solid color 
# -----------------
points = np.array([sol_J2.y[0], sol_J2.y[1], sol_J2.y[2]]).T.reshape(-1,1,3)
segments = np.concatenate([points[:-1], points[1:]], axis=1)
norm = plt.Normalize(0, len(sol_J2.t))
colors = cm.viridis(norm(range(len(sol_J2.t)-1)))
lc = art3d.Line3DCollection(segments, colors=colors, linewidths=2)
ax.add_collection(lc)
# -----------------
ax.plot_surface(x_earth, y_earth, z_earth, color='b', alpha=0.2)
ax.plot([0, 0], [0, 0], [-R-100, R+100], color='k', linewidth=2)
ax.set_xlabel("x (km)")
ax.set_ylabel("y (km)")
ax.set_zlabel("z (km)")
ax.set_title("Satellite trajectory in ECI Frame including J2 pertubations and atmospheric drag")
ax.set_box_aspect([1,1,1])
plt.grid(True)
plt.show()
# ------------------------------ J2 3D Plot -----------------------------

# ------------------------------ Decay visualization ------------------------------
altitudes = np.linalg.norm(sol_J2.y[:3,:], axis=0) - R

orbits = int(t_final / T)
perigee_altitudes = []
for k in range(orbits):
    start_idx = int(k*T/dt)
    end_idx = int((k+1)*T/dt)
    perigee_altitudes.append(np.min(altitudes[start_idx:end_idx]))
time_perigee = np.arange(orbits) * T / 3600

plt.figure(figsize=(10,5))
plt.plot(time_perigee, perigee_altitudes, color='red', linewidth=2)
plt.xlabel("Time [hours]")
plt.ylabel("Perigee Altitude [km]")
plt.title("Satellite Perigee Decay due to Drag and J2")
plt.grid(True)
plt.show()

# ------------------------------ Decay visualization ------------------------------