from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
import math as m

from scipy import integrate

import central_mass as cm

import matplotlib.animation as anim

# Constants
G = 6.6743*pow(10, -11)

# Pertubations that are taken into account
is_J2 = True
is_solar_pres = False
is_third_body = True
is_atm_drag = True

# All units in SI
r_lat_earth = 6357
sF_earth = 6378.14 / 6357
r_moon = 1737.4
d_moon = 384400
m_Earth = 5.972*pow(10, 24)

# Given radius on the lateral axis, stretch factor, and phi, find r
def R_Oblate(r_lat, sF, phi):
    return np.sqrt((pow(r_lat, 2) * pow(r_lat*sF, 2)) / (pow(r_lat,2)*pow(np.sin(phi), 2) + pow(r_lat * sF, 2)*pow(np.cos(phi), 2)))


# This is a secondary body that has some orbit about the primary body. Does not impact other bodies gravitationally.
class Body:
    def __init__(self, name, mass, orbit, r_lat, stretch_factor = 1):
        # Name object
        self.name = name
        # The orbit
        self.orbit = orbit
        # Create mgrid for theta and phi values
        theta, phi = np.mgrid[0:2 * np.pi:120j, 0:np.pi:80j]
        # Define needed values to get a random point on the surface
        self.r_lat = r_lat
        self.stretch_factor = stretch_factor
        self.mass = mass
        # Use theta and phi values to create an object surface in cartesian from spherical
        if r_lat > 0:
            self.isPoint = False
            self.x = (R_Oblate(r_lat, stretch_factor, phi) * np.sin(theta) * np.sin(phi))
            self.y = (R_Oblate(r_lat, stretch_factor, phi) * np.cos(theta) * np.sin(phi))
            self.z = (R_Oblate(r_lat, stretch_factor, phi) * np.cos(phi))
        elif r_lat == 0:
            self.isPoint = True
        else:
            self.isPoint = True
            raise Exception("Object \'" + name + "\' has invalid radius value")
    
    def Plot(self, ax, color_name, zorder = 0, alpha = 1):
        center_pos = self.orbit.get_init_cond()
        if self.isPoint == False:
            ax.plot_surface(self.x + center_pos[0], self.y + center_pos[1], self.z + center_pos[2], color = color_name, zorder = zorder, alpha = alpha)
        else:
            ax.scatter(center_pos[0], center_pos[1], center_pos[2], color = color_name, s=10, zorder = zorder, alpha = alpha)
        self.orbit.Plot_Orbit(ax, color_name, zorder, alpha)

    def Get_Mu(self):
        return G*self.mass

    def Get_Dist(self, point):
        return m.dist([self.xPos, self.yPos, self.zPos], point)

# ORBIT CLASS: non-inertial object may have an orbit instance, which contains keplerian elements of the orbit and how to get cartesian elements from true anomaly
class Orbit:
    # Init using keplerian
    def __init__(self, a, e, inc, RAAN, argp, TA_epoch):
        # Assign parameters
        self.a = a
        self.e = e
        self.inc = inc
        self.RAAN = RAAN
        self.argp = argp
        # Get period
        self.period = (2*np.pi*m.pow(a, 3/2)) / cm.Get_Mu()
        # Define epoch true anomaly
        self.TA_epoch = TA_epoch

    # Get initial conditions
    def get_init_cond(self):
        return kep_2_cart(self.a, self.e, self.inc, self.RAAN, self.argp, self.TA_epoch)
    
    # Plot the orbit trajectory over a period
    def Plot_Orbit(self, ax, color_name, zorder, alpha):
        # Define x, y, z, numPoints vars
        numPoints = 200
        x = np.zeros(numPoints)
        y = np.zeros(numPoints)
        z = np.zeros(numPoints)
        # Run through each TA
        index = 0
        for i in np.linspace(self.TA_epoch, self.TA_epoch + 2*np.pi, numPoints):
            [x_it, y_it, z_it] = kep_2_cart(self.a, self.e, self.inc, self.RAAN, self.argp, i)[0:3]
            x[index] = x_it
            y[index] = y_it
            z[index] = z_it
            index += 1
        ax.plot(x, y, z, color=color_name, zorder = zorder, alpha = alpha)

    # Class method to get the orbit via cartesian inputs
    @classmethod       
    def from_cartesian(cls, mu, x, y, z, vx, vy, vz):
        [a, e, inc, RAAN, argp, TA] = cart_2_kep(mu, x, y, z, vx, vy, vz)
        return [a, e, inc, RAAN, argp, TA]

# Will convert cartesian to keplerian
def kep_2_cart(a, e, inc, RAAN, argp, TA): 
    Ecc_Anom = np.arccos((e + np.cos(TA)) / (1 + e*np.cos(TA)))
    # Express r and v in terms of radius, theta, local z direction
    r_mag = (a*(1-e*np.cos(Ecc_Anom)))
    r = r_mag * np.array([np.cos(TA), np.sin(TA), 0])
    v = (np.sqrt(cm.Get_Mu()*a) / r_mag) * np.array([-np.sin(Ecc_Anom), np.sqrt(1 - e**2) * np.cos(Ecc_Anom), 0])
    # Get rotation
    A = np.dot(np.dot(DCM(3, RAAN), DCM(1, inc)), DCM(3, argp))
    # Convert r and v to new cartesian coordinate system with matrices
    r_cart = np.dot(A,r)
    v_cart = np.dot(A,v)
    return [r_cart[0], r_cart[1], r_cart[2], v_cart[0], v_cart[1], v_cart[2]]

# Will convert keplerian to cartesian
def cart_2_kep(mu, x, y, z, vx, vy, vz):
    # radius
    r = np.sqrt(x**2 + y**2 + z**2)
    # angular momentum, its magnitude, and unit vector
    h = [y*vz - z*vy, vx*z - x*vz, x*vy - y*vx]
    mh = np.sqrt(h[1]**2 + h[2]**2 + h[3]**2)
    uh = h/mh
    # eccentricity vector and eccentricity value
    ev = [(vy*h[3] - vz*h[2])/mu - x/r, (vz*h[1] - vx*h[3])/mu - y/r, (vx*h[2] - vy*h[1])/mu - z/r]
    ecc = np.sqrt(ev(1)**2 + ev(2)**2 + ev(3)**2)
    ne = ev/ecc
    # semi-major axis
    a = (mh**2/mu)/(1 - ecc**2)
    # inclination
    inc = np.arccos(h[3]/mh)
    # line of nodes
    n = [-uh[2], uh[1], 0]
    n = n/np.sqrt(n[1]**2 + n[2]**2)
    # raan
    raan = np.arctan2(uh[1],-uh[2])
    # argp
    cosw = (ne[1]*n[1] + ne[2]*n[2])
    sinw = ((n[2]*ne[3]-n[3]*ne[2])*uh[1] + (n[3]*ne[1]-n[1]*ne[3])*uh[2] + (n[1]*ne[2]-n[2]*ne[1])*uh[3])
    argp = np.mod(np.arctan2(sinw,cosw), 2*np.pi)
    # true anomaly
    cosf = (ne[1]*x + ne[2]*y + ne[3]*z)/r
    sinf = ((ne[2]*z-ne[3]*y)*uh[1] + (ne[3]*x-ne[1]*z)*uh[2] + (ne[1]*y-ne[2]*x)*uh[3])/r
    f = np.arctan2(sinf,cosf)
    return [a, ecc, inc, raan, argp, f]

# Propogate an object's orbit
# def propogator(orbit):
#     r = np.sqrt(x**2 + y**2 + z**2)
#     d2x_dt = -cm.mu * (x / r**3)
#     d2y_dt = -cm.mu * (y / r**3)
#     d2z_dt = -cm.mu * (z / r**3)
#     if (is_J2):
#         d2x_dt += (1/2) * (cm.mu * cm.J2 * cm.r**2 * ((15*x*z**2 / r**7) - (3*x / (r**5))))
#         d2y_dt += (1/2) * (cm.mu * cm.J2 * cm.r**2 * ((15*y*z**2 / r**7) - (3*y / (r**5))))
#         d2z_dt += (1/2) * (cm.mu * cm.J2 * cm.r**2 * ((15*z**3 / r**7) - (9*z / (r**5))))

# DCM creator
def DCM(axis, angle):
    # Get rotational matrices
    if(axis == 1 or axis == 'x'):
        return [[1, 0, 0], [0, np.cos(angle), -np.sin(angle)], [0, np.sin(angle), np.cos(angle)]]
    elif (axis == 2 or axis == 'y'):
        return [[np.cos(angle), 0, np.sin(angle)], [0, 1, 0], [-np.sin(angle), 0, np.cos(angle)]]
    elif (axis == 3 or axis == 'z'):
        return [[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]]
    else:
        raise Exception(f'DCM:axis -- Variable {axis} out of range', axis)
   

Satellite = Body("Probe", 100, Orbit(7000, 0, np.radians(30), np.radians(30), np.radians(30), 0), 0)
fig = plt.figure()
ax = plt.axes(projection='3d', computed_zorder=False)
ax.set_aspect('equal', adjustable='box')
cm.Plot_CM(ax, "Green")
Satellite.Plot(ax, "Red")
#ax.set_xlim3d(-r_lat_earth * sF_earth, r_lat_earth * sF_earth)
#ax.set_ylim3d(-r_lat_earth * sF_earth, r_lat_earth * sF_earth)
#ax.set_zlim3d(-r_lat_earth * sF_earth, r_lat_earth * sF_earth)
plt.title("Orbit Analysis")
plt.show()

def SavePlot():
    plt.savefig('map.png')