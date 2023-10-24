from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
import math as m

import matplotlib.animation as anim

# Constants
G = 6.6743*10^-11

# All units in SI
r_lat_earth = 6357
sF_earth = 6378.14 / 6357
r_moon = 1737.4
d_moon = 384400
m_Earth = 5.972*10^24

# Given radius on the lateral axis, stretch factor, and phi, find r
def R_Oblate(r_lat, sF, phi):
    return np.sqrt((pow(r_lat, 2) * pow(r_lat*sF, 2)) / (pow(r_lat,2)*pow(np.sin(phi), 2) + pow(r_lat * sF, 2)*pow(np.cos(phi), 2)))

# Plot body
def Plot_Body(ax, body, color_name, zorder = 0, alpha = 1):
    ax.plot_surface(body.x, body.y, body.z, color = color_name, zorder = zorder, alpha = alpha)

# Plot launch site (r_eff extends effective radius to make point more visible)
def Plot_LS(ax, Launch_Site, color_name, size, r_eff = 1, zorder = 10, alpha = 1):
    ax.scatter(Launch_Site.x * r_eff, Launch_Site.y * r_eff, Launch_Site.z * r_eff, color = color_name, s=size, zorder = zorder, alpha = alpha)

# BODY CLASS
class Body:
    def __init__(self, mass, r_lat, stretch_factor, xWorld = 0, yWorld = 0, zWorld = 0):
        # Create mgrid for theta and phi values
        theta, phi = np.mgrid[0:2 * np.pi:120j, 0:np.pi:80j]
        # Define needed values to get a random point on the surface
        self.r_lat = r_lat
        self.stretch_factor = stretch_factor
        self.xPos = xWorld
        self.yPos = yWorld
        self.zPos = zWorld
        self.mass = mass
        # Use theta and phi values to create an object
        self.x = (R_Oblate(r_lat, stretch_factor, phi) * np.sin(theta) * np.sin(phi)) + xWorld
        self.y = (R_Oblate(r_lat, stretch_factor, phi) * np.cos(theta) * np.sin(phi)) + yWorld
        self.z = (R_Oblate(r_lat, stretch_factor, phi) * np.cos(phi)) + zWorld

    def Get_XYZ(self, phi, theta):
        x = (R_Oblate(self.r_lat, self.stretch_factor, phi) * np.sin(theta) * np.sin(phi)) + self.xPos
        y = (R_Oblate(self.r_lat, self.stretch_factor, phi) * np.cos(theta) * np.sin(phi)) + self.yPos
        z = (R_Oblate(self.r_lat, self.stretch_factor, phi) * np.cos(phi)) + self.zPos
        return (x, y, z)
    
    def Get_Mu(self):
        return G*self.mass

    def Get_Dist(self, point):
        return m.dist([self.xPos, self.yPos, self.zPos], point)

# LAUNCH CLASS
class LaunchSite:
    def __init__(self, Body, latitude, longitude):
        (self.x, self.y, self.z) = Body.Get_XYZ(latitude, longitude)
    
    def Get_XYZ(self):
        return [self.x, self.y, self.z]

# SATELLITE CLASS
class Satellite:
    def __init__(self, mass):
        self.mass = mass

# ORBIT CLASS
class Orbit:
    # Define an orbit given a launch site, deltaV, launch angle, and azimuth angle
    def __init__(self, Body, Launch_Site, deltaV, launchAngle, azimuth):
        # Force of gravity acting at launch site
        Fg = Body.Get_Mu*() / np.pow(Body.Get_Dist(Launch_Site.Get_XYZ), 2)


Earth = Body(m_Earth, r_lat_earth, sF_earth)
sat = Satellite(20000)
launch = LaunchSite(Earth, np.pi/4, np.pi/4)

fig = plt.figure()
ax = plt.axes(projection='3d')
Plot_Body(ax, Earth, 'green', alpha = 0.7, zorder = 0)
Plot_LS(ax, launch, 'red', 20, r_eff = 1.01, zorder = 10)
ax.set_xlim3d(-r_lat_earth * sF_earth, r_lat_earth * sF_earth)
ax.set_ylim3d(-r_lat_earth * sF_earth, r_lat_earth * sF_earth)
ax.set_zlim3d(-r_lat_earth * sF_earth, r_lat_earth * sF_earth)
plt.title("Orbit Analysis")
plt.show()