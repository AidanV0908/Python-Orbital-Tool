import numpy as np

# Constants
G = 6.6743*pow(10, -11)

# Central mass variables (Earth by default)
name = "Earth"
r = 6378.14
r_pol = 6357
m = 5.972*pow(10, 24)

J2 = 0.001082635

isPoint = False

# Given radius on the lateral axis, stretch factor, and phi, find r
def R_Oblate(r_lat, sF, phi):
    return np.sqrt((pow(r_lat, 2) * pow(r_lat*sF, 2)) / (pow(r_lat,2)*pow(np.sin(phi), 2) + pow(r_lat * sF, 2)*pow(np.cos(phi), 2)))

def redefine_central_mass(bName, bR_Lat, bR_pol, bM):
    # Name object
    name = bName
    # Create mgrid for theta and phi values
    theta, phi = np.mgrid[0:2 * np.pi:120j, 0:np.pi:80j]
    # Define needed values to get a random point on the surface
    r_lat = bR_Lat
    bSF = bSF
    m = bM
    # Use theta and phi values to create an object surface in cartesian from spherical
    if r_lat > 0:
        isPoint = False
    elif r_lat == 0:
        isPoint = True
    else:
        isPoint = True
        raise Exception("Object \'" + name + "\' has invalid radius value")
        
# Plot central mass
def Plot_CM(ax, color_name, zorder = 0, alpha = 1):
    # Stretch factor
    sf = r / r_pol
    # Create mgrid for theta and phi values
    theta, phi = np.mgrid[0:2 * np.pi:120j, 0:np.pi:80j]
    # x, y, z of objects
    x = (R_Oblate(r, sf, phi) * np.sin(theta) * np.sin(phi))
    y = (R_Oblate(r, sf, phi) * np.cos(theta) * np.sin(phi))
    z = (R_Oblate(r, sf, phi) * np.cos(phi))
    ax.plot_surface(x, y, z, color = color_name, zorder = zorder, alpha = alpha)

def Get_Mu():
    return G*m