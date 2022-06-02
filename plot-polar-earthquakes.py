import numpy as np
import matplotlib.pyplot as plt

def realpolar():
    indata = np.loadtxt('eq.txt', dtype={'names': ('Magnitude','Depth','Latitude','Longitude'), 'formats': ('f','f','f','f')}, delimiter=',')
    m,d,r,theta = list(zip(*indata))
    area = np.square(m)
    fig = plt.figure(figsize=(12,8))
    fig.suptitle('Earthquake epicentres detected from Cornwall, UK\nby back azimuth and angular distance, coloured by depth (km)\nsized by magnitude squared.', fontsize=16, y=1)
    ax = fig.add_subplot(projection='polar')
    ax.set_theta_direction(-1)
    ax.set_theta_offset(np.pi / 2.0)
    plt.ylim([0,180])
    plt.yticks(np.arange(0, 180, 20))
    ctf = ax.scatter(theta, r, c=d, s=area, cmap='rainbow', alpha=0.75)
    plt.colorbar(ctf)
    plt.savefig('eqbyazimuth.png')
    plt.show()
    
realpolar()