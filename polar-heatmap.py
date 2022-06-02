# this code is heavily dependent on this page: https://towardsdatascience.com/polar-heatmaps-in-python-with-matplotlib-d2a09610bc55
# I have abandoned this development pathway since I couldn't marry the plot to regular polar axes or colorbar
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection

def meshpolar():
    indata = np.loadtxt('eq.txt', dtype={'names': ('Magnitude','Depth','Latitude','Longitude'), 'formats': ('f','f','f','f')}, delimiter=',')
    m,d,r,theta = list(zip(*indata))
    patches = []
    eq_count = []
    cell_color = []
    
    outlist=[[0 for y in range(36)] for x in range(18)]
    
    for eq in indata:
        x = int(eq[2]//10)
        temp = 450 - eq[3]*180/3.14159
        y = int((temp%360)//10)
        outlist[x][y] += 1
    
    for x in range(17,-1,-1):
        for y in range(36):
            eq_count.append(outlist[x][y])
            wedge = mpatches.Wedge(0, x*10, y*10, (y+1)*10)
            patches.append(wedge)
            
    maxcolor = max(eq_count)
    print("The most earthquakes in a cell is " + str(maxcolor))
    
    for eq in eq_count:
        if eq == 0:
            cell_color.append('#FFFFFF')
        else:
            cell_color.append(cm.rainbow(eq/maxcolor))
    
    collection = PatchCollection(patches, linewidth=0.1,
        edgecolor=['#000000' for e in eq_count], facecolor=[c for c in cell_color])
            
    fig = plt.figure(figsize=(12,8), edgecolor='w',facecolor='w')
    
    fig.suptitle('Count of earthquake epicentres detected from Cornwall, UK\nby back azimuth and angular distance\ncoloured by count in 10x10 square.', fontsize=16, y=1)
    
    ax = fig.add_subplot()
    ax.add_collection(collection)
    plt.axis('equal')
    plt.axis('off')
    plt.tight_layout()
    
    plt.savefig('colormap.png')

    plt.show()

meshpolar()
