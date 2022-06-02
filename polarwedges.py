## extensively uses https://stackoverflow.com/questions/54395307/how-to-add-a-wedge-sector-onto-a-polar-matplotlib-plot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
import matplotlib.cm as cm
from matplotlib import gridspec

def perp( a ) :
    ##from https://stackoverflow.com/a/3252222/2454357
    b = np.empty_like(a)
    b[0] = -a[1]
    b[1] = a[0]
    return b

def seq_intersect(a1,a2, b1,b2) :
    ##from https://stackoverflow.com/a/3252222/2454357
    da = a2-a1
    db = b2-b1
    dp = a1-b1
    dap = perp(da)
    denom = np.dot( dap, db)
    num = np.dot( dap, dp )
    return (num / denom.astype(float))*db + b1

def angle(a1, a2, b1, b2):
    ##from https://stackoverflow.com/a/16544330/2454357
    x1, y1 = a2-a1
    x2, y2 = b2-b1
    dot = x1*x2 + y1*y2      # dot product between [x1, y1] and [x2, y2]
    det = x1*y2 - y1*x2      # determinant
    return np.arctan2(det, dot)  # atan2(y, x) or atan2(sin, cos)

def draw_cells(ax, cellsize, cmap, bgcolor):
    ## load earthquake data, note the azimuth is in radians
    indata = np.loadtxt('eq.txt', dtype={'names': ('Magnitude','Depth','Distance','Azimuth (rad)'), 'formats': ('f','f','f','f')}, delimiter=',')    
    m,d,R,T = list(zip(*indata))
    ## check that the cell sizes work and adjust if necessary
    if cellsize > 180:
        cellsize = 180
        print("Rescaling cellsize to " + str(cellsize))
    elif 180%cellsize != 0:
    ## rescaling cellsize to create no fractional cells
        cellsize = int(180/(180//cellsize))
        print("Rescaling cellsize to " + str(cellsize))
    ## set up lists for the count of quakes and colours    
    outlist=[[0 for y in range(360//cellsize)] for x in range(180//cellsize)]
    colors=[[0 for y in range(360//cellsize)] for x in range(180//cellsize)]
    ## process the earthquake data to get counts of quakes for each cell
    maxcolor = 0
    for eq in indata:
        x = int(eq[2]//cellsize)
        temp = eq[3]*180/3.14159-90
        y = int((temp%360)//cellsize)
        outlist[x][y] += 1
        if outlist[x][y]>maxcolor:
            maxcolor = outlist[x][y]
    print("The largest number of quakes in any cell is " + str(maxcolor))
    
    ## set up the color for each cell based on the earthquake count
    for r in range(180//cellsize):
        for t in range(360//cellsize):
            if outlist[r][t]==0:
    ## background color for no data 
                colors[r][t]=bgcolor
            else:
    ## this uses a logarithmic scale. For a linear scale, use outlist[r][t]/maxcolor
                colors[r][t]=cmap(np.log(outlist[r][t])/(np.log(maxcolor)))
    
    ##compute the corner points of the wedge:
    for r in range(180//cellsize):
        for t in range(360//cellsize):
            axtmin = 0
            r_min = (r+1)*cellsize
            r_max = (r)*cellsize
            t_min = (t*cellsize+cellsize)/180*np.pi
            t_max = t*cellsize/180*np.pi
            axtmin = 0
            rs = np.array([r_min,  r_max,  r_min, r_max, r_min, r_max])
            ts = np.array([axtmin, axtmin, t_min, t_min, t_max, t_max])
            ## from https://matplotlib.org/users/transforms_tutorial.html
            trans = ax.transData + ax.transAxes.inverted()
            ## convert to figure cordinates
            xax, yax = trans.transform([(t,r) for t,r in zip(ts, rs)]).T
            ## compute the angles of the wedge:
            tstart = np.rad2deg(angle(*np.array((xax[[0,1,2,3]],yax[[0,1,2,3]])).T))
            tend = np.rad2deg(angle(*np.array((xax[[0,1,4,5]],yax[[0,1,4,5]])).T))
            ## the center is where the two wedge sides cross (maybe outside the axes)
            center=seq_intersect(*np.array((xax[[2,3,4,5]],yax[[2,3,4,5]])).T)
            ## compute the inner and outer radii of the wedge:
            rinner = np.sqrt((xax[1]-center[0])**2+(yax[1]-center[1])**2)
            router = np.sqrt((xax[2]-center[0])**2+(yax[2]-center[1])**2)
            ## show the colored cell, alpha 0.8 means you can also see the data points
            wedge = Wedge(center, router, tstart, tend, width=router-rinner,
                          transform=ax.transAxes, lw=3,
                          facecolor=colors[r][t], alpha=0.8)
            ax.add_artist(wedge)
    ## plot the earthquake epicentres
    ax.scatter(T,R, c='blue', alpha=1.0)
    return maxcolor

if __name__ == "__main__":
    ## change the cellsize in degrees. Must be a factor of 180
    cellsize=20
    ## pick a matplotlib colormap: viridis, plasma, inferno, magma and cividis are perceptually uniform
    ## see: https://matplotlib.org/stable/tutorials/colors/colormaps.html
    cmap = getattr(cm, "plasma")
    ## choose a color for the background of the plot #FFFFFF = white, #000000 = black, #FF0000 = red, #00FF00 = green
    bgcolor = "#001133"
    
    ## set up the figure
    fig = plt.figure(figsize=(10,8))
    ## include a space for a customisable colorbar to the right of the polar plot
    spec = gridspec.GridSpec(ncols=2, nrows=1, width_ratios=[20,1])
    ## add a main title
    titletext = 'Earthquake epicentres detected from Cornwall, UK\nby back azimuth and angular distance.\nColours show count of quakes in ' + str(cellsize) + 'x' + str(cellsize) + ' wedges.'
    fig.suptitle(titletext, fontsize=16, y=1)
    
    ## set up the polar plot on ax1
    ax1=fig.add_subplot(spec[0], projection='polar')
    ## the radius is angular distance up to 180 degrees.
    ax1.set_ylim([0,180])
    ## azimuth angle increases clockwise, like a compass.
    ax1.set_theta_direction(-1)
    ## set the startpoint to be vertically up.
    ax1.set_theta_offset(np.pi / 2.0)
    
    ## see https://stackoverflow.com/a/41823326/2454357
    fig.canvas.draw()
    ## read the earthquake data and draw the colored cells, return the maximum count, which will be color 1.0 on the colorbar
    maxcolor = draw_cells(ax1, cellsize, cmap, bgcolor)

    ## set up the colorbar on ax2
    ## parameters for the color image across the entire colormap
    a = np.outer(np.linspace(0, 1, 1000),np.ones(1000))
    ax2=fig.add_subplot(spec[1])
    extent = (0, 1, 0, 1)
    ax2.imshow(a, aspect='auto', cmap=cmap, origin="lower", extent=extent, interpolation="nearest")
    ## format the tickmarks
    ax2.yaxis.set_label_position("right")
    ax2.yaxis.tick_right()
    ## map the actual earthquake count onto the colorbar labels - this creates integer values on the colorbar y-axis
    i = 1
    yticks=[]
    ylabels=[]
    while i <= maxcolor:
        ylabels.append(i)
        yticks.append(np.log(i)/(np.log(maxcolor)))
        i=i*2
    ## show the ticks and labels for the data on the colorbar
    ax2.set_yticks(yticks)
    ax2.set_yticklabels(ylabels)
    ## switch off ticks on the x-axis of the colorbar
    ax2.xaxis.set_ticklabels([])
    ax2.axes.xaxis.set_visible(False)

    plt.show()