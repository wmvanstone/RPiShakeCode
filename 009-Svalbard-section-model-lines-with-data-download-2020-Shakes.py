#!/usr/bin/env python

"""
Script to download and plot RaspberryShake station data
Also computes and plots theoretical phase arrival times

See https://docs.obspy.org/packages/obspy.taup.html for more info on 
Earth models and phase-naming nomenclature.

Copy this file and ShakeNetwork2020.csv into the same folder

Mark Vanstone
Truro School, Cornwall
Feb 2020
"""
from obspy.clients.fdsn import Client
from obspy import UTCDateTime, Stream, read
from obspy.taup import TauPyModel
from obspy.geodetics.base import locations2degrees
from matplotlib.transforms import blended_transform_factory
from os import path
from matplotlib import cm
from numpy import linspace
import matplotlib.pyplot as plt
from geopy.geocoders import Nominatim
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
client = Client('http://fdsnws.raspberryshakedata.com')

# Start of parameters to define
# Things that change for each earthquake event
LABEL = "M5.2 North of Svalbard" # Title to plot on figure
EVT_TIME = "2020-02-18T07:29:39"  # origin time of event
EVT_LAT = 85.740 # Latitude of event
EVT_LON = 22.699 # Longitude of event
EVT_Z = 10 # Depth of event in km
FILE_STEM = 'Svalbard' # folder name for saving mseed data files
DECIMATION = 1 # decimation factor, can be up to 15 to keep 1 sample point in every 15, reduces memory usage on computer
# Things to change once for your station
NETWORK = 'AM'   # AM = RaspberryShake network
STATION = "R5DDF"  # Station code of local station to plot
STA_LAT = 30.09009009  # Latitude of local station  
STA_LON = -95.445128  # Longitude of local station
CHANNEL = 'EHZ'  # channel to grab data for (e.g. EHZ, SHZ, EHE, EHN)
# configuration of the section plot
DURATION = 1000  # Length in seconds of data to plot after origin time
MIN_DIST = 0 # minimum distance for a seismometer in degrees
MAX_DIST = 90 # maximum distance in degrees
STEP = 0.5 # step in degrees between seismometers
ANGLE_DX = 10 # angle between tick marks on x-axis of section
PHASES = ["P", "pP", "PP", "S", "Pdiff", "PKP", "PKIKP", "PcP", "ScP", "ScS", "PKiKP", "SKiKP",
            "SKP", "SKS"] # list of phases for which to compute theoretical times
PHASE_PLOT = "spots" # choose lines or spots for phases
DPI = 80 # dpi for plot
FIG_SIZE = (1920/DPI, 1080/DPI) # size of figure in inches. Slightly bigger than PLOT_SIZE
PLOT_SIZE = (FIG_SIZE[0]*DPI*0.75,FIG_SIZE[1]*DPI*0.75) # plot size in pixels with borders for labels
F1 = 0.5  # High-pass filter corner for seismometers up to 90 degrees
F2 = 1.0  # Low-pass filter corner 
F3 = 0.5  # High-pass filter corner for seismometers from 90 degrees
F4 = 1.0  # Low-pass filter corner 
MODEL = 'iasp91'  # Velocity model to predict travel-times through
EXCLUDE = ['R6F29', 'R4355', 'R6324', 'RE063', 'RAE6A', 'REB59', 'R7143', 'R6F15', 'R49B6', 'RFF8B', 'R026F',
           'RBD93', 'R8D5C', 'R4FF5', 'RA211', 'RDD01', 'R2DB4', 'R16F7', 'RE8ED', 'REFF7', 'R9DAA', 'R6924',
           'RA482', 'REFFB', 'RCC45', 'R9627', 'R433E', 'R3EE8', 'R5DFC', 'R9F18', 'RA419', 'RA60A'] # noisy or mis-located seismometers
NETWORK_DATA = "ShakeNetwork2020.csv" # data file containing seismometer names and locations, different format to 2019 data
GLOBE_PHASES = [
    # Phase, distance
    ('P', 26), 
    ('PP', 60),
    ('p', 6),
    ('pP', 45),
    ('PcP', 80),
    ('PKIKP', 150),
    ('PKiKP', 100),
    ('S', 65),
    ('ScS', -60),
    ('SKS', -82),
    ('ScP', -40),
    ('PKIKP', -150),
    ('Pdiff', -120)
]
# Calculated constants
STA_DIST = locations2degrees(EVT_LAT, EVT_LON, STA_LAT, STA_LON) # distance to local station
EVT_LATLON = (EVT_LAT, EVT_LON)
START_TIME=UTCDateTime(EVT_TIME)
END_TIME=START_TIME+DURATION
COLORS = [ cm.plasma(x) for x in linspace(0, 0.8, len(PHASES)) ] # colours from 0.8-1.0 are not very visible
# End of parameters to define

# utility subroutines for data handling and plotting
def plottext(xtxt, ytxt, phase, vertical, horizontal, color, textlist):
    clash = True
    while clash == True:
        clash = False
        for i in range(len(textlist)):
            while textlist[i][0] > (xtxt - 3) and textlist[i][0] < (xtxt + 3) and textlist[i][1] > (ytxt - 24) and textlist[i][1] < (ytxt + 24):
                clash = True
                ytxt -= 2
    plt.text(xtxt, ytxt, phase, verticalalignment=vertical, horizontalalignment=horizontal, color=color, fontsize=10)
    textlist.append((xtxt, ytxt))

def parse(string):
    out = [""]
    counter = 0
    for l in string:
        if l.upper() in "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789.:-,":
            if ord(l) == 44:
                counter += 1
                out.append("")
            else:
                out[counter] += l
    return out

# Function to read a file using the context manager
def readFile(networkdata):
    # Read list of lines
    out = [] # list to save lines
    with open (networkdata, "r") as rd:
        # Read lines in loop
        for line in rd:
            # All lines (besides the last) will include  newline, so strip it
            out.append(parse(line.strip()))
    return out

def justnum(string):
    out = ""
    for l in string:
        if l in "0123456789":
            out += l
    return out

def nospaces(string):
    out = ""
    for l in string.upper():
        if l in "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789":
            out += l
        else:
            out += "_"
    return out

# Process seismic data
# Get the list of seismometers from file and sort them
allStations = readFile(NETWORK_DATA) # changed to work with 2020 stations list from http://www.fdsn.org/networks/detail/AM/
# work out the distance of each seismometer from the epicentre and sort the list by distance from event
seismometers = [] # empty list of seismometers
for station in allStations:
    if station[2] != "Latitude": # ignore the header line
        distance = locations2degrees(EVT_LAT, EVT_LON, float(station[2]), float(station[3]))
        if (distance <= MAX_DIST and distance >= MIN_DIST and station[0] not in EXCLUDE):
            seismometers.append([station[0], round(float(station[2]),4), round(float(station[3]),4), round(distance,4)])
seismometers.sort(key = lambda i: i[3]) # Contains 'Code', Latitude, Longitude, distance (degrees), sort by distance

# read in seismic traces
waveform = Stream() # set up a blank stream variable
dist = 0 # record of how far away the next seismometer is that we are looking for
readit = False # flag to say that we have successfully read the data
loaded_stations = [] # list of stations successfully loaded
filtertext1 = ""
filtertext2 = ""

geolocator = Nominatim(user_agent="Raspberry Shake section plotter") # tool for getting place names
for station in seismometers:
    if station[0] == STATION or ((not((station[3] > STA_DIST-STEP) and (station[3] < STA_DIST+STEP))) and station[3] >= dist):
        # read in Raspberry Shake
        if readit == False:
            try:
                # Download and filter data
                st = client.get_waveforms(NETWORK, station[0], "00", CHANNEL, starttime=START_TIME, endtime=START_TIME + DURATION)
                st.merge(method=0, fill_value='latest')
                st.detrend(type='demean')
                if station[3] <= 90:
                    st.filter('bandpass', freqmin=F1, freqmax=F2)
                    filtertext1 = "Bandpass Filter: freqmin="+str(F1)+", freqmax="+str(F2)+" at up to 90 degrees."
                else:
                    st.filter('bandpass', freqmin=F3, freqmax=F4)
                    filtertext2 = "Bandpass Filter: freqmin="+str(F3)+", freqmax="+str(F4)+" at greater than 90 degrees."
                st.decimate(DECIMATION)
                test = st.slice(START_TIME, END_TIME)
                if len(test) > 0:
                    if station[0] == STATION: # make an additional copy to blacken plot
                        waveform += st.slice(START_TIME, END_TIME)    
                        loaded_stations.append(station)
                        sta_x = len(loaded_stations)-1
                    waveform += st.slice(START_TIME, END_TIME)
                    readit = True
                else:
                    readit = False
            except:
                readit = False
    if readit == True:
        # locate the seismometer and add this to the station record
        try:
            location = geolocator.reverse(str(station[1]) + ", " + str(station[2]))
            address_list = location.address.split(",")
        except:
            location = "Unknown"
            address_list = location
        if len(address_list) > 4:
            address = ",".join(address_list[-1*(len(address_list)-2):]).strip() # remove the two most specific parts
        else:
            address = ",".join(address_list[:]).strip() # use the whole address
        station.append(address)
        loaded_stations.append(station)
        print(station[0], "Lat:", station[1], "Lon:", station[2], "Dist:", station[3], "degrees, Address:", station[4])
        # reset the search to look beyond the current station by STEP
        readit = False
        if dist <= station[3]:
            dist = station[3] + STEP

if len(waveform) == len(loaded_stations):
# add station details to the waveforms and print out the details for interest
    for y in range(len(waveform)):
        waveform[y].stats["coordinates"] = {} # add the coordinates to the dictionary, needed for the section plot
        waveform[y].stats["coordinates"]["latitude"] = loaded_stations[y][1]
        waveform[y].stats["coordinates"]["longitude"] = loaded_stations[y][2]
        waveform[y].stats["distance"] = loaded_stations[y][3]
        # Set the abbreviated name of the location in the network field for the plot title
        waveform[y].stats.network = NETWORK
        waveform[y].stats.station = loaded_stations[y][0]

        waveform[y].stats.location = loaded_stations[y][4]
        waveform[y].stats.channel = CHANNEL
        print("--------------------------------------------------------------------------------------------------------")
        print(loaded_stations[y][0], "Lat:", loaded_stations[y][1], "Lon:", loaded_stations[y][2], "Dist:",
              loaded_stations[y][3], "degrees, Address:", loaded_stations[y][4])
        print("\n", waveform[y].stats, "\n")

# Create the section plot..
    fig = plt.figure(figsize=FIG_SIZE, dpi=DPI)
    plt.title('Section plot for '+LABEL+' '+EVT_TIME+" lat:"+str(EVT_LAT)+" lon:"+str(EVT_LON), fontsize=12, y=1.07)
    # set up the plot area
    plt.xlabel("Angle (degrees)")
    plt.ylabel("Elapsed time (seconds)")
    plt.suptitle("Modelled arrival times for phases using " + MODEL + " with a focal depth of " + str(EVT_Z) + "km", fontsize=10)
    # plot the waveforms
    waveform.plot(size=PLOT_SIZE, type='section', recordlength=DURATION, linewidth=1.5, grid_linewidth=.5, show=False, fig=fig, color='black', method='full', starttime=START_TIME, plot_dx=ANGLE_DX, ev_coord = EVT_LATLON, dist_degree=True, alpha=0.50, time_down=True)
    ax = fig.axes[0]
    transform = blended_transform_factory(ax.transData, ax.transAxes)
    for tr in waveform:
        ax.text(float(tr.stats.distance), 1.0, tr.stats.station, rotation=270,
                va="bottom", ha="center", transform=transform, zorder=10, fontsize=8)
    # Add the local station name in red
    ax.text(float(waveform[sta_x].stats.distance), 1.0, waveform[sta_x].stats.station, rotation=270,
            va="bottom", ha="center", transform=transform, zorder=10, fontsize=8, color = 'red')

    # print out the filters that have been used
    plt.text(0, DURATION *1.05, filtertext1)
    plt.text(0, DURATION *1.07, filtertext2)

# Print the coloured phases over the seismic section
    textlist = [] # list of text on plot, to avoid over-writing
    for j, color in enumerate(COLORS):
        phase = PHASES[j]
        x=[]
        y=[]
        model = TauPyModel(model=MODEL)
        for dist in range(MIN_DIST, MAX_DIST+1, 1): # calculate and plot one point for each degree from 0-180
            arrivals = model.get_travel_times(source_depth_in_km=EVT_Z,distance_in_degree=dist, phase_list=[phase])
            printed = False
            for i in range(len(arrivals)):
                instring = str(arrivals[i])
                phaseline = instring.split(" ")
                if phaseline[0] == phase and printed == False and int(dist) > 0 and int(dist) < 180 and float(phaseline[4])>0 and float(phaseline[4])<DURATION:
                    x.append(int(dist))
                    y.append(float(phaseline[4]))
                    printed = True
            if PHASE_PLOT == "spots":
                plt.scatter(x, y, color=color, alpha=0.1, s=1)
            else:
                plt.plot(x, y, color=color, linewidth=0.3, linestyle='solid', alpha=0.1)
        if len(x) > 0 and len(y) > 0: # this function prevents text being overwritten
            plottext(x[0], y[0], phase, 'top', 'right', color, textlist)
            plottext(x[len(x)-1], y[len(y)-1], phase, 'top', 'left', color, textlist)
            
# Add the inset picture of the globe at x, y, width, height, as a fraction of the parent plot
    ax1 = fig.add_axes([0.63, 0.44, 0.4, 0.4], polar=True)
    # Plot all pre-determined phases
    for phase, distance in GLOBE_PHASES:
        arrivals = model.get_ray_paths(142.6, distance, phase_list=[phase])
        ax1 = arrivals.plot_rays(plot_type='spherical',
                                legend=False, label_arrivals=False,
                                plot_all=True,
                                show=False, ax=ax1)
        
    # Annotate regions of the globe
    ax1.text(0, 0, 'Solid\ninner\ncore',
            horizontalalignment='center', verticalalignment='center',
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.7))
    ocr = (model.model.radius_of_planet -
           (model.model.s_mod.v_mod.iocb_depth +
            model.model.s_mod.v_mod.cmb_depth) / 2)
    ax1.text(np.deg2rad(180), ocr, 'Fluid outer core',
            horizontalalignment='center',
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.7))
    mr = model.model.radius_of_planet - model.model.s_mod.v_mod.cmb_depth / 2
    ax1.text(np.deg2rad(180), mr, 'Solid mantle',
            horizontalalignment='center',
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.7))
    rad = model.model.radius_of_planet*1.15
    for phase in GLOBE_PHASES:
        ax1.text(np.deg2rad(phase[1]), rad, phase[0],
            horizontalalignment='center',
            bbox=dict(facecolor='white', edgecolor='none', alpha=0))

# Save the plot to file, then print to screen
    plt.savefig(nospaces(LABEL)+'-Section.png')
    plt.show()

# report an error if the number of seismometers does not match the number of waveforms
else:
    print("Number of waveforms does not match number of seismometers, a seismometer may not have data in the required range")
    print(len(waveform), "waveforms", len(seismometers), "seismometers")