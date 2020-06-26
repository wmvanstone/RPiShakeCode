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
from obspy import UTCDateTime, Stream, read, read_inventory
from obspy.taup import TauPyModel
from obspy.geodetics.base import locations2degrees
from matplotlib.transforms import blended_transform_factory
from os import path
import matplotlib.pyplot as plt
from geopy.geocoders import Nominatim
import numpy as np
from matplotlib.cm import get_cmap
from obspy.geodetics import gps2dist_azimuth
DATA_PROVIDER = "RASPISHAKE"

client=Client(DATA_PROVIDER)

# Event details
URL = 'https://earthquake.usgs.gov/earthquakes/eventpage/us6000ah9t/executive'
EQNAME = 'M7.4 - 12km SSW of Santa María Zapotitlán'
EQLAT = 16.029
EQLON = -95.901
EQZ = 26.3
EQTIME = '2020-06-23 15:29:05'
FILE_STEM = 'Mexico-2020-06-23'
MAGNITUDE = "M7.4"

DECIMATION = 1 # decimation factor, can be up to 15 to keep 1 sample point in every 15, reduces memory usage on computer
# Things to change once for your station
NETWORK = 'AM'   # AM = RaspberryShake network
STATION = "R7FA5"  # Station code of local station to plot
STA_LAT = 50.2609  # Latitude of local station  
STA_LON = -5.0434  # Longitude of local station
CHANNEL = 'EHZ'  # channel to grab data for (e.g. EHZ, SHZ, EHE, EHN)
found_home = False
# configuration of the section plot
DURATION = 1600  # Length in seconds of data to plot after origin time
TDOWNLOAD = 1800 # data length to download and write to file
MIN_DIST = 0 # minimum distance for a seismometer in degrees
MAX_DIST = 180 # maximum distance in degrees
STEP = 1.0 # step in degrees between seismometers
ANGLE_DX = 5 # angle between tick marks on x-axis of section
PHASES = sorted(["P", "pP", "PP", "S", "Pdiff", "PKP", "PKIKP", "PcP", "ScP", "ScS", "PKiKP", "SKiKP", "SKP", "SKS", "sP"]) # list of phases for which to compute theoretical times
#PHASES = sorted(["P", "pP", "PP", "S", "PcP", "ScP", "ScS", "PKiKP", "pPKiKP", "SKiKP", "sP"]) # list of phases for which to compute theoretical times

PHASE_PLOT = "spots" # choose lines or spots for phases
DPI = 80 # dpi for plot
FIG_SIZE = (1920/DPI, 1080/DPI) # size of figure in inches. Slightly bigger than PLOT_SIZE
PLOT_SIZE = (FIG_SIZE[0]*DPI*0.75,FIG_SIZE[1]*DPI*0.75) # plot size in pixels with borders for labels
F1 = 0.8  # High-pass filter corner for seismometers up to 90 degrees
F2 = 2  # Low-pass filter corner 
F3 = 0.8  # High-pass filter corner for seismometers from 90 degrees
F4 = 2  # Low-pass filter corner 
MODEL = 'iasp91'  # Velocity model to predict travel-times through
#MODEL = 'ak135'  # Velocity model to predict travel-times through

EXCLUDE=['RB305','R25AC','R64DA','RFD01','R25AC','R79D5','RE81A','R7363','R04C9','RCD29','RA71B','R0BEF','R50BE','RD2D3','RBD5C']
EXCLUDE+=['R2683','R2323','R623F','R51F6','RC131','RBA56','R328C','R2370','R3F1B','R5FE9','RE1CC','RE684','RDE9F','RDEE3','R2FF2']
EXCLUDE+=['R2F34','RFB1A','R95B0','R359E','R5C4C','R5D53','RCE32', 'R10DB','R1E9D','RF212','RBF42','R03D7','R8D5C','S4458','R21C3']
EXCLUDE+=['RBF59','R2547','R24AE','RAA7F','R1033','R06C4','R8118','R77AA','R3CC7','RB7BB','R9AC0','R4186','R5A10','R1B4E','RF082']
EXCLUDE+=['RBFE8','RFCBE','R0128','RD897','R0ODC','R00DC','RO0DC', 'R8FE8','RA666','R8A6C','R21E0','R617F','RB822','R013B']
EXCLUDE+=['R7BC1','R7A15','REA0F','R4A43','R4F38','R62B3','RB9CC','RD184','RBD3F','RCB48','RC848','R00A4','RF7E5','R2E8D','RAC91']
EXCLUDE+=['RC574','R9CDF','RA877','R8C2E','R2942','R3DD0','R9627','RFE9F','R86F8','R3EE8','S63B4','RBE3F','S1822','RE4AA','R240F']
EXCLUDE+=['RDCA5','RFBE4','RBCF7','RB669','R1516','REFF7','RD886','R548C','RD29A', 'R02B7','RF356','R29F6','R3302','RBC34']
EXCLUDE+=['R9F18','R319F','R2716','R3C2D','R7078','RBE0A','R1BEF','R0D06','R198D','R84BE','RFD43','RBF0E','R87FD']
EXCLUDE+=['RF9F7','RA17E','R52CE','REAF3','R26B4','R21D1','R9F13','R7D6B','R8DE6','RE1D6','RD14A','RB5E8','R12EB','RA150','RED95']
EXCLUDE+=['R4601','R1EE2','R35BB','RAC7E','RBFC6','R78B8','RCEAB','R03EA','RD540','R38DC','R34FA','R0CC5','R96A1','R99D0','R77C3','R11F5']
EXCLUDE+=['R7182','RD666','R9F16','R7184','R05A5','RA4CA','RA77C','R4D07','R69FC','R4D07','R44DC','R76F0','R9081','R7F9F','RA4D3','RD3AF','RB9B4','R08CF','RC01C']
EXCLUDE+=['R5954','RBE49','R0176','R4B31','RE443','R4FB1','RB30F','R830F','R138B','R1388','R13BB','R530B','R20A1','R74D9','R914C','R6B6E']
EXCLUDE+=['RA0BA','R6F15','R026F','RED42','R40AC','R84CF','RDC44','R96DC','R21DC','SA211', 'R5010','RC28B','R5FB4','R10A9','R2310', 'R5DDF','RE1A7']
EXCLUDE+=['R7ED6','R3DC0','RCC68','R75C2','RB537','R5698','R2004','R17FC','R4277','R4DFA','R78EC','R87F8','R2A2B','R38D9','RD4B6','R1650','R1OA9','RF920','R47BE']
EXCLUDE+=['R5CD0','R976A','R7D6A','RFFD2','REFFB','REFF8','R923A','RAB7B','RB77C','R6B1E','R0318','R6C38','R7813','R44D3','R3868']
EXCLUDE+=['R1ADD','R6324','R349C','R00D1','RF737','RE8DB','R834D','R2472','R55F8','R3E8F','R3B68','RCCE6','R0FFA','RD0BD','R3950']
EXCLUDE+=['R6707','RAA90','RA14E','R0BBA','R60D1','RE583','REC8D']
NETWORK_DATA = "ShakeNetwork2020.csv" # data file containing seismometer names and locations, different format to 2019 data
GLOBE_PHASES = sorted([
    # Phase, distance
    ('P', 26),
    ('PP', 60),
    ('pP', 45),
    ('PcP', 80),
    ('PKIKP', 150),
    ('PKiKP', 100),
    ('S', 65),
    ('sP', -30),
    ('ScS', -60),
    ('SKS', -82),
    ('ScP', -40),
    ('Pdiff', -120),
    ('PKP', -160),
    ('SKiKP', -100),
    ('SKP', -140)
])
# Calculated constants
STA_DIST = locations2degrees(EQLAT, EQLON, STA_LAT, STA_LON) # distance to local station
EQLATLON = (EQLAT, EQLON)
BUFFER = 60 # time before and after plot for taper data
START_TIME=UTCDateTime(EQTIME)
END_TIME=START_TIME+DURATION
# Pretty paired colors. Reorder to have saturated colors first and remove
# some colors at the end. This cmap is compatible with obspy taup
cmap = get_cmap('Paired', lut=12)
COLORS = ['#%02x%02x%02x' % tuple(int(col * 255) for col in cmap(i)[:3]) for i in range(12)]
COLORS = COLORS[1:][::2][:-1] + COLORS[::2][:-1]
#COLORS = [ cm.plasma(x) for x in linspace(0, 0.8, len(PHASES)) ] # colours from 0.8-1.0 are not very visible
# End of parameters to define

# utility subroutines for data handling and plotting
def plottext(xtxt, ytxt, phase, vertical, horizontal, color, textlist):
    clash = True
    while clash == True:
        clash = False
        for i in range(len(textlist)):
            while textlist[i][0] > (xtxt - 3) and textlist[i][0] < (xtxt + 3) and textlist[i][1] > (ytxt - 14) and textlist[i][1] < (ytxt + 14):
                clash = True
                ytxt -= 1
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
        distance = locations2degrees(EQLAT, EQLON, float(station[2]), float(station[3]))
        if (distance <= MAX_DIST and distance >= MIN_DIST and station[0] not in EXCLUDE):
            seismometers.append([station[0], round(float(station[2]),4), round(float(station[3]),4), round(distance,4)])
seismometers.sort(key = lambda i: i[3]) # Contains 'Code', Latitude, Longitude, distance (degrees), sort by distance

# read in seismic traces
waveform = Stream() # set up a blank stream variable
loaded = []
dist = 0 # record of how far away the next seismometer is that we are looking for
readit = False # flag to say that we have successfully read the data
loaded_stations = [] # list of stations successfully loaded
filtertext1 = ""
filtertext2 = ""

geolocator = Nominatim(user_agent="Raspberry Shake section plotter") # tool for getting place names
fp = open(FILE_STEM+"_stations.txt", "w")
for station in seismometers:
    distAB, azAB, azBA = gps2dist_azimuth(station[1], station[2], EQLAT, EQLON)   # Station dist in m from epicentre
    if station[0] == STATION or ((not((station[3] > STA_DIST-STEP) and (station[3] < STA_DIST+STEP))) and station[3] >= dist):
        # read in Raspberry Shake
        # Use file if it is there
        try:
            fileflag = "N"
            infile = "../" + FILE_STEM + "/"+station[0]+".mseed"
            invfile = "../" + FILE_STEM + "/"+station[0]+".xml"
            errfile = "../"+FILE_STEM+"/"+station[0]+".error"
            if path.isfile(infile) and path.isfile(invfile):
                if path.getsize(infile) > 0 and path.getsize(invfile) > 0:
                    st = read(infile)
                    inventory = read_inventory(invfile)
                    fileflag = "F"
        except:
            print("Failed to read file")
        # If file has not been loaded and there is no error file, download data
        if fileflag == "N" and not path.isfile(errfile):
            try:
            # Download and filter data
                st = client.get_waveforms(NETWORK, station[0], "00", CHANNEL, starttime=START_TIME-BUFFER, endtime=START_TIME+TDOWNLOAD+BUFFER,attach_response=True)
                inventory = read_inventory('https://fdsnws.raspberryshakedata.com/fdsnws/station/1/query?network=AM&station=%s&level=resp&format=xml&nodata=404&starttime=%s' % (station[0], str(START_TIME)))
                fileflag = "D"
                st.merge(method=0, fill_value='latest')
                st.write(infile, format='MSEED')
                inventory.write(invfile, format="STATIONXML")
            except:
                print("Failed to download trace " + station[0])
        # in the case that either a file has been loaded or trace downloaded                    
        if fileflag != "N":
            try:
                st.detrend(type='demean')
                pre_filt = [0.001, 0.005, 45, 50]
                st.remove_response(inventory=inventory, pre_filt=pre_filt,output="DISP",water_level=60,taper=True)#, plot=True)
                if station[3] <= 90:
                    st.filter('bandpass', freqmin=F1, freqmax=F2)
                    filtertext1 = "Bandpass Filter: freqmin="+str(F1)+", freqmax="+str(F2)+" at up to 90 degrees."
                else:
                    st.filter('bandpass', freqmin=F3, freqmax=F4)
                    filtertext2 = "Bandpass Filter: freqmin="+str(F3)+", freqmax="+str(F4)+" at greater than 90 degrees."
                st.decimate(DECIMATION)
                test = st.slice(START_TIME, END_TIME)
                if len(test[0]) > 0:
                    if station[0] == STATION: # make an additional copy to blacken plot
                        waveform += st.slice(START_TIME, END_TIME)    
                        loaded_stations.append(station)
                        loaded.append(station)
                        sta_x = len(loaded_stations)-1
                        found_home = True
                    waveform += st.slice(START_TIME, END_TIME)
                    outstring = str(station[1])+"\t"+str(station[2])+"\t"+str(station[3])+"\t"+str(abs(waveform[-1].max()))+"\t"+ station[0]+ "\t"+str( distAB)+ "\t"+str( azAB) +  "\t" +str( azBA) + "\n"                        
                    fp.write(outstring)
                    print(fileflag + "\t" + outstring, end="")
                    readit = True
            except:
                print("Failed to remove response " + station[0])
        else:
            readit = False
                
        if readit == False:
            ef = open(errfile, "w")
            ef.write("Error reading file")
            ef.close()
                
    if readit == True:
        # locate the seismometer and add this to the station record
        try:
            location = geolocator.reverse(str(station[1]) + ", " + str(station[2]))
            address_list = location.address.split(",")
            if len(address_list) > 4:
                address = ",".join(address_list[-1*(len(address_list)-2):]).strip() # remove the two most specific parts
            else:
                address = ",".join(address_list[:]).strip() # use the whole address
        except:
            address = "Unknown"
        station.append(address)
        loaded_stations.append(station)
        loaded.append(station) # Add the details of loaded seismometers to the loaded list
        readit = False
        if dist <= station[3]:
            dist = station[3] + STEP
            
fp.close()
if len(waveform) == len(loaded_stations):
# add station details to the waveforms and print out the details for interest
    opfile = open(nospaces(FILE_STEM)+'-Stations.txt', "w")
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
        #        opfile.write("--------------------------------------------------------------------------------------------------------\n")
        try:
            opfile.write(str(y) + " " + loaded_stations[y][0] + " Lat: " + str(loaded_stations[y][1]) + " Lon: " + str(loaded_stations[y][2]) + " Dist: " +
                str(loaded_stations[y][3]) + " degrees, Address: " + loaded_stations[y][4] + "\n")
        except:
            opfile.write(str(y) + " " + loaded_stations[y][0] + " Lat: " + str(loaded_stations[y][1]) + " Lon: " + str(loaded_stations[y][2]) + " Dist: " +
                str(loaded_stations[y][3]) + " degrees, Address: Not writable \n")
        print("--------------------------------------------------------------------------------------------------------")
        print(loaded_stations[y][0], "Lat:", loaded_stations[y][1], "Lon:", loaded_stations[y][2], "Dist:",
              loaded_stations[y][3], "degrees, Address:", loaded_stations[y][4])
        print("\n", waveform[y].stats, "\n")
    opfile.close()

# Create the section plot..
    fig = plt.figure(figsize=FIG_SIZE, dpi=DPI)
    plt.title('Section plot for '+FILE_STEM+' '+EQTIME+" lat:"+str(EQLAT)+" lon:"+str(EQLON), fontsize=12, y=1.07)
    # set up the plot area
    plt.xlabel("Angle (degrees)")
    plt.ylabel("Elapsed time (seconds)")
    plt.suptitle("Modelled arrival times for phases using " + MODEL + " with a focal depth of " + str(EQZ) + "km", fontsize=10)
    # plot the waveforms
    waveform.plot(size=PLOT_SIZE, type='section', recordlength=DURATION, linewidth=1.5, grid_linewidth=.5, show=False, fig=fig, color='black', method='full', starttime=START_TIME, plot_dx=ANGLE_DX, ev_coord = EQLATLON, dist_degree=True, alpha=0.50, time_down=True)
    ax = fig.axes[0]
    transform = blended_transform_factory(ax.transData, ax.transAxes)
    for tr in waveform:
        ax.text(float(tr.stats.distance), 1.0, tr.stats.station, rotation=270,
                va="bottom", ha="center", transform=transform, zorder=10, fontsize=8)
    # Add the local station name in red
    if found_home == True:
        ax.text(float(waveform[sta_x].stats.distance), 1.0, waveform[sta_x].stats.station, rotation=270,
            va="bottom", ha="center", transform=transform, zorder=10, fontsize=8, color = 'red')

    # print out the filters that have been used
    plt.text(0, DURATION *1.05, filtertext1)
    plt.text(0, DURATION *1.07, filtertext2)

# Print the coloured phases over the seismic section
    textlist = [] # list of text on plot, to avoid over-writing
    for j, phase in enumerate(PHASES):
        color = COLORS[PHASES.index(phase) % len(COLORS)]
        x=[]
        y=[]
        model = TauPyModel(model=MODEL)
        for dist in range(MIN_DIST, MAX_DIST+1, 1): # calculate and plot one point for each degree from 0-180
            arrivals = model.get_travel_times(source_depth_in_km=EQZ,distance_in_degree=dist, phase_list=[phase])
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
    ax1 = fig.add_axes([0.63, 0.48, 0.4, 0.4], polar=True)
    # Plot all pre-determined phases
    for phase, distance in GLOBE_PHASES:
        arrivals = model.get_ray_paths(EQZ, distance, phase_list=[phase])
        ax1 = arrivals.plot_rays(plot_type='spherical',
                                legend=False, label_arrivals=False,
                                plot_all=True, phase_list = PHASES,
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
    plt.savefig(nospaces(FILE_STEM)+'-Section.png')
    plt.show()

# report an error if the number of seismometers does not match the number of waveforms
else:
    print("Number of waveforms does not match number of seismometers, a seismometer may not have data in the required range")
    print(len(waveform), "waveforms", len(seismometers), "seismometers")