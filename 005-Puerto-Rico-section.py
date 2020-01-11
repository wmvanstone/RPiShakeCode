# Designed to be run with Programs, Data and Plots subfolders under a Seismic folder
# Copy this file into the Programs folder
# Copy ShakeNetwork2019.csv into the Data folder
# Developed by Mark Vanstone using Thonny, a free Python IDE designed for new programmers
#
from datetime import datetime
from obspy import UTCDateTime, Stream
from obspy.geodetics import gps2dist_azimuth
from obspy.clients import fdsn, iris
import matplotlib.pyplot as plt
from matplotlib.transforms import blended_transform_factory
client = fdsn.Client()
irisclient = iris.Client()

# EQ details and paramaters for data selection and plotting
eqname = "M6.0 Puerto Rico"
eqlat = 17.8694
eqlon = -66.8088
eqlatlon=(eqlat, eqlon)
eqtime = "2020-01-11 12:54:45"
# Plot parameters
plots=['normal', 'section', 'distance'] # choose normal or section, section by distance or angle, also see sortkey below
sectiondx=1e5 # distance between tick marks on x-axis of distance section
angledx=2 # angle between tick marks on x-axis of section by angle
sortkey = 0 # 0 = sort by distance from epicentre, 2 = sort by azimuth, which is effective for the normal plot
# Vancouver Island excludes - noisy or geographically misplaced recorders
exclude=['RB293','R93B1','R7813','R7783','R5A78','RDCBA','R1E5E','R923A','RE650','R37BE','RB0B5','R3D81','R6392','RCD29','RCE32','R6324','R8473','R4FF5','R7C93','R7143','R21B0','REFFB','RCC6E','R1B4E','R7F9F','RA877','R976A','R13D2','R08CF','R17FC','R3F1B','R5C4C','R95B0','RE859','RA78B','R3F39','REB59','RF566','RF212','R78B7','R5768','R550B','R23FF', 'RC98E', 'R89A5', 'RDE9B','RDE9D', 'R6EE9', 'R49FD', 'REEAA', 'RA89C', 'RF71F', 'RD2A1', 'R0318', 'R0F1A', 'R79EC', 'R87F8','R70B6', 'R7363', 'RFC18','RABDC']
# peru excludes - noisy or geographically misplaced recorders
exclude+=['R440A','R8282','R16F7','R74D9','R1108','RA13A', 'R39DB']
maxdist = 2000000 # maximum distance in metres, ignore any seismometers beyond this distance
step = 15000 # step in metres between seismometers
# selecting by bearing allows a narrow window to be chosen by compass bearing
minbearing = 0 # this mechanism needs recoding, it does not work for degrees on either side of 360
maxbearing = 360
bearingstep = 2.5 # allows a spread of data in the normal plot
duration = 900 # duration of the record
timeoffset = 0 # time added before start of general plot
networkdata = "../Data/ShakeNetwork2019.csv"

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

def date2nthDay(date_string, format='%Y%m%d'):
    date_object = datetime.strptime(date_string, format)
    newYearDay = datetime(year=date_object.year, month = 1, day = 1)
    return(date_object - newYearDay).days + 1

def pad(num):
    out = str(num)
    while len(out) < 3:
        out = "0" + out
    return out

def bubblesort(ind, dist, col):
    # ignore the first line for the header row
    for i in range(len(ind)-1,0,-1):
        for j in range(1,i):
            if float(dist[ind[j]][col])>float(dist[ind[j+1]][col]):
                temp = ind[j]
                ind[j] = ind[j+1]
                ind[j+1] = temp
    return ind

# read station list from network data file
allStations = readFile(networkdata)
# remove duplicates and keep the latest in the list. Duplicates occur when stations have been moved. 
cutStations = []
for x in range(1, len(allStations)):
    if allStations[x][0] != allStations[x-1][0]:
        cutStations.append(allStations[x-1])
cutStations.append(allStations[x])
# set the data window
dayno = date2nthDay(justnum(eqtime[0:10]))
year=eqtime[0:4]
starttime=UTCDateTime(eqtime)
endtime=starttime+duration
# work out the distance of each seismometer from the epicentre
distance = [["Distance(m)", "Azimuth A-B", "Azimuth B-A"]]
for x in range(1, len(cutStations)):
    distance.append(gps2dist_azimuth(float(cutStations[x][4]), float(cutStations[x][5]), eqlat, eqlon))
index=[x for x in range(len(distance))]
index = bubblesort(index, distance, sortkey)# col 0=dist 1=azimuth A-B 2=azimuth B-A
# set up empty lists
seislist = []
seismometers = []
locations = []
latitudes = []
longitudes = []
distances = []
counter = 1
# run through the list and select stations that are within the required window - needs refining for bearings near 360
while counter < len(cutStations):
    x = index[counter]
    if distance[x][2] > minbearing and distance[x][2] < maxbearing and distance[x][0] <= maxdist and (cutStations[x][0] not in exclude):
        seislist.append(cutStations[x][0])
        locations.append(cutStations[x][0])
        latitudes.append(float(cutStations[x][4]))
        longitudes.append(float(cutStations[x][5]))
        distances.append(distance[x])
    counter += 1
# start reading the data
# read in seismometers
waveform = Stream() # set up a blank stream variable
readit = False
dist = 0
bearing = minbearing
for x in range (len(seislist)):
    print(x, int(distances[x][0]/1000), "km", round(distances[x][2],1), "degrees, looking for greater than", round(bearing, 1), "degrees")
    if (sortkey == 0 and distances[x][0] >= dist) or (sortkey == 2 and distances[x][2] >= bearing):        
        # read in Raspberry Shake
        client=fdsn.Client(base_url='https://fdsnws.raspberryshakedata.com/')
        readit = False
        sensorlist=["EHZ", "SHZ", "EHN", "EHE", "EH[Z,N,E]"]
        for sensor in sensorlist:
            if readit == False:
                try:
                    if len(seismometers) == 0:
                        copied = False
                        st=client.get_waveforms('AM',seislist[x],'00',sensor,starttime-timeoffset,endtime)
                        st.merge(method=0, fill_value='latest')
                        st.detrend(type='demean')
                        if distances[x][0] <= 5000000:
                            st.filter('bandpass', freqmin=0.1, freqmax=2)
                        else:
                            st.filter('bandpass', freqmin=0.1, freqmax=1)
                        print(len(st), locations[x])
                        if len(st) > 1:
                            for s in range(len(st)):
                                if starttime >= st[s].stats.starttime and starttime <= st[s].stats.endtime and copied == False:
                                    waveform = st[s].slice(starttime, endtime)
                                    copied = True
                            if copied == False:
                                for s in range(len(st)):
                                    if endtime > st[s].stats.starttime and endtime <= st[s].stats.endtime and copied == False:
                                        waveform = st[s].slice(starttime, endtime)
                                        copied = True
                        else:
                            waveform = st
                            copied = True
                        if copied == True:
                            seismometers.append(x)
                        readit = True
                    else:
                        copied = False
                        st=client.get_waveforms('AM',seislist[x],'00',sensor,starttime-timeoffset,endtime)
                        st.merge(method=0, fill_value='latest')
                        st.detrend(type='demean')
                        if distances[x][0] <= 5000000: # I have chosen to apply different filters depending on distance
                            st.filter('bandpass', freqmin=0.1, freqmax=2)
                        else:
                            st.filter('bandpass', freqmin=0.1, freqmax=1)
                        print(len(st), locations[x])
                        if len(st) > 1:
                            for s in range(len(st)):
                                if starttime >= st[s].stats.starttime and starttime <= st[s].stats.endtime and copied == False:
                                    waveform += st[s].slice(starttime, endtime)
                                    copied = True
                            if copied == False:
                                for s in range(len(st)):
                                    if endtime > st[s].stats.starttime and endtime <= st[s].stats.endtime and copied == False:
                                        waveform += st[s].slice(starttime, endtime)
                                        copied = True
                        else:
                            waveform += st
                            copied = True
                        if copied == True:
                            seismometers.append(x)
                        readit = True
                except:
                    print(seislist[x], sensor, "off line")
    if readit == True:
        readit = False
        if sortkey == 0:
            while dist <= distances[x][0]:
                dist = distances[x][0] + step
        elif sortkey == 2:
            while bearing <= distances[x][2]:
                bearing = float(distances[x][2]) + bearingstep

for y in range(len(seismometers)):
    waveform[y].stats["coordinates"] = {} # add the coordinates to the dictionary, needed for the section plot
    waveform[y].stats["coordinates"]["latitude"] = latitudes[seismometers[y]]
    waveform[y].stats["coordinates"]["longitude"] = longitudes[seismometers[y]]
    if 'angle' in plots: #put central angle in the distance
        result = round(irisclient.distaz(stalat=latitudes[seismometers[y]], stalon=longitudes[seismometers[y]], evtlat=eqlat, evtlon=eqlon)['distance']*1.0,1)
    else: # put distance in km in the distance attribute
        result = int(distances[seismometers[y]][0]/1000)
    waveform[y].stats["distance"] = str(result)
    # Set the abbreviated name of the location in the network field for the plot title
    waveform[y].stats.network = locations[seismometers[y]]
    waveform[y].stats.station = ""
    waveform[y].stats.location = "" #str(round(distances[seismometers[y]][0]/1000,1)) + "km"
    if 'angle' in plots:
        waveform[y].stats.channel = str(round(distances[seismometers[y]][2],1)) # Azimuth B->A
    else:
        waveform[y].stats.channel = str(int(distances[seismometers[y]][0]/1000)) # Distance in km
    print("\n", waveform[y].stats, "\n")

# very big data sets may need to be decimated, otherwise the plotting routine crashes
# if duration * len(waveform) > 6000:
#    waveform.decimate(10, strict_length=False, no_filter=True)
if 'normal' in plots:
    fig1 = plt.figure(figsize=(18, 13))
    #plt.title('Plot for epicentre near Taunton and rupture at UTC 22:49:19.5', fontsize=12, y=1.07)
    if sortkey == 0:
        plt.title('Plot for '+eqname+' sorted by distance. '+eqtime+" lat:"+str(eqlat)+" lon:"+str(eqlon), fontsize=12, y=1.07)
        plt.xticks([])
        plt.yticks([])
        waveform.plot(size=(1024,750),type='normal', automerge=False, equal_scale=False, fig=fig1, starttime=starttime-timeoffset, endtime=endtime-timeoffset, outfile='../Plots/' + nospaces(eqname) + '-normal-ByDistance.png', fontsize=10)
    elif sortkey == 2:
        plt.title('Plot for '+eqname+' sorted by azimuth. '+eqtime+" lat:"+str(eqlat)+" lon:"+str(eqlon), fontsize=12, y=1.07)
        plt.xticks([])
        plt.yticks([])
        waveform.plot(size=(1024,750),type='normal', automerge=False, equal_scale=False, fig=fig1, starttime=starttime-timeoffset, endtime=endtime-timeoffset, outfile='../Plots/' + nospaces(eqname) + '-normal-ByAzimuth.png', fontsize=10)

# Create the section plot..
# Resample data if necessary
#if duration * len(waveform) > 6000:
#    waveform.decimate(20, strict_length=False, no_filter=True)
if 'section' in plots:
    fig = plt.figure(figsize=(18, 13))
    plt.title('Section plot for '+eqname+' '+eqtime+" lat:"+str(eqlat)+" lon:"+str(eqlon), fontsize=12, y=1.07)
    if 'distance' in plots:
        waveform.plot(size=(1024,750), type='section', recordlength=duration, time_down=True, linewidth=1.5, grid_linewidth=.5, show=False, fig=fig, color='black', method='full', starttime=starttime, plot_dx=sectiondx, ev_coord = eqlatlon, dist_degree=False, alpha=0.50)
        # Plot customization: Add station labels to offset axis
        ax = fig.axes[0]
        transform = blended_transform_factory(ax.transData, ax.transAxes)
        for tr in waveform:
            ax.text(float(tr.stats.distance)/1e3, 1.0, tr.stats.network + " " + tr.stats.station, rotation=270,
                    va="bottom", ha="center", transform=transform, zorder=10, fontsize=10)
        plt.savefig('../Plots/'+nospaces(eqname)+'-Section-by-distance.png')
    elif 'angle' in plots:
        waveform.plot(size=(1024,750), type='section', recordlength=duration, time_down=True, linewidth=1.5, grid_linewidth=.5, show=False, fig=fig, color='black', method='full', starttime=starttime, plot_dx=angledx, ev_coord = eqlatlon, dist_degree=True, alpha=0.50)

        # Plot customization: Add station labels to offset axis
        ax = fig.axes[0]
        transform = blended_transform_factory(ax.transData, ax.transAxes)
        for tr in waveform:
            ax.text(float(tr.stats.distance), 1.0, tr.stats.network + " " + tr.stats.station, rotation=270,
                    va="bottom", ha="center", transform=transform, zorder=10, fontsize=10)
    #plt.show()
        plt.savefig('../Plots/'+nospaces(eqname)+'-Section-by-angle.png')
