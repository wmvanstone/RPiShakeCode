from obspy.clients.fdsn import Client
from obspy import UTCDateTime, Stream
import numpy as np
import matplotlib.pyplot as plt
from obspy.taup import TauPyModel
from matplotlib.cm import get_cmap
from obspy.geodetics.base import locations2degrees
from obspy.geodetics import gps2dist_azimuth
import matplotlib.transforms as transforms
from matplotlib.ticker import FormatStrFormatter
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

def nospaces(string):
    out = ""
    for l in string.upper():
        if l in "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789":
            out += l
        else:
            out += "_"
    return out

model = TauPyModel(model='iasp91')
#PHASES = sorted(["P", "pP", "PP"]) # list of phases for which to compute theoretical times
#PHASES = sorted(["PKP", "PKIKP", "PKiKP", "pPKiKP", "SKiKP", "SKP"]) # list of phases for which to compute theoretical times
PHASES = sorted(["P", "pP", "PP", "S", "Pdiff", "PKP", "PKIKP", "PcP", "ScP", "ScS", "PKiKP", "SKiKP", "SKP", "SKS"]) # list of phases for which to compute theoretical times
#PHASES = sorted(["P", "pP"])
#PHASES = sorted(["PKiKP", "PKIKP"])
#PHASES = sorted(["PKP", "PKIKP", "PKiKP", "pPKIKP", "pPKiKP"]) # list of phases for which to compute theoretical times
#PHASES = sorted(["PKP", "PKIKP", "PKiKP"]) # list of phases for which to compute theoretical times
#PHASES = sorted(["P", "pP", "PP", "PPP"])

DATA_PROVIDER = "RASPISHAKE"

# Event details
URL = 'https://earthquake.usgs.gov/earthquakes/eventpage/us6000ah9t/executive'
EQNAME = 'M7.4 - 12km SSW of Santa María Zapotitlán'
EQLAT = 16.029
EQLON = -95.901
EQZ = 26.3
EQTIME = '2020-06-23 15:29:05'
FILE_STEM = 'Mexico-2020-06-23'
MAGNITUDE = "M7.4"
MAPFILE='Mexico-2020-06-23-map.png'
LOGOS='logos.png'
LABEL=nospaces(EQNAME)
F1=0.03
F2=1.0

# set the data window
STARTTIME=UTCDateTime(EQTIME)
DURATION = 5000
PSTART = 0
PEND = 5000
SAMESCALE = True
YMAX = 7e-5

# Home station
NETWORK = 'AM'   # AM = RaspberryShake network
STATION = "R7FA5"  # Station code of local station to plot
STA_LAT = 50.2609  # Latitude of local station  
STA_LON = -5.0434  # Longitude of local station
CHANNEL = 'EHZ'  # channel to grab data for (e.g. EHZ, SHZ, EHE, EHN)
LOCATION = "TRURO SCHOOL"
DISTANCE=locations2degrees(EQLAT, EQLON, STA_LAT, STA_LON) # Station dist in degrees from epicentre
STA_DIST, _, _ = gps2dist_azimuth(STA_LAT, STA_LON, EQLAT, EQLON)   # Station dist in m from epicentre

# list of stations
SEISLIST=['RB30C','RB5E8','RD93E','R82BD','R7FA5','R0353','R9FEE', 'R303A']
LOCATIONS=['Falmouth','Penzance','Redruth','Richard Lander','Truro School','Penair','Truro High', 'Constantine']
LATITUDES=[50.1486,50.1179833,50.2344,50.2596,50.2609,50.2673,50.2570,50.117]
LONGITUDES=[-5.0945,-5.5391226,-5.2384,-5.1027,-5.0434,-5.0299,-5.0566,-5.175]

# list of stations
#SEISLIST=['RB30C','RB5E8','RD93E','R7FA5','R0353','R9FEE', 'R303A']
#LOCATIONS=['Falmouth','Penzance','Redruth','Truro School','Penair','Truro High', 'Constantine']
#LATITUDES=[50.1486,50.1179833,50.2344,50.2609,50.2673,50.2570,50.117]
#LONGITUDES=[-5.0945,-5.5391226,-5.2384,-5.0434,-5.0299,-5.0566,-5.175]

# fill the list of seismometers and sort it by distance from the epicentre
seismometers = [] # empty list of seismometers
for stationID in range(len(SEISLIST)):
    distance = locations2degrees(EQLAT, EQLON, LATITUDES[stationID], LONGITUDES[stationID])
    seismometers.append([SEISLIST[stationID], round(LATITUDES[stationID],4), round(LONGITUDES[stationID],4), round(distance,4), LOCATIONS[stationID]])
seismometers.sort(key = lambda i: i[3]) # Contains 'Code', Latitude, Longitude, distance (degrees), location name: sort by distance from epicentre

# Pretty paired colors. Reorder to have saturated colors first and remove
# some colors at the end. This cmap is compatible with obspy taup
CMAP = get_cmap('Paired', lut=12)
COLORS = ['#%02x%02x%02x' % tuple(int(col * 255) for col in CMAP(i)[:3]) for i in range(12)]
COLORS = COLORS[1:][::2][:-1] + COLORS[::2][:-1]

# read in list of Raspberry Shakes by looping through the list from the second seismometer
client=Client(DATA_PROVIDER)
# create a new variable of type "obspy stream" containing no data
waveform = Stream()
loaded = []
# loop through the list of seismometers. The index for a list starts at 0 in Python
for x in range (0,len(seismometers)):
    # add to the waveform with additional channels. If the station has no data, this will skip to the except option
    try:
        st=client.get_waveforms('AM',seismometers[x][0],'00','EHZ',STARTTIME,STARTTIME+DURATION,attach_response=True)
        st.merge(method=0, fill_value='latest')
        st.detrend(type='demean')
        pre_filt = [0.01, 0.02, 12.0, 15]
        st.remove_response(pre_filt=pre_filt,output="VEL",water_level=40,taper=True)#, plot=True)
#        st[0].trim(starttime=STARTTIME+PSTART, endtime=STARTTIME+PEND, pad=False, nearest_sample=True, fill_value=None)
        waveform += st
        loaded.append(seismometers[x]) # Add the details of loaded seismometers to the loaded list
    except:
        print(seismometers[x][0], "from", seismometers[x][4], "not returned from server, skipping this trace.")

# Filter the data
waveform.filter("bandpass", freqmin=F1, freqmax = F2, corners=4)
FILTERLABEL = '"bandpass", freqmin=' + str(F1) + ', freqmax=' + str(F2) + ', corners=4'

# Find the maximum y value to scale the axes
ylim = 0
for trace in waveform:        # find maximum value from traces for plotting when SAMESCALE is True
    print(trace.max())
    if abs(trace.max()) * 1.1 > ylim:
        ylim = abs(trace.max()) * 1.1
ylim = YMAX

# Plot figure with subplots of different sizes
fig = plt.figure(1, figsize=(19, 11))
# set up the page margins, height and spacing of the plots for traces
margin = 0.06
twidth = 0.60
spacing = 0.00
theight = ((1.0 - 2 * margin) / waveform.count()) - spacing
axes = [] # A list for the trace axis graphs
# set up an axis covering the whole plot area for annotations and turn off the axes and ticks
annotations = fig.add_axes([0, 0, 1, 1])
annotations.spines['top'].set_visible(False)
annotations.spines['right'].set_visible(False)
annotations.spines['bottom'].set_visible(False)
annotations.spines['left'].set_visible(False)
annotations.get_xaxis().set_ticks([])
annotations.get_xaxis().set_ticks([])
# Add descriptive labels to the plot
networksummary = "from the Raspberry Shake Network."
titletxt = EQNAME + " " + str(EQTIME[0:10]) + "\n" + networksummary + "\nSee: " + URL + " for more info."
annotations.text(0.36, 0.99, titletxt, transform=annotations.transAxes, horizontalalignment='center', verticalalignment='top')
#annotations.text(0.01, 0.50, 'Amplitude [m]', rotation=90, horizontalalignment='center', verticalalignment='center')
annotations.text(0.01, 0.50, 'Velocity [m/s]', rotation=90, horizontalalignment='center', verticalalignment='center')
annotations.text(0.36, 0.01, 'Time since earthquake [s]', horizontalalignment='center', verticalalignment='bottom')
annotations.text(0.01, 0.01, "Channel: " + CHANNEL + ", filter: " + FILTERLABEL, horizontalalignment='left', verticalalignment='bottom')
annotations.text(0.67, 0.51, "Epicentre: Lat: " + str(round(EQLAT, 2)) + "°  Lon: " + str(round(EQLON,2)) + "°  Depth: " + str(round(EQZ, 1)) + "km  Time: " + EQTIME + "UTC", horizontalalignment='left', verticalalignment='center')
annotations.text(0.67, 0.16, "Station: " + STATION + "-" + LOCATION + "  Lat: " + str(round(STA_LAT, 2)) + "°  Lon: " + str(round(STA_LON,2)) + "°  Dist: " + str(STA_DIST//1000) + "km", horizontalalignment='left', verticalalignment='center')
annotations.text(0.67, 0.56, "Ray paths and arrivals\ncalculated from iasp91 model.", horizontalalignment='left', verticalalignment='center') 
# Add the logos and map in the right-hand panel
# logo
arr_logos = plt.imread(LOGOS)
imagebox = OffsetImage(arr_logos, zoom=0.5)
ab = AnnotationBbox(imagebox, (0.830, 0.08), frameon=False)
annotations.add_artist(ab)
# map - generated from folium
arr_map = plt.imread(MAPFILE)
mapbox = OffsetImage(arr_map, zoom=0.24)
mb = AnnotationBbox(mapbox, (0.83, 0.335), frameon=False)
annotations.add_artist(mb)
# set up stream plots for each seismometer
for i in range(waveform.count()):
    if i == 0: # bottom axis first, including tick labels
        axes.append(fig.add_axes([margin, margin + ((theight + spacing) * i), twidth, theight]))
    else: # tick labels switched off for subsequent axes
        axes.append(fig.add_axes([margin, margin + ((theight + spacing) * i), twidth, theight], sharex=axes[0]))
        plt.setp(axes[i].get_xticklabels(), visible=False)
    # Create a list of elapsed times for the stream file (what does this do in the case of data gaps?
    time = np.arange(0, waveform[i].count()/100, 0.01)
    # set up plot parameters
    axes[i].set_xlim([PSTART,PEND])
    if not SAMESCALE:
        ylim = abs(waveform[i].max()) * 1.1
    axes[i].set_ylim([-1*ylim,ylim])
    axes[i].yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
    axes[i].tick_params(axis="both", direction="in", which="both", right=True, top=True)
    # plot the data
    axes[i].plot(time, waveform[i], linewidth=0.75)
    textstr = loaded[i][0] + " " + str(round(loaded[i][3], 2)) + "° " + loaded[i][4]
    plt.text(0.5, 0.96, textstr, transform=axes[i].transAxes, horizontalalignment='center', verticalalignment='top')
    # Add the arrival times and vertical lines
    plotted_arrivals = []
    for j, phase in enumerate(PHASES):
        color = COLORS[PHASES.index(phase) % len(COLORS)]
        arrivals = model.get_travel_times(source_depth_in_km=EQZ, distance_in_degree=loaded[i][3], phase_list=[phase])
        printed = False
        if arrivals:
            for k in range(len(arrivals)):
                instring = str(arrivals[k])
                phaseline = instring.split(" ")
                phasetime = float(phaseline[4])
                if phaseline[0] == phase and printed == False and (PSTART <  phasetime < PEND):
                    plotted_arrivals.append(tuple([round(float(phaseline[4]), 2), phaseline[0], color]))
                    printed = True
    if plotted_arrivals:
        plotted_arrivals.sort(key=lambda tup: tup[0])   #sorts the arrivals to be plotted by arrival time
        trans = transforms.blended_transform_factory(axes[i].transData, axes[i].transAxes)
        phase_ypos = 0

        for phase_time, phase_name, color_plot in plotted_arrivals:
            axes[i].vlines(phase_time, ymin=0, ymax=1, alpha=.50, color=color_plot, ls='--', zorder=1, transform=trans)
            plt.text(phase_time, 0, phase_name+" ", alpha=.50, c=color_plot, fontsize=11, horizontalalignment='right', verticalalignment='bottom', zorder=1, transform=trans)
# Add the inset picture of the globe at x, y, width, height, as a fraction of the parent plot
ax1 = fig.add_axes([0.64, 0.545, 0.38, 0.38], polar=True)
# Plot all pre-determined phases
for phase in PHASES:
    arrivals = model.get_ray_paths(EQZ, DISTANCE, phase_list=[phase])
    try:
        ax1 = arrivals.plot_rays(plot_type='spherical', legend=True, label_arrivals=False, plot_all=True, phase_list = PHASES, show=False, ax=ax1)
    except:
        print(phase + " not present.")
    
# Annotate regions of the globe
ax1.text(0, 0, 'Solid\ninner\ncore', horizontalalignment='center', alpha=0.5, verticalalignment='center', bbox=dict(facecolor='white', edgecolor='none', alpha=0.5))
ocr = (model.model.radius_of_planet - (model.model.s_mod.v_mod.iocb_depth + model.model.s_mod.v_mod.cmb_depth) / 2)
ax1.text(np.deg2rad(180), ocr, 'Fluid outer core', alpha=0.5, horizontalalignment='center', bbox=dict(facecolor='white', edgecolor='none', alpha=0.5))
mr = model.model.radius_of_planet - model.model.s_mod.v_mod.cmb_depth / 2
ax1.text(np.deg2rad(180), mr, 'Solid mantle', alpha=0.5, horizontalalignment='center', bbox=dict(facecolor='white', edgecolor='none', alpha=0.5))
rad = model.model.radius_of_planet*1.15
ax1.text(np.deg2rad(DISTANCE), rad, " " + STATION, horizontalalignment='left', bbox=dict(facecolor='white', edgecolor='none', alpha=0))
ax1.text(np.deg2rad(0), rad, EQNAME + "\nEpicentral distance: " + str(round(DISTANCE, 1)) + "°" + "\nEQ Depth: " + str(EQZ) + "km",
        horizontalalignment='left', bbox=dict(facecolor='white', edgecolor='none', alpha=0))

# Save the plot to file, then print to screen
plt.savefig(nospaces(LABEL)+'-Summary.png')
#plt.tight_layout()
plt.show()

print(EQNAME + " at " + EQTIME + " UTC recorded on the Cornwall Schools @uniteddownsgeo @raspishake network in Cornwall, UK, processed by @TruroSchool."
      + " For more info see: " + URL + ". Thanks to @obspy and @matplotlib for software libraries.")