# importing the required libraries for requesting the data
import obspy
from obspy.clients.fdsn import Client
from obspy import UTCDateTime, Stream
import os
import wmv_utils
import requests
from csv import DictReader
from obspy.geodetics.base import locations2degrees
from obspy.geodetics import gps2dist_azimuth
from obspy.taup import TauPyModel
from matplotlib.cm import get_cmap
import matplotlib.transforms as transforms
import numpy as np
import matplotlib.pyplot as plt

model = TauPyModel(model='iasp91')
client=Client("RASPISHAKE")

# URL for USGS Earthquake data
DATA_URL = 'https://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/4.5_day.csv'
#DATA_URL = 'http://192.168.1.171/seismic/4.5_day.csv'
DURATION = 400
# St Day details
seismometer = 'RAD67'
STA_LAT = 50.2385
STA_LON = -5.1822
client=Client("RASPISHAKE")
PHASES = sorted(["P", "pP", "pPP", "PP", "PPP", "S", "Pdiff", "PKP", "PKIKP", "PcP", "ScP", "ScS", "PKiKP", "SKiKP", "SKP", "SKS"])
CMAP = get_cmap('Paired', lut=12)
COLORS = ['#%02x%02x%02x' % tuple(int(col * 255) for col in CMAP(i)[:3]) for i in range(12)]
COLORS = COLORS[1:][::2][:-1] + COLORS[::2][:-1]

DATA_PROVIDER = "RASPISHAKE"
F1=0.9
F2=6.0
YAXIS = "VEL"

print("Downloading", DATA_URL)
resp = requests.get(DATA_URL)
quakes = list(DictReader(resp.text.splitlines()))
detected = []

for i, quake in enumerate(quakes):
    if float(quake['depth']) < 0:
        eqz = "0"
        quake['depth'] = "0"
    else:
        eqz = quake['depth']
    eqmag = quake['mag']
    eqloc = quake['place']
    eqtime = quake['time']
    url = "'https://earthquake.usgs.gov/earthquakes/eventpage/"+ quake['id'] + "/executive'"
    eqname = "'M" + quake['mag'] + " - " + quake['place'] + "'"
    eqlat = float(quake['latitude'])
    eqlon = float(quake['longitude'])
    eqz = float(quake['depth'])
    #eqtime = "'" + quake['time'][0:10] + " " + quake['time'][11:23] + "'"
    file_stem = "'" + quake['place'].split(", ")[-1] + "-" + quake['time'][0:10] + "'"
    distance=locations2degrees(float(quake['latitude']), float(quake['longitude']), STA_LAT, STA_LON) # Station dist in degrees from epicentre
    sta_dist, _, _ = gps2dist_azimuth(STA_LAT, STA_LON, float(quake['latitude']), float(quake['longitude']))   # Station dist in m from epicentre
    print(quake['mag'] + "\t" + quake['place'] + "\t" + str(quake['latitude']) + "\t" + str(quake['longitude']) + "\t" + str(quake['depth']) + "\t" + quake['time'][0:10] + " " + quake['time'][11:23] + "\t" + str(sta_dist/1000) + "\t" + str(distance) + "\t" + "https://earthquake.usgs.gov/earthquakes/eventpage/"+ quake['id'] + "/executive")
    arrivals=model.get_travel_times(source_depth_in_km=float(quake['depth']), distance_in_degree=distance)
    first_arrival = float(str(arrivals[0]).split(" ")[4])
    # Read the seismic stream
    try:
        eqstart = UTCDateTime(quake['time']) + first_arrival - 60
        waveform = client.get_waveforms('AM', seismometer, '00', 'EHZ', eqstart, eqstart+DURATION, attach_response=True)
        waveform.merge(method=0, fill_value='latest')
        pre_filt = [0.01, 0.1, 12.0, 15]
        waveform.remove_response(pre_filt=pre_filt,output=YAXIS,water_level=40,taper=True)#, plot=True)
        waveform.filter("bandpass", freqmin=F1, freqmax = F2, corners=4)
        time = np.arange(0, waveform[0].count()/100, 0.01)
        fig, ax = plt.subplots(1, 1)
        ax.plot(time, waveform[0], linewidth=0.75)
        plotted_arrivals = []
        for j, phase in enumerate(PHASES):
            color = COLORS[PHASES.index(phase) % len(COLORS)]
            #arrivals = model.get_travel_times(source_depth_in_km=eqz, distance_in_degree=distance, phase_list=[phase])
            printed = False
            if arrivals:
                for k in range(len(arrivals)):
                    instring = str(arrivals[k])
                    phaseline = instring.split(" ")
                    phasetime = float(phaseline[4])
                    if phaseline[0] == phase and printed == False and ((phasetime - first_arrival + 60) < DURATION):
                        plotted_arrivals.append(tuple([round(float(phaseline[4]), 2), phaseline[0], color]))
                        printed = True
        if plotted_arrivals:
            plotted_arrivals.sort(key=lambda tup: tup[0])   #sorts the arrivals to be plotted by arrival time
            trans = transforms.blended_transform_factory(ax.transData, ax.transData)
            phase_ypos = 0

            for phase_time, phase_name, color_plot in plotted_arrivals:
                ax.vlines(phase_time - first_arrival + 60, ymin=-1e-7, ymax=1e-7, alpha=.50, color=color_plot, ls='--', zorder=1, transform=trans)
                ax.text(phase_time - first_arrival + 60, -1e-7, phase_name+" ", alpha=.50, c=color_plot, fontsize=11, horizontalalignment='right', verticalalignment='bottom', zorder=1, transform=trans)
        plt.show()
        canyouseeit = input("Can you see a response")
        if canyouseeit.lower() == "y":
            detected.append(1)
        else:
            detected.append(0)
    except:
        print("No data returned")
        detected.append(0)
print("\nDetected Quakes\n")
for i, quake in enumerate(quakes):
    if detected[i] == 1:
        sta_dist, _, _ = gps2dist_azimuth(STA_LAT, STA_LON, float(quake['latitude']), float(quake['longitude']))   # Station dist in m from epicentre
        distance=locations2degrees(float(quake['latitude']), float(quake['longitude']), STA_LAT, STA_LON) # Station dist in degrees from epicentre
        arrivals=model.get_travel_times(source_depth_in_km=float(quake['depth']), distance_in_degree=distance)
        first_arrival = float(str(arrivals[0]).split(" ")[4])
        print(quake['mag'] + "\t" + quake['place'] + "\t" + str(quake['latitude']) + "\t" + str(quake['longitude']) + "\t" + str(quake['depth']) + "\t" + quake['time'][0:10] + " " + quake['time'][11:23] + "\t" + str(sta_dist/1000) + "\t" + str(distance) + "\t" + str(first_arrival) + "\t" + "https://earthquake.usgs.gov/earthquakes/eventpage/"+ quake['id'] + "/executive")