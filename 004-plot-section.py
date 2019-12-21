# Program to plot a section using seismometer locations, epicentre lat/long and rupture time
# import obspy libraries
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy import read, Stream
# import geodetics library for lat/long manipulation
from obspy.geodetics import gps2dist_azimuth
# import matplotlib libraries for extra plotting capabilities
import matplotlib.pyplot as plt
from matplotlib.transforms import blended_transform_factory

# station statistics in four separate lists. CCA1 is a BGS seismometers, the others are Raspberry Shakes
seislist=['CCA1','RB30C','RB5E8','RD93E','R82BD','R7FA5']
locations=['CM','FL','PZ','RD','RL','TS']
# Approximate latitudes and longitudes obtained from BGS web site and Raspberry Shake instrument response files on the shake map.
# https://earthquakes.bgs.ac.uk/monitoring/broadband_stationbook.html for the BGS stations. Click on the station ID to get the info page. 
# https://raspberryshake.org/community/station-view/ for the Raspberry Shakes.
# Double click on the triangle station marker, then download instrument response and look there for the latitude and longitude.
latitudes=[50.1867,50.1486,50.1179833,50.2344,50.2596,50.2609]
longitudes=[-5.2273,-5.0945,-5.5391226,-5.2384,-5.1027,-5.0434]

# Earthquakes' epicentre - by trial and error. This is not the same as the BGS published location.
# I used Google Earth Pro to get the lat and long, adjusting my position until the traces liend up on the section plot.
eq_lat = 50.054
eq_lon = -5.276

# set the data window from the likely rupture time
start_1="2019-08-08 16:52:05.09"
starttime=UTCDateTime(start_1)
endtime=starttime+8

# start reading the data
# create a new variable of type "obspy stream" containing no data
# see: https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.html#obspy.core.stream.Stream
waveform=Stream()
# read in BGS seismometer at Carnmenellis from IRIS. Not all BGS seismometers transmit data to IRIS. 
client=Client("IRIS")
waveform+=client.get_waveforms('GB','CCA1','--','HHZ',starttime,endtime)

# read in list of Raspberry Shakes by looping through the list from the second seismometer
client=Client(base_url='https://fdsnws.raspberryshakedata.com/')
# remember, Python lists are indexed from zero, so this looks from the second item in the list
for x in range (1,len(seislist)):
    waveform+=client.get_waveforms('AM',seislist[x],'00','EHZ',starttime,endtime)
    
# place the mean value at zero on the axes
waveform.detrend(type='demean')

# work out the distance of each seismometer from the epicentre and add the values to the stats for each wave form
for y in range(len(seislist)):
    waveform[y].stats["coordinates"] = {} # add the coordinates to your dictionary, needed for the section plot
    waveform[y].stats["coordinates"]["latitude"] = latitudes[y]
    waveform[y].stats["coordinates"]["longitude"] = longitudes[y]
    waveform[y].stats["distance"] = gps2dist_azimuth(waveform[y].stats.coordinates.latitude, waveform[y].stats.coordinates.longitude, eq_lat, eq_lon)[0]
    # Set the abbreviated name of the location in the network field for the plot title
    waveform[y].stats.network = locations[y]

# automerge=False stops all of the traces with the same ID from being merged (and printed in a different order)
# method='full' allows zooming of the matplotlib graph when there are lots of data points, but slows down the performance
# Print to file
waveform.plot(size=(1024,800),type='normal', automerge=False, equal_scale=False, starttime=starttime, endtime=endtime, outfile='CornishQuake.png')
# Print to screen
waveform.plot(size=(1024,800),type='normal', automerge=False, equal_scale=False, starttime=starttime, endtime=endtime)

# filter the data, if necessary
#waveform.filter('bandpass', freqmin=0.1, freqmax=10)

# Create the section plot using matplotlib. This is slightly different to the previous plot.
# Define the size of the canvas
fig = plt.figure(figsize=(14, 14))
# Add a title with a vertical offset of 0.08 so the main title does not obscure the station labels
plt.title('Section plot for epicentre at the Loe and rupture at UTC 16:52:05.09', fontsize=12, y=1.08)
# Create the plot 600pxx600px, 8sec record, 5km x tick marks, 
waveform.plot(size=(600,600), type='section', recordlength=8, plot_dx=5e3, time_down=True, linewidth=1.5, grid_linewidth=.5, show=False, fig=fig, starttime=starttime, endtime=endtime, color='station')

# Plot customization: Add station labels to offset axis
ax = fig.axes[0]
transform = blended_transform_factory(ax.transData, ax.transAxes)
for tr in waveform:
    ax.text(tr.stats.distance / 1e3, 1.0, tr.stats.network + " " + tr.stats.station, rotation=270,
            va="bottom", ha="center", transform=transform, zorder=10, fontsize=10)
# print to file
plt.savefig('CornishSection.png')
# print to screen
plt.show()
# end of the section program for the Helston earthquake of 2019.
# The method I have used to find the epicentre and rupture time takes no account of the depth of the earthquake focus.