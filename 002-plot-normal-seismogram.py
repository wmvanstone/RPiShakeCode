#!/usr/bin/env python3
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
import matplotlib.pyplot as plt

record_start="2019-08-08 16:52:05"
# convert this into UTC Date Time format
starttime = UTCDateTime(record_start)
# The end time is 8 seconds later
endtime = starttime + 8
# define which client will be the source of the data, this is the client for the Raspberry Shakes
client = Client(base_url='https://fdsnws.raspberryshakedata.com/')
# find the name of the seismometer from stationview at https://raspberryshake.net/stationview/
# for a station with a short period seismometer EHZ. Some Shakes have an SHZ sensor instead.
waveform = client.get_waveforms('AM', 'RB5E8', '00', 'EHZ', starttime, endtime)
#waveform = client.get_waveforms('AM', 'RB5E8', '00', 'SHZ', starttime, endtime)
# place the mean value at zero on the axes
waveform.detrend(type='demean')
# To apply a bandpass filter to cut out low-frequency and high-frequency noise uncomment the next line
# In the example case, of the Helston earthquake on 08/08/2019 with seismometer RB5E8, the data is fine unfiltered
#waveform.filter("bandpass", freqmin=0.7, freqmax = 20, corners=4)
# Create a normal plot in a file for the seismometer. This will overwrite any previous versions of Normal-plot-Penzance.png.
filename = "Normal-plot-Penzance.png"
# print to file
waveform.plot(outfile=filename,size=(1920,350),type='normal')
# print to screen
waveform.plot(size=(1200,400),type='normal')
# End of the program for plotting a short section of one seismometer

