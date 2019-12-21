# simple script to plot on day's data from a Raspberry Shake
# import functions from the obspy library - you may need to add the obspy library to your ide (Thonny) before you start
# the code does give a warning about deprecated features in numpy when it runs - these can be ignored
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
# define the date-time for the start of the record
record_start="2019-12-20 00:00:00"
# convert this into UTC Date Time format
starttime = UTCDateTime(record_start)
# The end time is one day later, there are 86400 seconds in a day
endtime = starttime + 86400
# define which client will be the source of the data, this is the client for the Raspberry Shakes
client = Client(base_url='https://fdsnws.raspberryshakedata.com/')
# find the name of the seismometer from stationview at https://raspberryshake.net/stationview/
# for a station with a short period seismometer EHZ. Some Shakes have an SHZ sensor instead.
st = client.get_waveforms('AM', 'R7FA5', '00', 'EHZ', starttime, endtime)
#st = client.get_waveforms('AM', 'R7FA5', '00', 'SHZ', starttime, endtime)
# place the mean value at zero on the axes
st.detrend(type='demean')
# Apply a bandpass filter to cut out low-frequency and high-frequency noise
st.filter("bandpass", freqmin=0.7, freqmax = 2, corners=4)
# Create a dayplot in a file for the seismometer with worldwide events of M5 or more labelled. This will overwrite any previous versions of Daily-plot.png.
# The horizontal axis contains 60 minutes of data
st.plot(type="dayplot", interval=60, right_vertical_labels=False, one_tick_per_line = True, color = ['k', 'r', 'b', 'g'], show_y_UTC_label=False, vertical_scaling_range=400, events={'min_magnitude': 5}, outfile="Daily-plot.png", size=(1920,1080))
# Print the plot to screen
st.plot(type="dayplot", interval=60, right_vertical_labels=False, one_tick_per_line = True, color = ['k', 'r', 'b', 'g'], show_y_UTC_label=False, vertical_scaling_range=400, events={'min_magnitude': 5}, size=(1200,800))
# The End - this was a simple program for printing dayplots to screen
