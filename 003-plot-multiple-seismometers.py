from obspy.clients.fdsn import Client
from obspy import UTCDateTime, read, Stream

# list of stations
seislist=['RB30C','RB5E8','RD93E','R480A','R82BD','R7FA5'] # Seismometer ID
locations=['FL','PZ','RD','CM','RL','TS'] # Names for the seismometers

# set the data window
#start_1="2019-08-08 16:52:06" 
record_start="2019-08-08 16:52:06"
starttime=UTCDateTime(record_start)
# end time is 8 seconds later
endtime=starttime+8

# read in list of Raspberry Shakes by looping through the list from the second seismometer
client=Client(base_url='https://fdsnws.raspberryshakedata.com/')
# create a new variable of type "obspy stream" containing no data
# see: https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.html#obspy.core.stream.Stream
waveform = Stream()
# loop through the list of seismometers. The index for a list starts at 0 in Python
for x in range (0,len(seislist)):
    # add to the waveform with additional channels. If the station has no data, this will skip to the except option
    try:
        waveform+=client.get_waveforms('AM',seislist[x],'00','EHZ',starttime,endtime)
    except:
        print(seislist[x], "from", locations[x], "not returned from server, skipping this trace.")
    
# place the mean value at zero on the axes
waveform.detrend(type='demean')

# un comment the following line to filter the data, if required
#waveform.filter('bandpass', freqmin=0.5, freqmax=10, corners=4)

# automerge=False stops all of the traces with the same ID from being merged (and printed in a different order)
# method='full' allows zooming of the matplotlib graph when there are lots of data points, but slows down the performance
# equal_scale=False means that each subplot has a different y-axis scale.
# Change this to equal_scale=True if you want to see how amplitude decreases with distance.  
# print to screen, the code continues after you close this plot
waveform.plot(size=(1024,800),type='normal', automerge=False, equal_scale=False, starttime=starttime, endtime=endtime)
#print to file
waveform.plot(size=(1024,800),type='normal', automerge=False, equal_scale=False, starttime=starttime, endtime=endtime, outfile='Multiple-seismometers.png')
# end of the program plotting data from multiple seismometers to the same plot as separate traces. 