# Open the read_inventory and UTCDateTime modules from the obspy libraries
from obspy import read_inventory, UTCDateTime
# Set up the list of raspberry shake seismometers
seismometers=['RB5E8','RD93E','R82BD','R7FA5','R0353','R9FEE', 'RAD67','R480A','S4BA9','REE67','R614C','RD232','R6C10','RD884','R0B2D','SDF0E','R92B8']
# add the time of the earthquake
eqtime = '2021-08-14 15:00:06'
# Loop through the seismometers in the list
for seismometer in seismometers:
    try:
        # Read the inventory information from fdsnws
        inventory = read_inventory('https://fdsnws.raspberryshakedata.com/fdsnws/station/1/query?network=AM&station=%s&level=resp&format=xml&nodata=404&starttime=%s'
                                   % (seismometer, str(UTCDateTime(eqtime))))
        # set up the name and path of the seismometer .xml file
        invfile = "../Data/" + seismometer + ".xml"
        # Write the file to disk, in the same folder as the downloaded .mseed data file
        inventory.write(invfile,format="STATIONXML")
        print("Wrote station.xml file for " + seismometer)
    except:
        print("Unable to load inventory for " + seismometer)
