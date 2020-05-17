import folium
from obspy.geodetics.base import locations2degrees
import math
import numpy as np
from numpy import arctan2, sqrt
import numexpr as ne
from sympy import *
from numpy import cross, eye, dot
from scipy.linalg import expm, norm

PI = 3.14159
R = 6367.5

def M(axis, theta):
    return expm(cross(eye(3), axis/norm(axis)*theta))

#v, axis, theta = [3,5,0], [4,4,1], 1.2
#M0 = M(axis, theta)

#print(dot(M0,v))
# [ 2.74911638  4.77180932  1.91629719]

#input("Hit enter to continue")

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

#v = [3, 5, 0]
#axis = [4, 4, 1]
#theta = 1.2 

#print(np.dot(rotation_matrix(axis, theta), v)) 
# [ 2.74911638  4.77180932  1.91629719]

#input("Hit enter to continue")

def appendSpherical_np(xyz):
    ptsnew = np.hstack((xyz, np.zeros(xyz.shape)))
    xy = xyz[:,0]**2 + xyz[:,1]**2
    ptsnew[:,3] = np.sqrt(xy + xyz[:,2]**2)
    ptsnew[:,4] = np.arctan2(np.sqrt(xy), xyz[:,2]) # for elevation angle defined from Z-axis down
    #ptsnew[:,4] = np.arctan2(xyz[:,2], np.sqrt(xy)) # for elevation angle defined from XY-plane up
    ptsnew[:,5] = np.arctan2(xyz[:,1], xyz[:,0])
    return ptsnew

def asCartesian(rthetaphi):
    #takes list rthetaphi (single coord)
    r       = rthetaphi[0]
    theta   = rthetaphi[1]*PI/180 # to radian
    phi     = rthetaphi[2]*PI/180
    x = r * math.sin( theta ) * math.cos( phi )
    y = r * math.sin( theta ) * math.sin( phi )
    z = r * math.cos( theta )
    return [x,y,z]

def asSpherical(xyz):
    #takes list xyz (single coord)
    x       = xyz[0]
    y       = xyz[1]
    z       = xyz[2]
    r       =  sqrt(x*x + y*y + z*z)
    theta   =  math.acos(z/r)*180/ PI #to degrees
    phi     =  math.atan2(y,x)*180/ PI
    return [r,theta,phi]

def cart2sph(x,y,z, ceval=ne.evaluate):
    """ x, y, z :  ndarray coordinates
        ceval: backend to use: 
              - eval :  pure Numpy
              - numexpr.evaluate:  Numexpr """
    azimuth = ceval('arctan2(y,x)')
    xy2 = ceval('x**2 + y**2')
    elevation = ceval('arctan2(z, sqrt(xy2))')
    r = eval('sqrt(xy2 + z**2)')
    return azimuth, elevation, r

#v = [-2.13091326,-0.0058279,0.83697319]
#x = v[0]
#y = v[1]
#z = v[2]

#r = math.sqrt(x*x + y*y + z*z)
#phi = math.atan2(y,x)
#theta = math.acos(z/r)
#input("Hit enter to continue")

#print(r, phi, theta)

#test = asCartesian(asSpherical([-2.13091326,-0.0058279,0.83697319]))

#print(test)
URL = 'https://earthquake.usgs.gov/earthquakes/eventpage/us70009hnm/executive'
EQNAME = 'M5.9 - 60km E of Lakatoro, Vanuatu'
EQLAT = -16.0936
EQLON = 167.9806
EQZ = 165.66
EQTIME = '2020-05-16T03:15:43.425Z'
FILE_STEM = 'Vanuatu-2020-05-16'

MAPFILE=FILE_STEM + '-map.png'
LOGOS='logos.png'

# Things to change once for your station
NETWORK = 'AM'   # AM = RaspberryShake network
STATION = "R7FA5"  # Station code of local station to plot
STA_LAT = 50.2609  # Latitude of local station  
STA_LON = -5.0434  # Longitude of local station
CHANNEL = 'EHZ'  # channel to grab data for (e.g. EHZ, SHZ, EHE, EHN)
LOCATION = "Truro School"
# Plot the epicentre and seismometer on the map
map = folium.Map(location=[EQLAT, EQLON],zoom_start=2,tiles='Stamen Terrain')
# The following arguments will centre the map on Caustic.
map = folium.Map(location=[47.8288, 1.9324],zoom_start=6,tiles='Stamen Terrain')

# theta is 90 - latitude
# phi is longitude but phi 180-360 becomes negative longitude
# r = radius of earth (6,378 EQ + 6,357 POLE) / 2 = 6367.5
dend = []
dstart = []

for phi in range(0, 360, 1):
    theta = 140
    dend.append(asCartesian([R, theta, phi]))
    omega = 96.7
    dstart.append(asCartesian([R, omega, phi]))

dendt = []
axis = [math.cos((90+EQLON)/180*PI), math.sin((90+EQLON)/180*PI), 0]
theta = (90-EQLAT)/180*PI
for i in range(len(dend)):
    dend[i]=np.dot(rotation_matrix(axis, theta), dend[i])
    send = asSpherical(dend[i])
    lat = 90 - (send[1])
    long = send[2]
    dendt.append([lat, long])

dstartt = []
axis = [math.cos((90+EQLON)/180*PI), math.sin((90+EQLON)/180*PI), 0]
theta = (90-EQLAT)/180*PI
for i in range(len(dstart)):
    dstart[i]=np.dot(rotation_matrix(axis, theta), dstart[i])
    sstart = asSpherical(dstart[i])
    lat = 90 - (sstart[1])
    long = sstart[2]
    dstartt.append([lat, long])

linepoints=[]
for lon in range(-180, 180, 10):
    # Northern hemisphere
    DISTANCE=locations2degrees(EQLAT, EQLON, 0, lon) # Station dist in degrees from epicentre
    bestfit = abs(DISTANCE-140)
    bestloc = [0, lon]
    plotit = False
    for lat in range(0, 91, 1):
        DISTANCE=locations2degrees(EQLAT, EQLON, lat, lon) # Station dist in degrees from epicentre
        if DISTANCE > 139.6 and DISTANCE < 140.4:
            plotit = True
            for latf in range((lat-1)*300, (lat+2)*300, 1):
                DISTANCE=locations2degrees(EQLAT, EQLON, latf/300, lon)
                fit = abs(DISTANCE-140)
                if fit < bestfit:
                    bestfit = fit
                    bestloc = [latf/300, lon]
    if plotit == True:
        linepoints.append(bestloc)
    # southern hemisphere
        DISTANCE=locations2degrees(EQLAT, EQLON, 0, lon) # Station dist in degrees from epicentre
    bestfit = abs(DISTANCE-140)
    bestloc = [0, lon]
    plotit = False
    for lat in range(-90, 1, 1):
        DISTANCE=locations2degrees(EQLAT, EQLON, lat, lon) # Station dist in degrees from epicentre
        if DISTANCE > 139.6 and DISTANCE < 140.4:
            plotit = True
            for latf in range((lat-1)*300, (lat+2)*300, 1):
                DISTANCE=locations2degrees(EQLAT, EQLON, latf/300, lon)
                fit = abs(DISTANCE-140)
                if fit < bestfit:
                    bestfit = fit
                    bestloc = [latf/300, lon]
    if plotit == True:
        linepoints.append(bestloc)

#folium.vector_layers.PolyLine(linepoints, color="green", weight=2, opacity=0.8, smooth_factor=0.5).add_to(map)
line = []
counter = 0
for i, dp in enumerate(dendt):
    if i > 0:
        if abs(dendt[i-1][0] - dendt[i][0]) > 10 or abs(dendt[i-1][1] - dendt[i][1]) > 10:
            folium.vector_layers.PolyLine(line, color="red", weight=3, opacity=1.0, smooth_factor=0.5).add_to(map)
            line = []
    if dp[0] < -90:
        dp[0] += 180
    if dp[0] > 90:
        dp[0] -= 180
    if dp[1] < -180:
        dp[1] += 360
    if dp[1] > 180:
        dp[1] -= 360
    line.append(dp)
folium.vector_layers.PolyLine(line, color="red", weight=3, opacity=1.0, smooth_factor=0.5).add_to(map)

line = []
counter = 0
for i, dp in enumerate(dstartt):
    if i > 0:
        if abs(dstartt[i-1][0] - dstartt[i][0]) > 10 or abs(dstartt[i-1][1] - dstartt[i][1]) > 10:
            folium.vector_layers.PolyLine(line, color="red", weight=3, opacity=1.0, smooth_factor=0.5).add_to(map)
            line = []
    if dp[0] < -90:
        dp[0] += 180
    if dp[0] > 90:
        dp[0] -= 180
    if dp[1] < -180:
        dp[1] += 360
    if dp[1] > 180:
        dp[1] -= 360
    line.append(dp)
folium.vector_layers.PolyLine(line, color="red", weight=3, opacity=1.0, smooth_factor=0.5).add_to(map)

'''
STATIONS = [
[39.8739,   25.2737,    139.0268,   8.08521143225e-08,  'R9AC0',  15464017.12610259,  63.8190132031462,   314.36286748274676],
[48.2162,   16.3481,    139.4912,   1.24012198719e-07,  'R21E0',  15512026.232576072, 45.386903080507025, 330.48708270854775],
[51.3423,   9.0719, 140.0754,   2.86703703369e-08,  'RC4FB',  15575527.635725247, 33.23599842700995,  339.1674835760409],
[40.7297,   22.9923,    140.1738,   4.13058393523e-08,  'RC574',  15591170.019331088, 60.498471853957795, 316.790179065272],
[51.4234,   6.8339, 140.7428,   6.0645088915e-08,   'RDE9F',  15649547.206057,    30.011949426045884, 341.0963191101199],
[46.973,    15.3948,    140.8208,   2.81219932953e-08,  'R830F',  15660220.298032654, 45.19373918516495,  329.8180786207821],
[37.973,    23.7496,    140.9189,   5.08230866469e-08,  'R7AD4',  15675231.197069323, 64.27509328069624,  312.51392119632123],
[46.3333,   15.9599,    140.9848,   1.8184929551e-08,   'RCF63',  15678722.622802122, 46.49579256213755,  328.66368400640465],
[50.955,    7.1663, 141.0401,   2.77531725314e-08,  'RB99F',  15682738.027138494, 30.777767617374888, 340.43823999631377],
[50.2252,   8.5907, 141.1751,   9.32369137867e-08,  'RA52C',  15698050.774803933, 33.293873827404056, 338.6062523520899],
[49.964,    8.8242, 141.3091,   7.73763462492e-08,  'R74D9',  15713030.505896715, 33.804801530906445, 338.1777091962111],
[48.1441,   12.287, 141.3893,   3.92489175662e-08,  'R06C4',  15722782.58651194,  40.02431899045048,  333.5340139097541],
[50.1081,   7.8955, 141.5161,   2.26880020602e-08,  'R0333',  15735927.126348088, 32.37494113266299,  339.1029048180517],
[39.5405,   21.7645,    141.5738,   8.58100145689e-08,  'RAC91',  15747320.91143851,  60.405997753038875, 315.8870849930575],
[50.9459,   5.3051, 141.6273,   1.07104225243e-07,  'R78B8',  15747881.05902008,  28.074279350863975, 342.06159555428354],
[48.1802,   11.566, 141.6691,   8.23380529745e-08,  'R4A43',  15753809.657778796, 38.99492989160851,  334.16664813071804],
[47.4775,   12.5303,    141.7881,   4.08355153448e-08,  'R6B73',  15767348.634716062, 40.900511430314495, 332.6404427704048],
[49.3243,   8.6952, 141.8854,   4.72325183687e-08,  'R7266',  15777262.278000155, 34.065110096124954, 337.7168934997417],
[52.4595,   -1.9369,    141.9411,   8.07108824554e-07,  'R8118',  15781960.760614118, 16.451417886967825, 349.67343004257185],
[39.6577,   20.8424,    142.1317,   4.0037490716e-08,   'R1388',  15809244.968678603, 59.23388036395095,  316.63182100347433],
[49.8739,   6.2069, 142.2818,   9.75732727564e-08,  'RCEAB',  15820966.248155689, 30.06173508736019,  340.40819745854367],
[51.3243,   1.413,  142.326,    9.34488989481e-07,  'R8C09',  15825170.411010163, 22.042485089244035, 345.9001941804129],
[51.6486,   -0.2178,    142.3846,   1.19991870698e-07,  'R0AEF',  15831509.249612723, 19.402947017522838, 347.6360727565328],
[51.6036,   -0.3337,    142.451,    1.32745781091e-07,  'R0887',  15838889.695852866, 19.24399094062836,  347.7227018187861],
[51.4595,   -0.0145,    142.52, 9.20167598181e-08,  'R2E51',  15846599.297615675, 19.796556176448803, 347.3380303820712],
[51.4234,   -0.0722,    142.5661,   4.41461436921e-07,  'R63BA',  15851735.500129828, 19.72348510103283,  347.3735289040536],
[51.3694,   -0.0722,    142.6169,   7.52574441491e-08,  'RCC45',  15857390.435926516, 19.746979319482367, 347.34371049929035],
[45.9369,   13.0459,    142.6926,   5.91814131326e-08,  'RC01C',  15868464.006681547, 42.92982544353731,  330.5364368375782],
[48.6577,   7.8152, 142.7607,   3.93093117343e-07,  'R4B4A',  15874660.529874427, 33.25702570840153,  337.9001004221063],
[48.6396,   7.7447, 142.8014,   5.59079951175e-08,  'R60B1',  15879189.883457486, 33.1661588739929,   337.94807612388104],
[48.5135,   7.6979, 142.9237,   3.20104435716e-07,  'RA7C1',  15892816.31329672,  33.18718121640123,  337.877203766338],
[47.4685,   8.9568, 143.3026,   3.76069604302e-08,  'R7B2E',  15935347.134961065, 35.819513671601946, 335.7417032994448],
[46.7748,   10.1561,    143.3614,   9.0894392508e-08,   'R3EDD',  15942228.94002849,  38.121506197903884, 333.9564403829112],
[47.3423,   8.5896, 143.5503,   3.38895483384e-08,  'R03EA',  15962889.109817965, 35.37698095181571,  335.95751799710655],
[47.8108,   7.3379, 143.6426,   7.88474960924e-08,  'R79D9',  15972870.846218215, 33.16126376874516,  337.57289695680373],
[47.6847,   7.4818, 143.6943,   6.89762168935e-08,  'R3774',  15978660.330086328, 33.469777065949536, 337.3204837508338],
[47.6937,   7.3481, 143.7365,   5.32719718545e-08,  'RF5BC',  15983344.958010405, 33.262782150363186, 337.45550813108116],
[47.6486,   7.2625, 143.8059,   7.43248010302e-08,  'R9DBD',  15991061.48833064,  33.16758638776435,  337.4953321114253],
[47.6126,   7.2303, 143.8479,   3.53749610299e-07,  'R9D74',  15995735.964670816, 33.14577901869697,  337.4928599005947],
[43.5946,   13.5099,    144.1286,   2.38406744292e-07,  'RF7E5',  16028986.757775346, 45.75965395991925,  327.40367424975835],
[46.6847,   7.683,  144.444,    1.84838254678e-07,  'R0C73',  16062300.992469292, 34.53495626041903,  336.1799603274634],
[50.2523,   -5.0309,    144.5607,   2.51145293762e-07,  'R7FA5',  16073291.109304361, 12.281288488578284, 351.87865795451887],
[50.2072,   -5.3075,    144.642,    6.40721692254e-08,  'R480A',  16082321.111570034, 11.84150678184185,  352.15963946080495],
[50.1081,   -5.1702,    144.7205,   2.64242225928e-07,  'R303A',  16091069.94561483,  12.095280303645783, 351.97633823192405],
[46.3063,   7.5646, 144.8011,   2.59958774876e-07,  'R2AEB',  16102087.825666783, 34.65440524888269,  335.9265473703645],
[46.5225,   6.5733, 145.0048,   2.74433362573e-08,  'RD14A',  16124529.877269402, 32.966298390593025, 337.1196870664719],
[45.5586,   8.158,  145.168,    4.60441850113e-07,  'R976C',  16143155.27180041,  36.168169327888016, 334.5907317093059],
[47.9009,   1.8814, 145.3455,   8.36002571101e-08,  'S7A70',  16161490.230750268, 24.604960440739557, 343.1485930842906],
[47.8288,   1.9324, 145.3965,   6.26512975665e-07,  'R51FD',  16167182.75630816,  24.730229658813684, 343.0415486458332],
[44.045,    10.0522,    145.5144,   1.19852495572e-07,  'RF212',  16182414.243399696, 40.360204199417744, 331.09795579490924],
[44.7838,   7.9972, 145.8573,   9.17646713394e-08,  'R4FB1',  16219982.196730752, 36.58369347328169,  333.94446752315616],
[44.7838,   6.8674, 146.3284,   8.8959736833e-08,   'R7CE4',  16272168.531511128, 34.82697852931433,  335.10896602059364],
[45.1802,   5.726,  146.4516,   1.97668119373e-08,  'R7A15',  16285573.409416324, 32.69184757052468,  336.71352442451627],
[45.7477,   3.0599, 146.9137,   1.78557865177e-07,  'R86F8',  16336383.342223246, 27.92976715054852,  340.15944210589885],
[44.2072,   5.0778, 147.519,    1.37479772081e-07,  'RC7DF',  16404354.204405338, 32.43963674649258,  336.46912822273646],
[15.5045,   20.4563,    148.3249,   2.34908562585e-08,  'RA77C',  16508080.183242608, 94.08223112632622,  274.724478184751],
[44.8468,   0.2542, 148.5904,   2.28835896918e-07,  'RB1D6',  16522453.564030897, 23.765449904187925, 342.7418089265241],
[43.3063,   -0.3591,    150.1771,   7.02172507817e-08,  'R7F64',  16698921.543941455, 23.696663499459525, 342.3232616513428],
[42.8739,   0.0123, 150.4587,   8.19871805497e-08,  'R0CA8',  16730349.43788516,  24.684836419617433, 341.47249041736734],
[43.2613,   -2.9321,    150.9065,   1.77826384688e-07,  'R12EB',  16779484.17465724,  18.912428472926837, 345.8141881933233],
[43.3694,   -4.115, 151.0687,   6.19362193259e-08,  'R9081',  16797301.650472563, 16.586000795386173, 347.55818235728253],
[43.3604,   -5.9107,    151.4153,   7.75969544797e-08,  'R0D06',  16835553.593549155, 13.092527692671533, 350.1542839395658],
[42.3243,   -8.652, 152.8045,   7.3528983592e-08,   'R29F5',  16989607.780029785, 7.90396847118016,   353.9409527560108],
[40.1982,   -8.4452,    154.884,    4.44281770378e-08,  'RC085',  17220640.77622039,  8.969649704782185,  352.8988608463395]
]
'''
# Now you can add markers to show each station in turn
# station is a simple list showing the stationID, location, lat, long, for each station in turn
#for station in STATIONS:
#    folium.Marker(location=[station[0], station[1]], popup=station[4], icon=folium.Icon(color='orange')).add_to(map)
# Finally, add a red marker for UDDGP, the deepest borehole in mainland UK
#folium.Marker(location=[47.8288, 1.9324], popup='Caustic', icon=folium.Icon(color='red')).add_to(map)

folium.Marker(location=[EQLAT, EQLON], popup='Earthquake_location',
    icon=folium.Icon(color='orange')).add_to(map)
folium.Marker(location=[STA_LAT, STA_LON], popup=LOCATION,
    icon=folium.Icon(color='red')).add_to(map)
map.save(FILE_STEM + '-earthquakemap.html')


'''
    if plotit == True:
        folium.Circle(
          location=[bestloc[0],bestloc[1]],
          popup=str(lat) + " " + str(lon) + " " + str(DISTANCE),
          radius=1,
          color='crimson',
          fill=True,
          fill_color='crimson'
           ).add_to(map)
    if plotit == True:
        folium.Circle(
          location=[bestloc[0],bestloc[1]+360],
          popup=str(lat) + " " + str(lon) + " " + str(DISTANCE),
          radius=1,
          color='crimson',
          fill=True,
          fill_color='crimson'
           ).add_to(map)
    if plotit == True:
        folium.Circle(
          location=[bestloc[0],bestloc[1]-360],
          popup=str(lat) + " " + str(lon) + " " + str(DISTANCE),
          radius=1,
          color='crimson',
          fill=True,
          fill_color='crimson'
           ).add_to(map)

    if plotit == True:
        folium.Circle(
          location=[bestloc[0],bestloc[1]],
          popup=str(lat) + " " + str(lon) + " " + str(DISTANCE),
          radius=1,
          color='crimson',
          fill=True,
          fill_color='crimson'
           ).add_to(map)
    if plotit == True:
        folium.Circle(
          location=[bestloc[0],bestloc[1]+360],
          popup=str(lat) + " " + str(lon) + " " + str(DISTANCE),
          radius=1,
          color='crimson',
          fill=True,
          fill_color='crimson'
           ).add_to(map)
    if plotit == True:
        folium.Circle(
          location=[bestloc[0],bestloc[1]-360],
          popup=str(lat) + " " + str(lon) + " " + str(DISTANCE),
          radius=1,
          color='crimson',
          fill=True,
          fill_color='crimson'
           ).add_to(map)
'''