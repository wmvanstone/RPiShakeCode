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

# Station details
FILE_STEM = "StDay"
MAPFILE=FILE_STEM + '-map.png'
LOGOS='logos.png'

# Things to change once for your station
NETWORK = 'AM'   # AM = RaspberryShake network
STATION = "RAD67"  # Station code of local station to plot
STA_LAT = 50.2385  # Latitude of local station  
STA_LON = -5.1822  # Longitude of local station
LOCATION = "St Day"
# Plot the epicentre and seismometer on the map
#map = folium.Map(location=[EQLAT, EQLON],zoom_start=2,tiles='Stamen Terrain')
# The following arguments will centre the map on Caustic.
map = folium.Map(location=[50.2385, -5.1822],zoom_start=6,tiles='Stamen Terrain')

# theta is 90 - latitude
# phi is longitude but phi 180-360 becomes negative longitude
# r = radius of earth (6,378 EQ + 6,357 POLE) / 2 = 6367.5
dend = []
dstart = []

for phi in range(0, 360, 1):
    theta = 140
    dend.append(asCartesian([R, theta, phi]))
    omega = 104
    dstart.append(asCartesian([R, omega, phi]))

dendt = []
axis = [math.cos((90+STA_LON)/180*PI), math.sin((90+STA_LON)/180*PI), 0]
theta = (90-STA_LAT)/180*PI
for i in range(len(dend)):
    dend[i]=np.dot(rotation_matrix(axis, theta), dend[i])
    send = asSpherical(dend[i])
    lat = 90 - (send[1])
    long = send[2]
    dendt.append([lat, long])

dstartt = []
axis = [math.cos((90+STA_LON)/180*PI), math.sin((90+STA_LON)/180*PI), 0]
theta = (90-STA_LAT)/180*PI
for i in range(len(dstart)):
    dstart[i]=np.dot(rotation_matrix(axis, theta), dstart[i])
    sstart = asSpherical(dstart[i])
    lat = 90 - (sstart[1])
    long = sstart[2]
    dstartt.append([lat, long])

linepoints=[]
for lon in range(-180, 190, 10):
    # Northern hemisphere
    DISTANCE=locations2degrees(STA_LAT, STA_LON, 0, lon) # Station dist in degrees from epicentre
    bestfit = abs(DISTANCE-140)
    bestloc = [0, lon]
    plotit = False
    for lat in range(0, 91, 1):
        DISTANCE=locations2degrees(STA_LAT, STA_LON, lat, lon) # Station dist in degrees from epicentre
        if DISTANCE > 139.6 and DISTANCE < 140.4:
            plotit = True
            for latf in range((lat-1)*300, (lat+2)*300, 1):
                DISTANCE=locations2degrees(STA_LAT, STA_LON, latf/300, lon)
                fit = abs(DISTANCE-140)
                if fit < bestfit:
                    bestfit = fit
                    bestloc = [latf/300, lon]
    if plotit == True:
        linepoints.append(bestloc)
    # southern hemisphere
        DISTANCE=locations2degrees(STA_LAT, STA_LON, 0, lon) # Station dist in degrees from epicentre
    bestfit = abs(DISTANCE-140)
    bestloc = [0, lon]
    plotit = False
    for lat in range(-90, 1, 1):
        DISTANCE=locations2degrees(STA_LAT, STA_LON, lat, lon) # Station dist in degrees from epicentre
        if DISTANCE > 139.6 and DISTANCE < 140.4:
            plotit = True
            for latf in range((lat-1)*300, (lat+2)*300, 1):
                DISTANCE=locations2degrees(STA_LAT, STA_LON, latf/300, lon)
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

QUAKES = [
['5.8 18km SSE of Lone Pine, CA',36.4468333,-117.9751667],
['5.9 38 km SE of Hasaki, Japan',35.4661,141.1005],
['5.4 11 km SSE of Azalp, Turkey',38.5583,44.0232],
['6.3 278 km SE of Hotan, China',35.5948,82.4158],
['5.2 13 km NW of Galmarmara, Turkey',38.7907,27.794],
['4.5 N Atlantic West of Lorient',47.53,-5.44],
['5.5 13 km S of Marmaris, Turkey',36.7307,28.2514],
['5.1 6 km W of Isangel, Vanuatu',-19.5429,169.2234],
['5.2 105 km E of Shikotan, Russia',43.7226,148.0257],
['5.5 219 km S of Sand Point, Alaska',53.3659,-160.4307],
['5.9 33 km ENE of Port-Olry, Vanuatu',-14.926,167.3646],
['6.7 98 km N of Batang, Indonesia',-5.5956,110.6952],
['5.1 108 km SSE of Akhiok, Alaska',56.0177,-153.6203],
['4.9 138 km ENE of Kurilsk, Russia',45.7491,149.4895],
['4.9 141 km SSW of Hihifo, Tonga',-17.1659,-174.2141],
['4.9 Hokkaido, Japan region',41.0224,143.0792],
['4.6 263 km SSE of Alo, Wallis and Futuna',-16.5685,-177.3411],
['5.4 6 km W of El Crucero, Nicaragua',11.9976,-86.3678],
['5 117 km E of Mutsu, Japan',41.2154,142.6111],
['5.1 73 km NW of Port-Vila, Vanuatu',-17.2821,167.8058],
['4.9 South of the Fiji Islands',-26.5734,178.3776],
['5.8 70 km WSW of Rio Grande, Panama',7.4708,-81.922],
['7 114 km NNW of Popondetta, Papua New Guinea',-7.8428,147.7656],
['5.9 2 km SSE of Iquique, Chile',-20.2353,-70.1399],
['5 51 km SSW of Kurilsk, Russia',44.8392,147.5213],
['5.7 246 km E of Port Blair, India',12.1027,94.9647],
['6 Fiji region',-20.8319,-178.5723],
['5.1 24 km ESE of Albania, Colombia',11.046,-72.3972],
['7.8 105 km SSE of Perryville, Alaska',55.0298,-158.5217],
['6.1 103 km ESE of Sand Point, Alaska',54.9322,-159.0423],
['5.7 123 km SE of Sand Point, Alaska',54.7132,-158.9076],
['5.1 90 km WSW of Port-Olry, Vanuatu',-15.3968,166.3129],
['5.8 72 km E of Sand Point, Alaska',55.2233,-159.3719],
['4.7 98 km SSE of Perryville, Alaska',55.0659,-158.6912],
['5.1 69 km S of Sand Point, Alaska',54.7104,-160.5043],
['6.3 western Xizang',33.1313,86.8397],
['5.5 76 km S of Ivanof Bay, Alaska',55.2191,-159.2945],
['5.1 western Xizang',33.1466,86.7977],
['4.9 56 km ESE of Sand Point, Alaska',55.2309,-159.6341],
['5 56 km SSW of Sand Point, Alaska',54.8383,-160.6905],
['5 47 km SSW of Sand Point, Alaska',54.9333,-160.739],
['5.4 Fiji region',-20.4685,-178.4061],
['4.9 61 km SSE of Sand Point, Alaska',54.8117,-160.2162],
['5.5 26 km S of Sand Point, Alaska',55.0965,-160.4806],
['6.1 66 km SW of Sand Point, Alaska',54.8674,-161.1436],
['5 39 km S of Sand Point, Alaska',54.9841,-160.3834],
['4.8 Fiji region',-18.7505,-177.7718],
['5.6 41 km SSE of Sand Point, Alaska',54.9728,-160.3332],
['4.9 53 km ENE of Oxapampa, Peru',-10.3158,-74.9872],
['4.8 75 km SSW of Sand Point, Alaska',54.7508,-161.0895],
['5 132 km SSE of Perryville, Alaska',54.845,-158.1992],
['5.7 66 km SE of Neiafu, Tonga',-19.0322,-173.4935],
['4.9 Fiji region',-20.6892,-178.4019],
['6.4 10 km WSW of Polloc, Philippines',7.3037,124.1419],
['5.6 211 km SSW of ‘Ohonua, Tonga',-23.0645,-175.8175],
['4.9 53 km E of Luganville, Vanuatu',-15.4489,167.659],
['6.4 Vanuatu',-16.1118,168.0816],
['5.2 88 km SE of King Cove, Alaska',54.5997,-161.2001],
['5.1 11 km NNW of Ojiya, Japan',37.3899,138.7273],
['4.9 59 km SSE of Vilyuchinsk, Russia',52.4764,158.8815],
['5.3 south of Alaska',54.8011,-158.2608],
['5.5 67 km ESE of King Cove, Alaska',54.8389,-161.3373],
['5.1 4 km SE of Sparta, North Carolina',36.4758333,-81.0933333],
['5.2 central Mid-Atlantic Ridge',0.9294,-27.733],
['5.2 120 km SE of Perryville, Alaska',55.0605,-157.9545],
['6 66 km ESE of Vikindu, Tanzania',-7.3283,39.8073],
['5 228 km SSW of Severo-Kuril’sk, Russia',48.81,154.8063],
['5.6 55 km SSW of Surab, Pakistan',28.0478,65.9911],
['5.7 158 km ESE of Akutan, Alaska',53.42,-163.6973],
['2.9 6 km NW of Urville-Nacqueville, France',49.72,-1.77],
['5.7 245 km ESE of Atka, Alaska',51.1005,-171.1262],
['5 51 km SSE of Ýdra, Greece',36.9398,23.7364],
['5 238 km E of Levuka, Fiji',-18.0716,-178.4276],
['5 127 km SE of Perryville, Alaska',55.0003,-157.9167],
['4.7 Fiji region',-21.3275,-179.2414],
['6.9 126 km WSW of Bengkulu, Indonesia',-4.2781,101.228],
['5.1 44 km SSE of Sand Point, Alaska',54.9499,-160.3011],
['6.9 220 km SSE of Katabu, Indonesia',-6.6704,123.4927],
['4.9 52 km NNW of Port-Vila, Vanuatu',-17.3361,168.0459],
['5 54 km NW of Port-Vila, Vanuatu',-17.3278,168.0206],
['5 30 km NW of Rumonge, Burundi',-3.7626,29.2663],
['6 3 km ESE of Jacó, Costa Rica',9.5976,-84.6014],
['5.2 71 km SW of ‘Ohonua, Tonga',-21.826,-175.3925],
['5.5 11 km WNW of Papayal, Peru',-4.0303,-80.8246],
['4.5 279 km ENE of Olonkinbyen, Svalbard and Jan Mayen',71.9229,-1.5186],
['4.8 123 km SW of Leava, Wallis and Futuna',-15.0243,-179.0258],
['5 42 km NE of Mohr, Iran',27.8326,53.1812],
['5 77 km NNW of Port-Vila, Vanuatu',-17.1501,167.9134],
['6.5 central Mid-Atlantic Ridge',0.8696,-29.7046],
['-0.1 PORTREATH,CORNWALL',50.28,-5.308],
['1.7 English Channel 18km S of Plymouth',50.206,-4.128],
['6.8 78 km NW of Vallenar, Chile',-28.0121,-71.2377],
['5.5 65 km SW of Palana, Russia',58.7546,159.0051],
['5.8 Pacific-Antarctic Ridge',-55.4953,-129.7095],
['4.9 32 km SE of Nemuro, Japan',43.1064,145.838],
['3.7 Spain',41.93,3.71],
['6.6 central Mid-Atlantic Ridge',7.6997,-37.2385],
['5.3 northern Mid-Atlantic Ridge',47.7463,-27.4915],
['6.2 101 km NW of Port-Vila, Vanuatu',-17.1497,167.5791],
['5.7 65 km NNE of Port-Vila, Vanuatu',-17.1673,168.4862],
['6.2 73 km NNE of Port-Vila, Vanuatu',-17.0936,168.4969],
['5.9 197 km SSE of Amahai, Indonesia',-4.9124,129.7627],
['5.3 31 km SSE of Karyes, Greece',39.9957,24.3777],
['5.6 127 km SSW of Nikolski, Alaska',51.9341,-169.7665],
['5.6 279 km E of Levuka, Fiji',-17.8809,-178.0547],
['6.3 83 km NNE of Tocopilla, Chile',-21.3928,-69.8943],
['5.9 84 km NW of Port-Vila, Vanuatu',-17.2562,167.6938],
['6.1 57 km SE of Ōfunato, Japan',38.7591,142.2473],
['5.2 209 km NE of Lospalos, Timor Leste',-7.3138,128.4577],
['3.1 France',49.56,-1.55],
['2.7 France',49.62,-1.54],
['2.8 France',49.62,-1.55],
['2.8 France',49.68,-1.54],
['2.7 France',49.6,-1.54],
['2.9 France',49.6,-1.57],
['4.6 Fiji region',-20.7211,-178.3549],
['6.4 13 km WNW of Esso, Russia',55.9569,158.4936],
['4.8 10 km NW of Aguas Verdes, Peru',-3.4104,-80.3066],
['5.3 11 km SW of Kodāri̇̄, Nepal',27.8713,85.8877],
['3.5 France',49.01,-1.64],
['5.6 southeast of the Loyalty Islands',-22.3257,171.4996],
['5.7 central Mid-Atlantic Ridge',7.9762,-37.0142],
['5.2 7 km S of Belisario Domínguez, Mexico',15.2344,-92.3866],
['5.6 259 km SSE of Alo, Wallis and Futuna',-16.5089,-177.2759],
['6.9 central Mid-Atlantic Ridge',0.9454,-26.8495],
['5.9 16 km E of Pýrgos, Greece',34.9972,25.3313],
['4.9 Banda Sea',-6.8617,129.3913],
['3 Leighton Buzzard',51.941,-0.662],
['3.5 Leighton Buzzard',51.926,-0.737],
['5.2 112 km SSW of Nikolski, Alaska',52.0953,-169.7613],
['6.4 37 km NNE of Pangai, Tonga',-19.495,-174.2254],
['4.5 0 km NNW of Huarte-Uharte, Spain',42.836,-1.5935],
['5.1 16 km W of Kéfalos, Greece',36.7633,26.781],
['5.3 72 km SW of Atka, Alaska',51.7326,-174.9395],
['3.9 PYRENEES',42.83,-1.49],
['5.3 eastern Greenland',69.2902,-29.7839],
['6 232 km E of Levuka, Fiji',-18.0249,-178.4906],
['5.9 68 km SE of Sand Point, Alaska',54.8444,-159.8598],
['5.7 47 km ESE of Nikolski, Alaska',52.7184,-168.2558],
['5.2 112 km WNW of Port-Vila, Vanuatu',-17.1897,167.4235],
['4.9 92 km WSW of Türkmenbaşy, Turkmenistan',39.5951,52.021],
['5.2 46 km N of Palekastro, Greece',35.6143,26.2215],
['5.1 45 km N of Palekastro, Greece',35.6061,26.2682],
['5.9 107 km SSE of Sand Point, Alaska',54.4243,-159.955],
['5.5 133 km SSE of Sand Point, Alaska',54.2393,-159.6532],
['7.5 91 km SE of Sand Point, Alaska',54.662,-159.6752],
['5.6 10 km WSW of Hafnarfjörður, Iceland',64.0196,-22.1159],
['5.6 128 km SSE of Sand Point, Alaska',54.2886,-159.6651],
['5.9 107 km SSE of Sand Point, Alaska',54.4243,-159.955],
['5.5 133 km SSE of Sand Point, Alaska',54.2393,-159.6532],
['5.7 92 km SE of Sand Point, Alaska',54.6663,-159.6393],
['7.6 97 km SSE of Sand Point, Alaska',54.6079,-159.6552],
['5.9 183 km ESE of Neiafu, Tonga',-19.3191,-172.3884],
['4.7 Reykjanes Ridge',52.689,-35.1277],
['5.8 148 km WNW of Haveluloto, Tonga',-20.8778,-176.6074],
['5.2 52 km SSW of Lithakiá, Greece',37.3249,20.5069],
['5.1 Reykjanes Ridge',52.6958,-34.9841],
['4.5 153 km SSE of Sand Point, Alaska',54.0495,-159.6577],
['5.1 102 km SSE of Sand Point, Alaska',54.5552,-159.662],
['5.3 Fiji region',-20.8916,-178.5061],
['7 15 km NNE of Néon Karlovásion, Greece',37.9175,26.7901],
]
# Now you can add markers to show each station in turn
# station is a simple list showing the stationID, location, lat, long, for each station in turn
for quake in QUAKES:
    folium.Marker(location=[quake[1], quake[2]], popup=quake[0], icon=folium.Icon(color='orange')).add_to(map)

folium.Marker(location=[STA_LAT, STA_LON], popup=LOCATION,
    icon=folium.Icon(color='red')).add_to(map)
map.save(FILE_STEM + '-earthquakemap.html')