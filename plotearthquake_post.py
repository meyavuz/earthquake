
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm
 
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.colors import Normalize
import numpy as np

def ParseInput():
  ''' Parse input arguments'''
  import argparse
  parser = argparse.ArgumentParser()

  parser.add_argument("-markersize", type=float, default=5, help="Marker size to represent the earthquake data on the map")
  parser.add_argument("-npoints", type=int, default=1, help="Total number of points in the graph [sweep the whole range for final video output]")
  parser.add_argument("-nsimpoints", type=int, default=1, help="Number of simultaneous points having different marker size")
  parser.add_argument("-usgsdata", type=str, help="Filename (e.g. usgs.csv) that contains the earthquake data as downloaded from USGS")

  args = parser.parse_args()

  return args

def WriteCityNames(m):
  ''' Write names of selected cities on the map'''

  # Lat/lon coordinates of several cities that lie in the map of interest
  lats = [41.00, 41.71, 35.12, 35.24, 37.04, 37.26, 39.90, 
         44.42, 44.78, 41.32, 36.89, 35.46, 31.20, 32.09,
         43.60, 33.89, 39.93, 42.13, 31.94, 45.04, 36.20,
         43.85, 41.11, 31.76, 29.87, 44.61, 38.50, 38.35,
         36.43, 45.65, 42.26, 38.42, 
         42.83]
  lons = [28.97, 44.82, 33.42, 24.80, 22.11, 35.39, 41.26, 
         26.10, 20.44, 19.81, 30.71, 44.38, 29.91, 20.18,
         39.73, 35.50, 32.85, 24.74, 35.92, 41.96, 37.13,
         18.41, 16.87, 35.21, 40.10, 33.52, 43.37, 38.33,
         28.21, 25.60, 42.71, 27.14,
         31.70]
  cities = ['Istanbul','Tblisi','Cyprus','Crete','Kalamata','Adana','Erzurum',
           'Bucharest','Belgrade','Tirana','Antalya','Kerkuk','Alexandria', 'Benghazi',
           'Sochi','Beirut','Ankara','Plovdiv','Amman','Stavropol','Aleppo',
           'Sarajevo','Bari','Jerusalem','Sakaka','Sevastopol','Van','Malatya',
           'Rhodes', 'Brasov','Kutaisi','Izmir',
           'B  L  A  C  K    S  E  A']

  # Compute the native map projection coordinates for cities.
  xc, yc = m(lons, lats)

  # Plot filled circles at the locations of the cities.
  m.plot(xc[:-1],yc[:-1],'bo')

  # Some certain city names need to be shifted for better visualization
  for name, xpt, ypt in zip(cities,xc,yc):
    if name == 'Alexandria' or name == 'Crete' or name == 'Van' or name == 'Malatya':
      plt.text(xpt+10000,ypt-20000,name,fontsize=9)
    elif name == 'Jerusalem':
      plt.text(xpt-40000,ypt-30000,name,fontsize=9)
    elif name == 'Kalamata':
      plt.text(xpt-30000,ypt+15000,name,fontsize=9)
    elif name == 'Istanbul' or name == 'Benghazi':
      plt.text(xpt+20000,ypt-10000,name,fontsize=9)
    elif name == 'B  L  A  C  K    S  E  A':
      plt.text(xpt+15000,ypt+10000,name,fontsize=11)
    else:
      plt.text(xpt+10000,ypt+10000,name,fontsize=9)

  return

def ReadAndGetData(filename):
  ''' 
  Reads the data file downloaded from the USGS site and filters as necessary 
  (e.g. only accept earthquakes with magnitude greater than 5 etc)
  '''

  # Use pandas to get the USGS data content
  df = pd.read_csv(filename)  # e.g. filename:  USGS_americas_1950_2017_over6.csv

  # Use hardcoded min/max Longitude/Latitude for the specific map output
  minLon, maxLon = (18.81  ,  51.327  )
  minLat, maxLat = (  29.155  ,  47.883  )
  
  # Filter data & get only thos larget than 5.0 magnitude
  minMagnitudeDesired = 5.0

  filteredDf = df[ (df['latitude'] >= minLat) & 
                   (df['latitude'] <= maxLat) &
                   (df['longitude'] >= minLon-3.4) &
                   (df['longitude'] <= maxLon-0.7) &
                   (df['mag'] >= minMagnitudeDesired) ]

  filteredDf = df
  #print "len(new.df) = ", len(filteredDf)

  # Filtered data is in panda dataframe format, use df.values.tolist() to convert to list  
  #lat = df.latitude  # Unfiltered
  #lon = df.longitude
  lat = (filteredDf.latitude).values.tolist()
  lon = (filteredDf.longitude).values.tolist()
  date = (filteredDf.time).values.tolist()
  magn = (filteredDf.mag).values.tolist()

  #print " lat = ", lat, " len(lat) = ", len(lat)
  #print " lat[0] =", lat[0]
  
  # If you don't want filtering, then min/max can be obtained from the read data
  #maxLon = max(lon)
  #minLon = min(lon)
  #maxLat = max(lat)
  #minLat = min(lat)

  midLat = 0.5*(maxLat + minLat)
  midLon = 0.5*(maxLon + minLon)

  print "Number of points = ", len(lon)
  print "Latitude (min,max)  = ( ", minLat, " , ", maxLat," )"
  print "Longitude (min,max) = ( ", minLon, " , ", maxLon," )"
  print "Mid (Longitude, Langitude)= ( ", midLon, " , ", midLat," )"

  return (lat, lon, date, magn, minLon, maxLon, minLat, maxLat, midLat, midLon)

def UseGMPLOTtoDumptoGoogleMap(lat, lon, midLat, midLon):
  ''' Convert the same earthquake data info to Google Map heat map format '''

  import gmplot
  gmap = gmplot.GoogleMapPlotter(midLat, midLon, 3)
  print "gmap = ", gmap

  #exit()
  #gmap.plot(latitudes, longitudes, 'cornflowerblue', edge_width=10)
  #gmap.plot(latitudes, longitudes, 'red', edge_width=8)
  #gmap.scatter(more_lats, more_lngs, '#3B0B39', size=40, marker=False)
  #gmap.scatter(marker_lats, marker_lngs, 'k', marker=True)
  #gmap.scatter(lat, lon, 'r', size=10, marker=False)

  gmap.heatmap(lat, lon)

  gmap.draw("earthquake_test.html")

  return

def PlotEarthquakeLocationsOnMap(m, lon, lat, plt, date, magn):

    # Default size for already displayed points
    pstart = ARGS.npoints - ARGS.nsimpoints
    pend = ARGS.npoints
    x, y = m(lon[0:pstart],lat[0:pstart])
    m.plot(x, y,  'ro', alpha=0.8, markersize=5 , markeredgecolor='red', 
              fillstyle='full', markeredgewidth=0.1)

    # Custom (most of the time bigger) font for the new point to be displayed
    x, y = m(lon[pstart:pend],lat[pstart:pend])
    m.plot(x, y,  'ro', alpha=0.8, markersize=65-ARGS.markersize , markeredgecolor='red', 
              fillstyle='full', markeredgewidth=0.1)

    day = (date[ARGS.npoints].split('T'))[0]
    #plt.title('Earthquakes above magnitude ' + str(minMagnitudeDesired)  - ' + day)
    plt.title('Earthquake on ' + day + ' - magnitude ' + str(magn[pend]))
    #plt.title('Earthquakes')

    return


def SaveSnapshotsToFile(bSaveFigs = True ):
  ''' Save the earthquale snapshots in time to file'''

  if bSaveFigs:
    OutFolder = 'Snapshots_'+str(ARGS.npoints-ARGS.npoints%100)
    if not os.path.exists(OutFolder):
      os.mkdir(OutFolder)

    plt.savefig(OutFolder+'/earthquakes_dpi240_'+ 
                str(ARGS.npoints)+'_'+ str(ARGS.markersize) +'.png',
                facecolor='w',dpi=240) 
      
  return

# TODO move midLat/midLon to inside here 
def DrawMap(minLat, maxLat, minLon, maxLon, midLat, midLon):

  #exit()
  #print "lat =", lat
  #print "lon =", lon
  #fig, ax = plt.subplots(figsize=(20,10))
  #fig, ax = plt.subplots(figsize=(15.6,11.4))
  fig, ax = plt.subplots(figsize=(16,9))
  fig.patch.set_facecolor('white') # Set white background

  ##m = Basemap(resolution='c', # c(crude), l(low), i)intermediate), h(high), f(full) or None
  #m = Basemap(resolution='h', # c(crude), l(low), i)intermediate), h(high), f(full) or None
              #projection='merc',
              ##lat_0=40.320373, lon_0=-74.43,
              ##llcrnrlon=-75.00, llcrnrlat= 40.0000, urcrnrlon=-72.00, urcrnrlat=42.000 )
              #lat_0=40.320373, lon_0=-74.43,
              #llcrnrlon=minLon, llcrnrlat= minLat, urcrnrlon=maxLon, urcrnrlat=maxLat )

  #m = Basemap(resolution='h', # c(crude), l(low), i)intermediate), h(high), f(full) or None
              #projection='ortho',
              #lat_0=midLat, lon_0=midLon)


  #m = Basemap(resolution='h', # c(crude), l(low), i)intermediate), h(high), f(full) or None
              #projection='gnom',
              #width=15.e6,height=15.e6,
              #lat_0=midLat, lon_0=midLon)

  #m = Basemap(height=1.9e6, width=3.8e6,
  #m = Basemap(height=1.9e6, width=2.6e6,
  m = Basemap(height=1.7e6, width=2.8e6,
            resolution='f',area_thresh=10.,projection='omerc',\
            lon_0=midLon,lat_0=midLat, \
            lon_1=minLon,lat_1=minLat, lon_2=maxLon,lat_2=maxLat)

  #m = Basemap(projection='mill', resolution='h',\
        #llcrnrlon=minLon, llcrnrlat= minLat, urcrnrlon=maxLon, urcrnrlat=maxLat, \
        #epsg = 4269)
  #m.arcgisimage(service='World_Physical_Map', xpixels = 5000, verbose= False)
    

  m.drawmapboundary(fill_color='#46bcec')
  m.fillcontinents(color='#f2f2f2',lake_color='#46bcec')
  m.drawcounties()
  m.drawcountries(linewidth=0.25)

  m.drawcoastlines()

  m.drawparallels(np.arange(minLat,maxLat,10.))
  m.drawmeridians(np.arange(minLon,maxLon,10.))

  return m


def main():

  # Refer this page: http://www.datadependence.com/2016/06/creating-map-visualisations-in-python/
  # Data downloaded from https://earthquake.usgs.gov/earthquakes/search/
  # A small tutorial : https://peak5390.wordpress.com/2012/12/08/matplotlib-basemap-tutorial-plotting-global-earthquake-activity/


  lat, lon, date, magn, minLon, maxLon, minLat, maxLat, midLat, midLon = ReadAndGetData(ARGS.usgsdata)

  m = DrawMap(minLat, maxLat, minLon, maxLon, midLat, midLon)

  WriteCityNames(m)

  bPlotPoints = True
  if bPlotPoints:
    PlotEarthquakeLocationsOnMap(m, lon, lat, plt, date, magn)

  plt.tight_layout()


  SaveSnapshotsToFile(True)
    
  #plt.show()

  # This is additional stuff - dumping the heatmap to google maps format 
  UseGMPLOTtoDumptoGoogleMap(lat, lon, midLat, midLon)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':

  ARGS = ParseInput()

  main()
