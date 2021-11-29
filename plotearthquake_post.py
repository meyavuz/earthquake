''' Plots selected earthquake data on the map of the world '''

import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import gmplot

# Update 2021_11: Basemap has been deprecated and not possible to install it on the Mac
# The recommendation is to utilize the Cartopy library instead. The installation of it was
# also problematic with pip. Recommendation is to utilize conda.
# However the code needs to be updated from Basemap to Cartopy which will take some time.

class EarthquakeData(object):
    ''' Class storing & manipulating earthquake data '''

    def __init__(self):
        self.latitude = []
        self.longitude = []
        self.minLatitude = 0.0
        self.maxLatitude = 0.0
        self.minLongitude = 0.0
        self.maxLongitude = 0.0
        self.midLatitude = 0.0
        self.midLongitude = 0.0
        self.date = ''
        self.magnitude = 0.0
        self.map = Basemap()

    def ReadAndGetData(self, filename):
        '''
        Reads the data file downloaded from the USGS site and filters as necessary
        (e.g. only accept earthquakes with magnitude greater than 5 etc)
        '''

        # Use pandas to get the USGS data content
        df = pd.read_csv(filename)  # e.g. filename:  USGS_americas_1950_2017_over6.csv

        # Use hardcoded min/max Longitude/Latitude for the specific map output
        # One can also use other min/max values (up to the user)
        self.minLongitude, self.maxLongitude = (18.81, 51.327)
        self.minLatitude, self.maxLatitude = (29.155, 47.883)

        # Filter data & get only thos larget than 5.0 magnitude
        minMagnitudeDesired = 5.0

        filteredDf = df[(df['latitude'] >= self.minLatitude) &
                        (df['latitude'] <= self.maxLatitude) &
                        (df['longitude'] >= self.minLongitude - 3.4) &
                        (df['longitude'] <= self.maxLongitude - 0.7) &
                        (df['mag'] >= minMagnitudeDesired)]

        filteredDf = df

        # Filtered data is in panda dataframe format, use df.values.tolist() to convert to list
        #self.latitude = df.latitude  # Unfiltered
        #self.longitude = df.longitude
        self.latitude = (filteredDf.latitude).values.tolist()
        self.longitude = (filteredDf.longitude).values.tolist()
        self.date = (filteredDf.time).values.tolist()
        self.magnitude = (filteredDf.mag).values.tolist()

        # If you don't want filtering, then min/max can be obtained from the read data
        #self.minLongitude = min(self.longitude)
        #self.maxLongitude = max(self.longitude)
        #self.minLatitude = min(self.latitude)
        #self.maxLatitude = max(self.latitude)

        # Calculate mid-point for longitudes and latitudes to center the map upon
        self.midLatitude = 0.5*(self.maxLatitude + self.minLatitude)
        self.midLongitude = 0.5*(self.maxLongitude + self.minLongitude)

        print("Quake info: \n", self.__str__())


    def DrawMap(self):
        ''' Draw the main background map where the earthquake data will be overlayed

        Some other options to check for:
        self.map = Basemap(resolution='h', # c(crude), l(low), i)intermediate), h(high), f(full) or None
                    projection='merc', # 'ortho', 'gnom', 'mill'
                    lat_0=40.320373, lon_0=-74.43,
                    llcrnrlon=minLon, llcrnrlat= minLat, urcrnrlon=maxLon, urcrnrlat=maxLat )

        # It is also possible to download arcgis images through the following command
        self.map.arcgisimage(service='World_Physical_Map', xpixels = 5000, verbose= False)
        '''

        fig, ax = plt.subplots(figsize=(16, 9))
        fig.patch.set_facecolor('white') # Set white background

        self.map = Basemap(height=1.7e6, width=2.8e6,
                  resolution='f', area_thresh=10., projection='omerc',
                  lon_0=self.midLongitude, lat_0=self.midLatitude,
                  lon_1=self.minLongitude, lat_1=self.minLatitude,
                  lon_2=self.maxLongitude, lat_2=self.maxLatitude)

        self.map.drawmapboundary(fill_color='#46bcec')
        self.map.fillcontinents(color='#f2f2f2', lake_color='#46bcec')
        self.map.drawcounties()
        self.map.drawcountries(linewidth=0.25)

        self.map.drawcoastlines()

        self.map.drawparallels(np.arange(self.minLatitude, self.maxLatitude, 10.))
        self.map.drawmeridians(np.arange(self.minLongitude, self.maxLongitude, 10.))

        plt.tight_layout()

        return

    def PlotEarthquakeLocationsOnMap(self, bPlotPoints):
        ''' Just plot the points where the earthquake occurred '''

        if bPlotPoints:
            # Default size for already displayed points (in a time-lapse fashion, some points
            # have already been displayed as shrinking points - at this stage only show shrunk versions)
            pstart = ARGS.npoints - ARGS.nsimpoints
            pend = ARGS.npoints
            x, y = self.map(self.longitude[0:pstart], self.latitude[0:pstart])
            self.map.plot(x, y, 'ro', alpha=0.8, markersize=5, markeredgecolor='red',
                          fillstyle='full', markeredgewidth=0.1)

            # Custom (most of the time bigger) font for the new point to be displayed
            x, y = self.map(self.longitude[pstart:pend], self.latitude[pstart:pend])
            self.map.plot(x, y, 'ro', alpha=0.8, markersize=65-ARGS.markersize, markeredgecolor='red',
                          fillstyle='full', markeredgewidth=0.1)

            day = (self.date[ARGS.npoints].split('T'))[0]
            magnitude = self.magnitude[pend]
            plt.title('Earthquake on {} - magnitude {}'.format(day, magnitude))

        return


    def WriteCityNamesOnTheMap(self):
        '''
        Write names of selected cities on the map after finding
        their corresponding latitude/longitude value
        '''

        # Lat/lon coordinates of several cities that lie in the map of interest
        lats = [41.00, 41.71, 35.12, 35.24, 37.04, 37.26, 39.90,
              44.42, 44.78, 41.32, 36.89, 35.46, 31.20, 32.09,
              43.60, 33.89, 39.93, 42.13, 31.94, 45.04, 36.20,
              43.85, 41.11, 31.76, 29.87, 44.61, 38.50, 38.35,
              36.43, 45.65, 42.26, 38.42, 42.83]

        lons = [28.97, 44.82, 33.42, 24.80, 22.11, 35.39, 41.26,
              26.10, 20.44, 19.81, 30.71, 44.38, 29.91, 20.18,
              39.73, 35.50, 32.85, 24.74, 35.92, 41.96, 37.13,
              18.41, 16.87, 35.21, 40.10, 33.52, 43.37, 38.33,
              28.21, 25.60, 42.71, 27.14, 31.70]

        cities = ['Istanbul', 'Tblisi', 'Cyprus', 'Crete', 'Kalamata', 'Adana',
                  'Erzurum', 'Bucharest', 'Belgrade', 'Tirana', 'Antalya',
                  'Kerkuk', 'Alexandria', 'Benghazi', 'Sochi', 'Beirut',
                  'Ankara', 'Plovdiv', 'Amman', 'Stavropol', 'Aleppo',
                  'Sarajevo', 'Bari', 'Jerusalem', 'Sakaka', 'Sevastopol',
                  'Van', 'Malatya', 'Rhodes', 'Brasov', 'Kutaisi', 'Izmir',
                  'B  L  A  C  K    S  E  A']

        # Compute the native map projection coordinates for cities.
        xc, yc = self.map(lons, lats)

        # Plot filled circles at the locations of the cities.
        self.map.plot(xc[:-1], yc[:-1], 'bo')

        # Some certain city names need to be shifted for better visualization
        for name, xpt, ypt in zip(cities, xc, yc):
            if name == 'Alexandria' or name == 'Crete' or name == 'Van' or name == 'Malatya':
                plt.text(xpt+10000, ypt-20000, name, fontsize=9)
            elif name == 'Jerusalem':
                plt.text(xpt-40000, ypt-30000, name, fontsize=9)
            elif name == 'Kalamata':
                plt.text(xpt-30000, ypt+15000, name, fontsize=9)
            elif name == 'Istanbul' or name == 'Benghazi':
                plt.text(xpt+20000, ypt-10000, name, fontsize=9)
            elif name == 'B  L  A  C  K    S  E  A':
                plt.text(xpt+15000, ypt+10000, name, fontsize=11)
            else:
                plt.text(xpt+10000, ypt+10000, name, fontsize=9) # Default visualization

        return

    def UseGMPLOTtoDumptoGoogleMap(self, htmlfilename):
        ''' Convert the same earthquake data info to Google Map heat map format
            NOTE: gmplot package needs to be pre-installed

        # Some other options using gmplot
        gmap.plot(latitudes, longitudes, 'cornflowerblue', edge_width=10)
        gmap.plot(latitudes, longitudes, 'red', edge_width=8)
        gmap.scatter(more_lats, more_lngs, '#3B0B39', size=40, marker=False)
        gmap.scatter(marker_lats, marker_lngs, 'k', marker=True)
        gmap.scatter(lat, lon, 'r', size=10, marker=False)
        '''

        gmap = gmplot.GoogleMapPlotter(self.midLatitude, self.midLongitude, 3) # lat/lon/google map zoom level

        gmap.heatmap(self.latitude, self.longitude)

        gmap.draw(htmlfilename)

        return

    def SaveSnapshotsToFile(self, bSaveFigs=True, snapshotfilename='earthquakes_dpi240'):
        '''
        Save the earthquake snapshots in time to file
        Note that for (16,9) sized figure, dpi=120 gives (16,9)*120 =[1920,1080] pixels png file
        Similarly, dpi=240 gives (16,9)*240 =[3840,2160] pixels png file
        '''

        if bSaveFigs:
            OutFolder = 'Snapshots_{}'.format(ARGS.npoints - ARGS.npoints%100)
            if not os.path.exists(OutFolder):
                os.mkdir(OutFolder)

            outpngfilename = '{}/{}_{}_{}.png'.format(OutFolder, snapshotfilename,
                                                      ARGS.npoints, ARGS.markersize)
            plt.savefig(outpngfilename, facecolor='w', dpi=240)
        else:
            plt.show()

        return


    def __str__(self):
        ''' Print some data on the earthquake class'''

        str1 = "Number of points = {}".format(len(self.longitude))
        str2 = "Latitude (min,max) = {}, {}".format(self.minLatitude, self.maxLatitude)
        str3 = "Longitude (min,max) = {}, {}".format(self.minLongitude, self.maxLongitude)
        str4 = "Mid points (Longitude, Langitude)= {}, {}".format(self.midLongitude, self.midLatitude)

        return '{}\n{}\n{}\n{}\n'.format(str1, str2, str3, str4)


def ParseInput():
    ''' Parse input arguments '''

    parser = argparse.ArgumentParser()

    parser.add_argument("-markersize", type=float, default=5, help="Marker size to represent the earthquake data on the map")
    parser.add_argument("-npoints", type=int, default=1, help="Total number of points in the graph [sweep the whole range for final video output]")
    parser.add_argument("-nsimpoints", type=int, default=1, help="Number of simultaneous points having different marker size")
    parser.add_argument("-usgsdata", type=str, help="Filename (e.g. usgs.csv) that contains the earthquake data as downloaded from USGS")

    args = parser.parse_args()

    return args


def testing():
    # Prettymaps
    from prettymaps import *
    # Vsketch
    # import vsketch
    # OSMNX
    import osmnx as ox
    # Matplotlib-related
    import matplotlib.font_manager as fm
    from matplotlib import pyplot as plt
    from descartes import PolygonPatch
    # Shapely
    from shapely.geometry import *
    from shapely.affinity import *
    from shapely.ops import unary_union

    def prettymaps_example():
        """
        Experimental:
        https://nbviewer.org/github/marceloprates/prettymaps/blob/main/notebooks/examples.ipynb
        """

        # Init matplotlib figure
        fig, ax = plt.subplots(figsize = (12, 12), constrained_layout = True)

        backup = plot(
            # Address:
            # 'Praça Ferreira do Amaral, Macau',
            'Çengelköy Mahallesi, Üsküdar, Istanbul, Marmara Region, 34680, Turkey',
            # Plot geometries in a circle of radius:
            radius = 1100,
            # Matplotlib axis
            ax = ax,
            # Which OpenStreetMap layers to plot and their parameters:
            layers = {
                    # Perimeter (in this case, a circle)
                    'perimeter': {},
                    # Streets and their widths
                    'streets': {
                        'width': {
                            'motorway': 5,
                            'trunk': 5,
                            'primary': 4.5,
                            'secondary': 4,
                            'tertiary': 3.5,
                            'residential': 3,
                            'service': 2,
                            'unclassified': 2,
                            'pedestrian': 2,
                            'footway': 1,
                        }
                    },
                    # Other layers:
                    #   Specify a name (for example, 'building') and which OpenStreetMap tags to fetch
                    'building': {'tags': {'building': True, 'landuse': 'construction'}, 'union': False},
                    'water': {'tags': {'natural': ['water', 'bay']}},
                    'green': {'tags': {'landuse': 'grass', 'natural': ['island', 'wood'], 'leisure': 'park'}},
                    'forest': {'tags': {'landuse': 'forest'}},
                    'parking': {'tags': {'amenity': 'parking', 'highway': 'pedestrian', 'man_made': 'pier'}}
                },
                # drawing_kwargs:
                #   Reference a name previously defined in the 'layers' argument and specify matplotlib parameters to draw it
                drawing_kwargs = {
                    'background': {'fc': '#F2F4CB', 'ec': '#dadbc1', 'hatch': 'ooo...', 'zorder': -1},
                    'perimeter': {'fc': '#F2F4CB', 'ec': '#dadbc1', 'lw': 0, 'hatch': 'ooo...',  'zorder': 0},
                    'green': {'fc': '#D0F1BF', 'ec': '#2F3737', 'lw': 1, 'zorder': 1},
                    'forest': {'fc': '#64B96A', 'ec': '#2F3737', 'lw': 1, 'zorder': 1},
                    'water': {'fc': '#a1e3ff', 'ec': '#2F3737', 'hatch': 'ooo...', 'hatch_c': '#85c9e6', 'lw': 1, 'zorder': 2},
                    'parking': {'fc': '#F2F4CB', 'ec': '#2F3737', 'lw': 1, 'zorder': 3},
                    'streets': {'fc': '#2F3737', 'ec': '#475657', 'alpha': 1, 'lw': 0, 'zorder': 3},
                    'building': {'palette': ['#FFC857', '#E9724C', '#C5283D'], 'ec': '#2F3737', 'lw': .5, 'zorder': 4},
                }
        )

        plt.show()
        return


    import plotly.graph_objects as go
    import plotly.express as px

    def plotly_maps():
        """
        This is a nice one with actual open street map with upload/download
        Need to improve this.
        https://medium.com/mlearning-ai/earthquakes-data-visualization-5da0d2abeea1
        """

        def fetch_date(time):
            return str(time).split(' ')[0]

        def fetch_weekday(time):
            date = fetch_date(time)
            return date + ' - ' + str(time.weekday())

        def fetch_hour(time):
            t = str(time).split(' ')
            return t[0] + ' - ' + t[1].split(':')[0]

        def subarea_extractor(x):
            return x[0]

        def area_extractor(x):
            return x[-1]

        def fetch_quakes_data(manner='daily', bbox='Worldwide', mag_thresh=2):
            url_link = 'https://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/{}.csv'

            if (manner == 'weekly'):
                data_link = url_link.format('all_week')
            elif (manner == 'monthly'):
                data_link = url_link.format('all_month')
            else:
                data_link = url_link.format('all_day')

            eqdf = pd.read_csv(filepath_or_buffer=data_link)
            eqdf = eqdf[['time', 'latitude', 'longitude', 'mag', 'place']]

            place_list = eqdf['place'].str.split(', ')
            eqdf['sub_area'] = place_list.apply(func=subarea_extractor)
            eqdf['area'] = place_list.apply(func=area_extractor)
            eqdf = eqdf.drop(columns=['place'], axis=1)

            if isinstance(mag_thresh, int) and (mag_thresh > 0):
                eqdf = eqdf[eqdf['mag'] >= mag_thresh]
            else:
                eqdf = eqdf[eqdf['mag'] > 0]

            if bbox in eqdf['area'].to_list():
                eqdf = eqdf[eqdf['area'] == bbox]
                max_mag = eqdf['mag'].max()
                center_lat = eqdf[eqdf['mag'] == max_mag]['latitude'].values[0]
                center_lon = eqdf[eqdf['mag'] == max_mag]['longitude'].values[0]
            else:
                center_lat, center_lon = [54, 15]

            eqdf['time'] = pd.to_datetime(eqdf['time'])

            if (manner == 'weekly'):
                new_col = 'weekday'
                eqdf[new_col] = eqdf['time'].apply(func=fetch_weekday)
            elif (manner == 'monthly'):
                new_col = 'date'
                eqdf[new_col] = eqdf['time'].apply(func=fetch_date)
            else:
                new_col = 'hours'
                eqdf[new_col] = eqdf['time'].apply(func=fetch_hour)

            eqdf = eqdf.sort_values(by='time')

            return eqdf, center_lat, center_lon

        def visualize_quakes_data(manner='daily', bbox='Worldwide', mag_thresh=2):
            eqdf, clat, clon = fetch_quakes_data(manner=manner, bbox=bbox, mag_thresh=mag_thresh)

            if (manner == 'monthly'):
                af = 'date'
            elif (manner == 'weekly'):
                af = 'weekday'
            else:
                af = 'hours'

            zoom = 3 if bbox != 'Worldwide' else 1

            fig = px.scatter_mapbox(
                data_frame=eqdf,
                lat='latitude',
                lon='longitude',
                center=dict(lat=clat, lon=clon),
                size='mag',
                color='mag',
                hover_name='sub_area',
                zoom=zoom,
                mapbox_style='stamen-terrain',
                # color_continuous_scale=px.colors.cyclical.IceFire,
                animation_frame=af,
                title='{} Earthquakes - {}; mag - {}'.format(manner.title(), bbox, mag_thresh)
            )

            fig.update_layout(
                margin=dict(l=0, r=0, t=40, b=0)
            )

            fig.show()
            return None

        visualize_quakes_data(manner='daily', mag_thresh=None)
        plt.show()
        return


def main():

    quake = EarthquakeData();
    quake.ReadAndGetData(ARGS.usgsdata)

    quake.DrawMap()
    quake.WriteCityNamesOnTheMap()
    quake.PlotEarthquakeLocationsOnMap(True)
    quake.SaveSnapshotsToFile(True)

    # This is a bonus feature - dumping the earthquake heatmap to google maps format [uses gmplot]
    quake.UseGMPLOTtoDumptoGoogleMap("earthquake_test.html")

    # Experimental - Basemap is deprecated. Until that is replaced by Cartopy or others
    # We will have the following functions for testing.
    # return
    # testing()
    # prettymaps_example()
    # plotly_maps()



if __name__ == '__main__': # standard boilerplate calling main()
    ARGS = ParseInput()
    main()
