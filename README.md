# Earthquake
Plotting earthquake data 

![alt tag](https://meyavuz.files.wordpress.com/2017/03/earthquakes_dpi240_1303_55-0.png)


<b>To check the code using pylint: </b>

<code> pylint --rcfile=./pylint.rc plotearthquake_post.py</code>

<b> How to run? </b>

<code>./plotearthquake_post.py -markersize 10 -npoints 300 -usgsdata "./USGSData/USGS_Turkey_1900_2017_over6.csv"</code>

To generate the time-lapse, need to follow :
1) Need to obtain snapshots for each time instant, i.e. npoints should be changing from <code>1</code> to <code>length(data)</code> where data is read from the <code>.csv</code> file 
2) Also, need to obtain snapshots for different markersizes so that a zooming in effect is visible, I used markersizes from 5 to 60 with steps of 5, (i.e. <code>5:5:60</code>) 
3) Since the names of the files are dumped accordingly, you can use Blender or any other time-lapse tool to create the final time-lapse
