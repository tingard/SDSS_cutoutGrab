## SDSS_cutoutGrab

This is a short Python script containing functions to get a fits file from the SDSS skyserver for an object with a given Ra and Dec.

###Contains:

- `getFitsFromCoord`

This function accepts a right ascension and declination in degrees, as well as a search threshold, and queries [http://skyserver.sdss.org/dr13/en/tools/search/](http://skyserver.sdss.org/dr13/en/tools/search/) for matching objects. It takes the run, camcol, and field of the returned object (if none or more than one are returned it adjusts the threshold accordingly and queries again) and obtains the relevant observation fits file from:

[http://data.sdss.org/sas/dr13/eboss/photoObj/frames/301/{run}/{camcol}/frame-{band}-{run:06d}-{camcol}-{field:04d}.fits.bz2](http://data.sdss.org/sas/dr13/eboss/photoObj/frames/301/{run}/{camcol}/frame-{band}-{run:06d}-{camcol}-{field:04d}.fits.bz2)

- `cutFits`

This function accepts a file path to a fits image, a right ascension and declination in degrees in the same system as the fits file being cut, and a pixel size (int or length-2 list/tuple. It opens the fits image and attempts to crop it at the given Ra and Dec to the given size.


### TODOs:
- better error handling
- better logging if desired
- accept command line arguments
- output cropped fits file?