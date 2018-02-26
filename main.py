import os
import requests
import json
import bz2
import re
import subprocess
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord, Angle
import astropy.units as u
from astropy.nddata import Cutout2D
from astropy.wcs import WCS

loc = os.path.abspath(os.path.dirname(__file__)) + '/'
def printWarning(s): print('\033[33m[WARNING] ' + s + '\033[0m')

def downloadFile(url, fname):
    global loc
    # check if we already have the file in the output file location
    if os.path.isfile(loc + fname):
        return True
    # check if directory tree needs updating
    if not os.path.isdir(os.path.abspath(loc+fname)):
        fstruc = fname.split('/')
        for i in range(len(fstruc)):
            if not os.path.isdir(loc + '/'.join(fstruc[:i])):
                os.mkdir(loc + '/'.join(fstruc[:i]))
    # start the download
    r = requests.get(url, stream=True)
    decompressor = bz2.BZ2Decompressor()
    with open(loc+fname, 'wb') as imgf:
        for chunk in r:
            dChunk = decompressor.decompress(chunk)
            imgf.write(dChunk)
    return True

def getFitsFromCoord(ra=None, dec=None, threshold=5E-6, **kwargs):
    # if we haven't been given any coordinates, return
    if ra == None or dec == None: return None
    dataUrl = "http://data.sdss.org/sas/dr13/eboss/photoObj/frames/301/{path}"
    sqlUrl = "http://skyserver.sdss.org/dr13/en/tools/search/x_sql.aspx?format=json&cmd={sql}"
    # choose the coords range to search in
    ramin = (ra - threshold) % 360
    ramax = (ra + threshold) % 360
    decmin = dec - threshold
    decmax = dec + threshold
    # define the sql to send to skyserver
    sql = "SELECT ra, dec, run, camcol, field FROM photoObjAll WHERE "
    s1 = "({0} between {1} and {2})"
    s2 = "(({0} < {1}) or ({0} > {2}))" # this is used when ramax overflows past 360
    s = "{ra} and {dec}".format(
        ra = (s1 if ramin < ramax else s2).format('ra', ramin, ramax),
        dec = s1.format('dec', decmin, decmax),
    )
    # send the request and read the results in as json (if it's a valid query)
    r = requests.get(sqlUrl.format(sql=sql+s).replace(' ', '%20'))
    if r.status_code != 200:
        return
    returned = json.loads(r.text)
    if len(returned[0]['Rows']) == 1:
        # we have only one row returned, makes life easier
        fitsFiles = []
        for band in ["u","g","r","i","z"]:
            fitsFiles.append(
                '{run}/{camcol}/frame-{band}-{run:06d}-{camcol}-{field:04d}.fits'.format(
                    band=band, **returned[0]['Rows'][0]
                )
            )
        failedFiles = []
        for f in fitsFiles:
            # download the image
            downloadFile(dataUrl.format(path=f + '.bz2'), 'images/' + f)
            # need to make sure the image existed, easiest way to do this is
            # check the returned file size (some html would be kB, rather than a
            # MB-size fits file)
            subprocess.call(
                'du -h {} > {}'.format(loc + 'images/' + f, loc + '_tmpSzlog.log'),
                shell=True
            )
            with open(loc+'_tmpSzlog.log', 'r') as fl:
                size = fl.read().strip().split('\t')[0]
                if not 'M' in size.split('\t')[0]:
                    # this file is too small, delete it (not Megabyte size)
                    printWarning('File too small: <{}>'.format(f))
                    failedFiles.append(f)
                    os.remove(loc + 'images/' + f)

        # remove failed download from the list of fits
        for f in failedFiles: fitsFiles.remove(f)
        os.remove(loc + '_tmpSzlog.log')
        return fitsFiles
    elif len(returned[0]['Rows']) > 1:
        # we have a few possible galaxies returned
        # TODO: use GC separation to find one nearest to provided ra and dec
        #       rather than spamming the SDSS API
        printWarning('Multiple returned, decreasing threshold to {}'.format(threshold/10))
        return getFitsFromCoord(
            ra=ra,
            dec=dec,
            threshold=threshold / 5
        )
    else:
        # didn't find anything, increase the search radius and try again
        printWarning('Threshold too small, increasing')
        return getFitsFromCoord(
            ra=ra,
            dec=dec,
            threshold=threshold*3
        )

def cutFits(f, ra, dec, size=(200, 200), sigma=False):
    if not os.path.isfile(f): return None
    if type(size) == type(1): size =  (size, size)
    try:
        if len(size) != 2: return
    except Exception as e:
        print("\033[31msize must be an int or length-2 list/tuple/array\033[0m")
    fFile = fits.open(f)
    frame = fFile[0].header['SYSTEM'].strip()
    p = SkyCoord(ra = ra*u.degree, dec=dec*u.degree, frame=frame.lower())
    size = u.Quantity((50, 50), u.arcsec)
    wcs = WCS(header=fFile[0].header)
    cutout = Cutout2D(fFile[0].data, p, size, wcs=wcs)
    if sigma:
        sigma = getSigma(fFile)
        sigma_cutout = Cutout2D(sigma_data, p, size, wcs=wcs)
        return cutout.data, sigma_cutout.data
    else:
        return cutout.data

def getSigma(fFile):
    img = fFile[0].data
    simg = getSky(fFile)
    cimg = getCalib(fFile)
    dn = img / cimg + simg
    camcol = self.header["CAMCOL"]
    filter = self.header["FILTER"]
    run = self.header["RUN"]
    darkVariance = get_darkvariance(camcol, band, run)
    gain = get_gain(camcol, band, run)
    dn_err= np.sqrt(dn / gain + darkVariance)
    img_err = dn_err * cimg
    return img_err

# Adapted from https://github.com/MickaelRigault/astrobject/blob/master/astrobject/instruments/sdss.py
def getSky(fFile):
    from scipy.interpolate import RectBivariateSpline as RBS
    data = fFile[2].data
    sky = data['allsky'][0]
    xi = np.arange(sky.shape[0])
    yi = np.arange(sky.shape[1])
    interp = RBS(xi, yi, sky)
    return interp(data['yinterp'], data['xinterp'])

def getCalib(fFile):
    return fFile[1].data * np.ones(fFile[0].shape)

# Thanks to https://github.com/MickaelRigault/astrobject/blob/master/astrobject/instruments/sdss.py
def get_darkvariance(camcol,band,run=None):
    """
    data.sdss3.org/datamodel/files/BOSS_PHOTOOBJ/frames/RERUN/RUN/CAMCOL/frame.html
    """
    DARK_VAR_CCD = {
            0:{"u":9.61,   "g":15.6025,"r":1.8225,
               "i":7.84,     "z":0.81},
            1:{"u":12.6025,"g":1.44,   "r":1.00,
               "i":[5.76,6.25],"z":1.0},
            2:{"u":8.7025, "g":1.3225, "r":1.3225,
               "i":4.6225,   "z":1.0},
            3:{"u":12.6025,"g":1.96,   "r":1.3225,
               "i":[6.25,7.5625],"z":[9.61,12.6025]},
            4:{"u":9.3025, "g":1.1025, "r":0.81,
               "i":7.84,     "z":[1.8225,2.1025]},
            5:{"u":7.0225, "g":1.8225, "r":0.9025,
               "i":5.0625,   "z":1.21}
            }

    dark = DARK_VAR_CCD[camcol-1][band]
    # ----------
    # - output
    if type(dark) == float:
        return dark
    if run is None:
        raise ValueError("there is two dark-variance possibilites for "+\
                         " *camcol* %d, *band* %s "%(
                            camcol-1,band) + "Please, provide a *run*")

    return dark[1] if run>1500 else dark[0]

# Thanks to https://github.com/MickaelRigault/astrobject/blob/master/astrobject/instruments/sdss.py
def get_gain(camcol,band,run=None):
    """
    data.sdss3.org/datamodel/files/BOSS_PHOTOOBJ/frames/RERUN/RUN/CAMCOL/frame.html
    """
    GAIN_CCD = {
            0:{"u":1.62, "g":3.32, "r":4.71,
               "i":5.165,"z":4.745},
            1:{"u":[1.595,1.825],"g":3.855,"r":4.6,
               "i":6.565,"z":5.155},
            2:{"u":1.59,"g":3.845,"r":4.72,
               "i":4.86,"z":4.885},
            3:{"u":1.6,"g":3.995,"r":4.76,
               "i":4.885,"z":4.775},
            4:{"u":1.47,"g":4.05,"r":4.725,
               "i":4.64,"z":3.48},
            5:{"u":2.17,"g":4.035,"r":4.895,
               "i":4.76,"z":4.69}}

    gain = GAIN_CCD[camcol-1][band]
    # ----------
    # - output
    if type(gain) == float:
        return gain
    if run is None:
        raise ValueError("there is two gain possibilites for *camcol* %d, *band* %s "%(
            camcol-1,band) + "Please, provide a *run*")

    return gain[1] if run>1100 else gain[0]

def main():
    ra = 20.91722848729911
    dec = 0.6823345549867197
    x = [loc+'/images/'+i for i in getFitsFromCoord(ra=ra, dec=dec, threshold=3E-6)]
    rFile = [i for i in x if 'frame-r' in i.split('/')[-1]][0]
    gFile = [i for i in x if 'frame-g' in i.split('/')[-1]][0]
    bFile = [i for i in x if 'frame-u' in i.split('/')[-1]][0]
    ims = []
    fig, ax = plt.subplots(3)
    c = [
        cutFits(rFile, ra, dec),
        cutFits(gFile, ra, dec),
        cutFits(bFile, ra, dec),
    ]
    for i in range(3):
        ims += [ax[i].imshow(c[i], cmap='bone')]
        fig.colorbar(ims[-1], ax=ax[i])
    plt.show()

if __name__ == "__main__": main()
