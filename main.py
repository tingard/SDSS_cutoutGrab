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

def cutFits(f, ra, dec, size):
    if not os.path.isfile(f): return None
    fFile = fits.open(f)
    frame = fFile[0].header['SYSTEM'].strip()
    p = SkyCoord(ra = ra*u.degree, dec=dec*u.degree, frame=frame.lower())
    size = u.Quantity((50, 50), u.arcsec)
    wcs = WCS(header=fFile[0].header)
    cutout = Cutout2D(fFile[0].data, p, (200, 200), wcs=wcs)
    return cutout.data

def main():
    ra = 20.91722848729911
    dec = 0.6823345549867197
    x = [loc+'/images/'+i for i in getFitsFromCoord(ra=ra, dec=dec, threshold=3E-6)]
    rFile = [i for i in x if 'frame-r' in i.split('/')[-1]][0]
    gFile = [i for i in x if 'frame-g' in i.split('/')[-1]][0]
    bFile = [i for i in x if 'frame-u' in i.split('/')[-1]][0]
    fig, (ax0, ax1, ax2) = plt.subplots(3)
    fig, ax = plt.subplots(3)
    r = cutFits(rFile, ra, dec, 10)
    im0 = ax0.imshow(r, cmap='bone')
    fig.colorbar(im0, ax=ax0)
    g = cutFits(gFile, ra, dec, 10)
    im1 = ax1.imshow(g, cmap='bone')
    fig.colorbar(im1, ax=ax1)
    b = cutFits(bFile, ra, dec, 10)
    im2 = ax2.imshow(b, cmap='bone')
    fig.colorbar(im2, ax=ax2)
    plt.show()
