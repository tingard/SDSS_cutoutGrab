# TODO: credit Hugh?
import os, json, bz2, re, subprocess
import requests
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord, Angle
from astropy.nddata.utils import NoOverlapError
import astropy.units as u
from astropy.nddata import Cutout2D
from astropy.wcs import WCS

def printWarning(s): print('\033[33m[sdssCutoutGrab Warning] ' + s + '\033[0m')

def printInfo(s): print('🤖 \033[32m[sdssCutoutGrab Info] ' + s + '\033[0m')

def downloadFile(url, fname, overwrite=False, decompress=True):
    try:
        # check if we already have the file in the output file location
        if os.path.isfile(fname) and not overwrite:
            return True
        # check if directory tree needs updating
        if not os.path.isdir(os.path.abspath(fname)):
            fstruc = fname.split('/')
            for i in range(1, len(fstruc)):
                if not os.path.isdir('/'.join(fstruc[:i])):
                    printInfo('Making directory {}'.format('/'.join(fstruc[:i])))
                    os.mkdir('/'.join(fstruc[:i]))
        # start the download
        r = requests.get(url, stream=True)
        if decompress:
            decompressor = bz2.BZ2Decompressor()
        printInfo('writing to {}'.format(fname))
        with open(fname, 'wb') as imgf:
            for chunk in r:
                if decompress:
                    dChunk = decompressor.decompress(chunk)
                else:
                    dChunk = chunk
                imgf.write(dChunk)
        return True
    except FileNotFoundError as fnfe:
        print(url, fname)
        raise(fnfe)

def queryFromRaDec(ra, dec, radius=0.5, limit=10):
    queryUrl = 'http://skyserver.sdss.org/dr13/en/tools/search/x_results.aspx'
    res = requests.get(queryUrl, params={
        'searchtool': 'Radial',
        'TaskName': 'Skyserver.Search.Radial',
        'whichphotometry': 'optical',
        'coordtype': 'equatorial',
        'ra': ra,
        'dec': dec,
        'radius': radius,
        'limit': limit,
        'format': 'json',
    })
    if res.status_code == 200:
        try:
            result = res.json()[0]['Rows']
            return result
        except json.decoder.JSONDecodeError:
            print(res.url)
            printWarning('Could not parse returned JSON: ' + res.url)
            return []
    else:
        printWarning('Could not connect to SDSS skyserver: ' + res.url)
        return []

def getBandFits(dataObj, band='r', outFile=None, overwrite=False):
    # make sure we have correct parameters in provided object
    if not all(dataObj.get(i, False) for i in ['run', 'camcol', 'field']):
        return
    # define strings and trigger file download
    dataUrl = "http://data.sdss.org/sas/dr13/eboss/photoObj/frames/301/{path}.bz2"
    queryParams = '{run}/{camcol}/frame-{band}-{run:06d}-{camcol}-{field:04d}.fits'
    sdssFilePath = dataUrl.format(path=queryParams.format(band=band, **dataObj))
    fileName = outFile if outFile else 'images/' + queryParams.format(band=band, **dataObj)
    return fileName if downloadFile(sdssFilePath, fileName, overwrite=overwrite) else False

def cutFits(f, ra, dec, size=(4 * u.arcsec, 4 * u.arcsec)):
    if not os.path.isfile(f): return None
    try:
        if len(size) != 2: return
    except Exception as e:
        print("\033[31msize must be an int or length-2 list/tuple/array\033[0m")
    fFile = fits.open(f)
    frame = fFile[0].header['SYSTEM'].strip()
    p = SkyCoord(ra = ra*u.degree, dec=dec*u.degree, frame=frame.lower())
    sizeQuant = u.Quantity(size, u.arcsec)
    wcs = WCS(header=fFile[0].header)
    try:
        cutout = Cutout2D(fFile[0].data, p, sizeQuant, wcs=wcs)
        return cutout.data
    except NoOverlapError:
        printWarning(
            'Ra, Dec not inside frame for Ra: {ra}, Dec: {dec}, and frame: {f}'
                .format(ra=ra, dec=dec, f=f
            )
        )
    return False

def montageFolder(imageFolder, band, objID, montageFolder=''):
    imF = imageFolder[:-1] if imageFolder[-1] == '/' else imageFolder
    monF = montageFolder + '/' if len(montageFolder) and montageFolder[-1] != '/' else montageFolder
    print(imF, monF)
    if not os.path.exists(folder):
        raise OSError('Folder provided to montageFolder does not exist')
    k = { 'imF': imF, 'monF': monF, 'band': band, 'objID': objID }
    commands = '\n'.join(map(lambda i: i.format(k) for i in (
        'cd {imF}',
        '{monF}mImgtbl ./ raw.tbl',
        '{monF}mMakeHdr raw.tbl template.hdr',
        '{monF}mProjExec -p ./ raw.tbl template.hdr proj stats.tbl',
        '{monF}mImgtbl proj images.tbl',
        '{monF}mOverlaps images.tbl diffs.tbl',
        '{monF}mDiffExec -p proj diffs.tbl template.hdr diff',
        '{monF}mFitExec diffs.tbl fits.tbl diff',
        '{monF}mBgModel images.tbl fits.tbl corr.tbl',
        '{monF}mBgExec -p proj images.tbl corr.tbl corr',
        '{monF}mAdd -p corr images.tbl template.hdr ../{band}{objID}.fits',
    )))
    print(commands)

if __name__ == "__main__":
    r = queryFromRaDec(20.85483159, 0.7553174294, radius=0.2)[0]
    print(r)
    fname = getBandFits(r, overwrite=True)
    r = cutFits(fname, 20.85487832, 0.7552652609);
    plt.imshow(r, cmap='gray')
    plt.axis('off')
    plt.show();
