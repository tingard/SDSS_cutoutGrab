# TODO: credit Hugh / edd edmonson
#       (https://github.com/eddedmondson/sdss_fits_cutout)
import os
import json
import bz2
import requests
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.nddata.utils import NoOverlapError
import astropy.units as u
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
import sdssCutoutGrab.sdss_psf as sdss_psf


def printWarning(s):
    print('\033[33m[sdssCutoutGrab Warning] ' + s + '\033[0m')


def printInfo(s):
    print('🤖 \033[32m[sdssCutoutGrab Info] ' + s + '\033[0m')


def downloadFile(url, fname, overwrite=False, decompress=True, verbose=False):
    try:
        # check if we already have the file in the output file location
        if os.path.isfile(fname) and not overwrite:
            return True
        # check if directory tree needs updating
        if not os.path.isdir(os.path.abspath(fname)):
            fstruc = fname.split('/')
            for i in range(1, len(fstruc)):
                if not os.path.isdir('/'.join(fstruc[:i])):
                    if verbose:
                        printInfo(
                            'Making directory {}'.format('/'.join(fstruc[:i]))
                        )
                    os.mkdir('/'.join(fstruc[:i]))
        # start the download
        r = requests.get(url, stream=True)
        if decompress:
            decompressor = bz2.BZ2Decompressor()
        if verbose:
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
        if verbose:
            print(url, fname)
        raise(fnfe)


def queryFromRaDec(ra, dec, radius=0.5, limit=10, verbose=False):
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
            resultJson = res.json()
            result = resultJson[0]['Rows']
            return sorted(
                result,
                key=lambda i: (
                    (i['ra'] - 160.65883)**2 + (i['dec'] - 23.95189)**2
                )
            )
            return result
        except json.decoder.JSONDecodeError:
            if versose:
                print(res.url)
                printWarning('Could not parse returned JSON: ' + res.url)
            return []
    else:
        if verbose:
            printWarning('Could not connect to SDSS skyserver: ' + res.url)
        return []


def getBandFits(frame, band='r', outFile=None, overwrite=False):
    # make sure we have correct parameters in provided object
    if not all(frame.get(i, False) for i in ['run', 'camcol', 'field']):
        return
    # define strings and trigger file download
    dataUrl = (
        'http://data.sdss.org/sas/dr13/eboss/photoObj/frames/301/'
        '{path}.bz2'
    )
    queryParams = (
        '{run}/{camcol}/'
        'frame-{band}-{run:06d}-{camcol}-{field:04d}.fits'
    )
    sdssFilePath = dataUrl.format(
        path=queryParams.format(band=band, **frame)
    )
    fileName = outFile if outFile else (
        'fits_images/' + queryParams.format(band=band, **frame)
    )
    return (
        fileName if downloadFile(sdssFilePath, fileName, overwrite=overwrite)
        else False
    )


def cutFits(fFile, ra, dec, size=(4, 4), sigma=False, verbose=False):
    if len(size) != 2:
        raise IndexError(
            "\033[31msize must be an int or length-2 list/tuple/array\033[0m"
        )
    try:
        frame = fFile[0].header['SYSTEM'].strip()
    except KeyError:
        frame = 'FK5'
    p = SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame=frame.lower())
    sizeQuant = u.Quantity(size, u.arcsec)
    wcs = WCS(header=fFile[0].header)
    try:
        cutout = Cutout2D(fFile[0].data, p, sizeQuant, wcs=wcs)
        if sigma:
            sigma = getSigma(fFile)
            sigma_cutout = Cutout2D(sigma, p, sizeQuant, wcs=wcs)
            return cutout, sigma_cutout
        else:
            return cutout
    except NoOverlapError:
        if verbose:
            print(
                'Ra, Dec not inside frame for Ra:'
                + ' {ra}, Dec: {dec}, and frame: {f}'.format(
                    ra=ra, dec=dec, f=fFile[0]
                )
            )
    return False


def getSigma(fFile):
    img = fFile[0].data
    hdr = fFile[0].header
    simg = getSky(fFile)
    cimg = getCalib(fFile)
    dn = img / cimg + simg
    camcol = hdr["CAMCOL"]
    band = hdr["FILTER"]
    run = hdr["RUN"]
    darkVariance = get_darkvariance(camcol, band, run)
    gain = get_gain(camcol, band, run)
    dn_err = np.sqrt(dn / gain + darkVariance)
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
def get_darkvariance(camcol, band, run=None):
    """
    data.sdss3.org/datamodel/files/BOSS_PHOTOOBJ/frames/RERUN/RUN/CAMCOL/frame.html
    """
    DARK_VAR_CCD = {
        0: {"u": 9.61,    "g": 15.6025, "r": 1.8225, "i": 7.84,
            "z": 0.81},
        1: {"u": 12.6025, "g": 1.44,    "r": 1.00,   "i": [5.76, 6.25],
            "z": 1.0},
        2: {"u": 8.7025,  "g": 1.3225,  "r": 1.3225, "i": 4.6225,
            "z": 1.0},
        3: {"u": 12.6025, "g": 1.96,    "r": 1.3225, "i": [6.25, 7.5625],
            "z": [9.61, 12.6025]},
        4: {"u": 9.3025,  "g": 1.1025,  "r": 0.81,   "i": 7.84,
            "z": [1.8225, 2.1025]},
        5: {"u": 7.0225,  "g": 1.8225,  "r": 0.9025, "i": 5.0625,
            "z": 1.21}
    }

    dark = DARK_VAR_CCD[camcol - 1][band]
    # ----------
    # - output
    if type(dark) == float:
        return dark
    if run is None:
        raise ValueError(
            'there is two dark-variance possibilites for '
            + (' *camcol* %d, *band* %s ' % (camcol - 1, band))
            + 'Please, provide a *run*'
        )

    return dark[1] if run > 1500 else dark[0]


def getPSF(galCoord, frame, fitsFile,
           fname='./tmpPsfFile.fit', deleteOnComplete=True):
    wcs = WCS(fitsFile[0].header)
    coords = wcs.wcs_world2pix([galCoord], 1)
    psfQueryUrl = 'https://data.sdss.org/sas/dr14/eboss/photo/redux/' + \
        '{rerun}/{run}/objcs/{camcol}/' + \
        'psField-{run:06d}-{camcol}-{field:04d}.fit'
    downloadFile(
        psfQueryUrl.format(**frame),
        fname,
        overwrite=True,
        decompress=False
    )
    psfield = fits.open(fname)
    bandnum = 'ugriz'.index('r')
    hdu = psfield[bandnum + 1]
    if deleteOnComplete:
        os.remove(fname)
    return sdss_psf.sdss_psf_at_points(hdu, *coords[0])


# Thanks to https://github.com/MickaelRigault/astrobject/blob/master/astrobject/instruments/sdss.py
def get_gain(camcol, band, run=None):
    """
    data.sdss3.org/datamodel/files/BOSS_PHOTOOBJ/frames/RERUN/RUN/CAMCOL/frame.html
    """
    GAIN_CCD = {
        0: {"u": 1.62, "g": 3.32, "r": 4.71,
            "i": 5.165, "z": 4.745},
        1: {"u": [1.595, 1.825], "g": 3.855, "r": 4.6,
            "i": 6.565, "z": 5.155},
        2: {"u": 1.59, "g": 3.845, "r": 4.72,
            "i": 4.86, "z": 4.885},
        3: {"u": 1.6, "g": 3.995, "r": 4.76,
            "i": 4.885, "z": 4.775},
        4: {"u": 1.47, "g": 4.05, "r": 4.725,
            "i": 4.64, "z": 3.48},
        5: {"u": 2.17, "g": 4.035, "r": 4.895,
            "i": 4.76, "z": 4.69}}

    gain = GAIN_CCD[camcol - 1][band]
    # ----------
    # - output
    if type(gain) == float:
        return gain
    if run is None:
        raise ValueError(
            "there is two gain possibilites for *camcol* %d, *band* %s " % (
                camcol - 1, band
            ) + "Please, provide a *run*"
        )

    return gain[1] if run > 1100 else gain[0]


def montageFolderCommand(imageFolder, band, objID, montageFolder=''):
    imF = imageFolder[:-1] if imageFolder[-1] == '/' else imageFolder
    monF = (
        montageFolder + '/'
        if len(montageFolder) and montageFolder[-1] != '/'
        else montageFolder
    )
    if not os.path.exists(imageFolder):
        raise OSError('Folder provided to montageFolder does not exist')
    k = {'imF': imF, 'monF': monF, 'band': band, 'objID': objID}
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
    return commands


if __name__ == "__main__":
    r = queryFromRaDec(20.85483159, 0.7553174294, radius=0.2)[0]
    print(r)
    fname = getBandFits(r, overwrite=True)
    ffile = fits.open(fname)
    r, sr = cutFits(ffile, 20.85487832, 0.7552652609, sigma=True)
    fig, ax = plt.subplots(1, 2)
    ax[0].imshow(r, cmap='gray')
    ax[1].imshow(sr, cmap='gray', vmin=0, vmax=r.max() / 10.0)
    ax[0].axis('off')
    ax[1].axis('off')
    plt.tight_layout()
    plt.show()
