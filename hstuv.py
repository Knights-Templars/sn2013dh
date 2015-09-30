'''
Created on Jun 17, 2014

@author: cmccully
'''
from glob import glob
from pyraf import iraf
iraf.stsdas()
iraf.hst_calib()
iraf.stis()
import pyfits
from numpy import isnan, min, logical_and, bitwise_and, linspace, median, arange, correlate, exp, average
from scipy.interpolate import interp1d
from scipy.integrate import trapz
import os
from scipy.optimize import curve_fit

import stistools
iraf.twodspec()
iraf.apextract()
iraf.set(clobber='yes')

def gauss(x, a, x0, sigma, sky):
    return a * exp(-(x - x0) ** 2 / (2 * sigma ** 2) + sky)


def tofits(filename, data, hdr=None, clobber=False):
    """simple pyfits wrapper to make saving fits files easier."""
    from pyfits import PrimaryHDU, HDUList
    hdu = PrimaryHDU(data)
    if hdr is not None:
        hdu.header = hdr
    hdulist = HDUList([hdu])
    hdulist.writeto(filename, clobber=clobber, output_verify='ignore')

os.environ['oref'] = '/Users/cmccully/Documents/sne/sn2013dh/cti_corrected/ref/references/hst/oref/'
os.environ['ctirefs'] = '/Users/cmccully/Documents/sne/sn2013dh/cti_corrected/ref'


def combine(doreduce=True, doshifts=True):

    if doreduce:
        ims = glob('oc01020[1-4]*_raw.fits')
        ims += glob('oc016*_raw.fits')
        # for each image
        for im in ims:
            print im
            stistools.basic2d.basic2d(im, im.replace('raw','flt'))
            im = im.replace('raw', 'flt')
            print im
            # Update the target position at 0.0
            for i in range(4):
                pyfits.setval(im, 'POSTARG2', value=0.0, ext=i)
            # reset the aperture table to the newer file (we maybe should check this)
            pyfits.setval(im, 'APERTAB', value='oref$y2r1559to_apt.fits')
            # Reset the wcs to have CRPIX2 along the trace

            # Run x2d on the flt frame

            stistools.x2d.x2d(input=im, output=im.replace('flt','x2d') )

            h = pyfits.open(im.replace('flt','x2d'), mode='update')

            # Replace all of the bad pixels in the image by -10000 based on the DQ array
            # save them to a new file
            # Throw away bad reference file pixels and saturated pixels. None of the other error codes
            # were in the first file so I haven't included them here, but we might want to
            d = h[3].data
            badpix = logical_and(bitwise_and(d, 256) == 256,
                                 bitwise_and(d, 512) == 512)
            h[1].data[badpix] = -10000
            h.flush()

            # Trim the image
            for i in range(1,4):
                h[i].data = h[i].data[100:-100, 120:-100]
                h[i].header['CRPIX1'] -= 120
                h[i].header['CRPIX2'] -= 100
            h.close()

    ims = glob('oc01020[1-4]*_x2d.fits')
    ims += glob('oc01610[1-2]*_x2d.fits')
    if doshifts:
        init_guesses = [501, 542, 522, 523, 541, 524]
        centroids = []
        for i, im in enumerate(ims):
            print(im)
            h = pyfits.open(im)
            d = average(h[1].data[:, 915:925], axis=1)
            popt, _pcov = curve_fit(gauss, arange(len(d)), d, p0=[10, init_guesses[i], 1.5, 0])
            centroids.append(popt[1])
            shift = centroids[0] - popt[1]
            from matplotlib import pyplot as pl
            pl.ion()
            pl.clf()
            pl.plot(arange(len(d)), d)
            pl.plot(arange(len(d)), gauss(arange(len(d)), popt[0], popt[1], popt[2], popt[3]))
            _w = raw_input('Press return to continue')
            # watch the sign convention
            # This gives what you need to shift the input image by to get to the reference image
            iraf.unlearn(iraf.imshift)
            iraf.imshift(im + '[1]', im[:-8] + 'shift1.fits', 0.0, shift,
                         interp_type='drizzle')

    # Run imcombine on the rectified (but not flux scaled) images with crreject
    iraf.unlearn(iraf.imcombine)
    imlist = ''
    for im in ims:
        imlist += im[:-8] + 'shift1.fits,'
    # drop the last comma
    imlist = imlist[:-1]
    iraf.imcombine(input=imlist, output='13dh_uv_com.fits', reject='none',
                   lthreshold=-20, hthreshold=300)

    # run apall on the combined flux scaled image
    iraf.unlearn(iraf.apall)
    iraf.apall(input='13dh_uv_com.fits', output='13dh_uv', review='no',
               line=1024, nsum=-50, b_order=2, b_function='legendre',
               b_niterate=30, b_naverage=-21, nfind=1, t_order=2,
               background='fit', weights='variance')
