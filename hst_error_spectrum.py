__author__ = 'cmccully'
from astropy.io import fits
from scipy.signal import savgol_filter
import numpy as np

def fitshdr_to_wave(hdr):
    crval = float(hdr['CRVAL1'])
    crpix = float(hdr['CRPIX1'])
    # Convert crpix to be zero indexed
    crpix -= 1
    if 'CDELT1' in hdr.keys():
        cdelt = float(hdr['CDELT1'])
    else:
        cdelt = float(hdr['CD1_1'])
    npix = float(hdr['NAXIS1'])
    lam = np.arange(crval - cdelt * crpix ,
                    crval + cdelt * (npix - crpix) - 1e-4,
                    cdelt)
    return lam

hdu = fits.open('13dh_hst.fits')
lam = fitshdr_to_wave(hdu[0].header)
# Estimate the noise
s = savgol_filter(hdu[0].data, 101, 5)
noise = hdu[0].data - s
# Do a moving standard deviation
err = np.zeros(lam.shape)
for i in range(50,len(err) - 50):
    err[i] = np.std(noise[i-50:i+51])
err[:50] = err[50]
err[-50:] = err[-51]

d = np.zeros((3, len(err)))
d[0] = lam
d[1] = hdu[0].data
d[2] = err

m = np.median(d[1])
d[1] /= m
d[2] /= m

np.savetxt('13dh_hst.dat', d.transpose())