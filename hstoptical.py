'''
Created on Jun 17, 2014

@author: cmccully
'''
from glob import glob

import numpy as np
import pyfits
from astropy.io import fits
from numpy import isnan,min,logical_and,bitwise_and,linspace,median,arange,correlate, average, exp
from scipy.interpolate import interp1d
from scipy.integrate import trapz
from scipy.optimize import curve_fit
import os
os.environ['oref'] = '/Users/cmccully/Documents/sne/sn2013dh/cti_corrected/ref/references/hst/oref/'
os.environ['ctirefs'] = '/Users/cmccully/Documents/sne/sn2013dh/cti_corrected/ref'

from astroscrappy import detect_cosmics

from pyraf import iraf
iraf.stsdas()
iraf.hst_calib()
iraf.stis()
iraf.mstools()
iraf.fitting()
iraf.longslit()

import stistools
iraf.twodspec()
iraf.apextract()
iraf.onedspec()
iraf.set(ctirefs='/Users/cmccully/Documents/sne/sn2013dh/cti_corrected/ref/')
iraf.set(clobber='YES')


def gauss(x, a, x0, sigma, sky):
    return a * exp(-(x - x0) ** 2 / (2 * sigma ** 2) + sky)

def tofits(filename, data, hdr=None,clobber=False):
    """simple pyfits wrapper to make saving fits files easier."""
    from pyfits import PrimaryHDU,HDUList
    hdu = PrimaryHDU(data)
    if not hdr is None: hdu.header = hdr
    hdulist = HDUList([hdu])
    hdulist.writeto(filename, clobber=clobber,output_verify='ignore')

def fitshdr_to_wave(hdr):
    crval = float(hdr['CRVAL1'])
    if 'CDELT1' in hdr.keys():
        cdelt = float(hdr['CDELT1'])
    else:
        cdelt = float(hdr['CD1_1'])
    crpix = float(hdr['CRPIX1'])
    nlam = float(hdr['NAXIS1'])
    lam = np.arange(crval - cdelt * (crpix - 1), crval - cdelt * (crpix - 1) + cdelt * nlam - 1e-4, cdelt)
    return lam

def stitch():
    #Find the mean 3100 and 3125 of blue and uv
    uvhdu = pyfits.open('13dh_uv.fits')
    uvlam = fitshdr_to_wave(uvhdu[0].header)
    w = np.logical_and(uvlam > 3050, uvlam < 3100)
    uvmean = np.median(uvhdu[0].data[1,0,w])
    uvhdu.close()
    print("uv mean %e"%uvmean)


    bluehdu = pyfits.open('13dh_blue.fits')
    bluelam = fitshdr_to_wave(bluehdu[0].header)
    w = np.logical_and(bluelam > 3050, bluelam < 3100)
    blueuvmean = np.median(bluehdu[0].data[1,0,w])
    print("blue mean uv %e"%blueuvmean)

    # Find mean of red and blue between 5375 and 5600
    w = np.logical_and(bluelam > 5375, bluelam < 5600)
    blueredmean = bluehdu[0].data[1,0,w].mean()
    bluehdu.close()
    print("blue mean red %e"%blueredmean)

    redhdu = pyfits.open('13dh_red.fits')
    redlam = fitshdr_to_wave(redhdu[0].header)
    w = np.logical_and(redlam > 5375, redlam < 5600)
    redmean = redhdu[0].data[1,0,w].mean()
    redhdu.close()
    print("red mean %e"%redmean)

    # trim uv at 3140
    iraf.unlearn(iraf.scopy)
    iraf.scopy('13dh_uv.fits', '13dh_uv_trim.fits', w1='INDEF', w2=3140, rebin='no')
    
    # trim blue at 3130 and 5600
    iraf.unlearn(iraf.scopy)
    iraf.scopy('13dh_blue.fits', '13dh_blue_trim.fits', w1=3130, w2=5600, rebin='no')
    
    # Trim red at 5375
    iraf.unlearn(iraf.scopy)
    iraf.scopy('13dh_red.fits', '13dh_red_trim.fits', w1=5375, w2='INDEF', rebin='no')

    # Copy the spectra from the second extraction to the first
    for im in ['13dh_uv_trim.fits', '13dh_blue_trim.fits', '13dh_red_trim.fits']:
        hdu = pyfits.open(im, mode='update')
        hdu[0].data[0, 0, :] = hdu[0].data[1, 0, :]
        hdu.flush()
        hdu.close()

    # Write out the scale factors
    lines = ['%f\n' % (blueuvmean / uvmean)**-1,'1.0\n','%f\n' % (blueredmean / redmean)**-1]

    f = open('scales.dat','w')
    f.writelines(lines)
    f.close()
    #Scombine with the correct scalings using average
    iraf.unlearn(iraf.scombine)
    iraf.scombine('13dh_uv_trim, 13dh_blue_trim, 13dh_red_trim', '13dh_hst.fits', scale='@scales.dat')
    return

def combine(do_cti=False, doreduce=True, doshifts=True):

    if do_cti:
        os.system('stis_cti --crds_update')
    if doreduce:
        # Defringing didn't seem to converge because of the low S/N
        stistools.ocrreject.ocrreject('oc0102070_flc.fits','oc0102070_crc.fits')
        iraf.normspflat(inflat='oc0102070_crc.fits',outflat='oc0102070_nsp.fits', do_cal='no')

        iraf.imcalc(input='oc0102070_nsp.fits', output='temp_nsp.fits', equals='if(x .lt. 250) then 1 else im1')
        iraf.imcopy('temp_nsp.fits[1][1:250,*]', 'oc0102070_nsp.fits[1][1:250,*]')

        #iraf.defringe('oc0102060_flc.fits', 'oc0102070_nsp.fits', 'oc0102060_dfr.fits')
        #for each image
        for im in ['oc0102050_flc','oc0102060_flc']:
            outbase = 'blue'
            if im[:-4] == 'oc0102060':
                outbase = 'red'
            #reset the aperture table to the newer file (we maybe should check this)
            pyfits.setval(im +'.fits','APERTAB',value='oref$y2r1559to_apt.fits')
            pyfits.setval(im +'.fits', 'SPTRCTAB', value='oref$qa31608go_1dt.fits')

            # fixpix any negative values. In principle some of this noise
            # could be real, but I have found that is often not the case
            hdu = fits.open(im+ '.fits')
            mask1 = hdu[1].data < -20
            mask2 = hdu[4].data < -20
            hdu.close()
            fits.writeto(outbase+'mask1.fits', mask1.astype('i'), clobber=True)
            fits.writeto(outbase+'mask2.fits', mask2.astype('i'), clobber=True)

            iraf.unlearn(iraf.fixpix)
            iraf.fixpix(im+'[1]', outbase+'mask1.fits')

            iraf.unlearn(iraf.fixpix)
            iraf.fixpix(im+'[4]', outbase+'mask2.fits')

            # Subtract off the median value
            hdu = fits.open(im+ '.fits', mode='update')
            hdu[1].data -= np.median(hdu[1].data)
            hdu[4].data -= np.median(hdu[4].data)

            readnoise1 = 1.4826 * np.median(np.abs(hdu[1].data))
            readnoise2 = 1.4826 * np.median(np.abs(hdu[4].data))

            # Cosmic ray reject both images using scrappy
            # Make sure to treat the noise in a sensible way
            crmask1, clean1 = detect_cosmics(hdu[1].data, readnoise=readnoise1,
                                             sigclip=5, objlim=5, sigfrac=0.8,
                                             fsmode='median', psfmodel='gaussy',
                                             psffwhm=2., cleantype='idw')

            crmask2, clean2 = detect_cosmics(hdu[4].data, readnoise=readnoise2,
                                             sigclip=5, objlim=5, sigfrac=0.8,
                                             fsmode='median', psfmodel='gaussy',
                                             psffwhm=2., cleantype='idw')

            hdu.flush()
            hdu.close()

            fits.writeto(outbase + '_crmask1.fits', crmask1.astype('i'), clobber=True)
            fits.writeto(outbase + '_crmask2.fits', crmask2.astype('i'), clobber=True)
            # Run fixpix on the frames
            iraf.unlearn(iraf.fixpix)
            iraf.fixpix(im+'[1]', outbase+'_crmask1.fits')

            iraf.unlearn(iraf.fixpix)
            iraf.fixpix(im+'[4]', outbase+'_crmask2.fits')

            if outbase=='red':
                iraf.mkfringeflat('oc0102060_flc.fits', 'oc0102070_nsp.fits', 'oc0102070_frr.fits',
                                  beg_scale=0.6, end_scale=1.5, scale_step=0.01,
                                  beg_shift=-3.0, end_shift=3.0,shift_step=0.05)
                iraf.defringe('oc0102060_flc.fits', 'oc0102070_frr.fits', 'oc0102060_dfr.fits')
                #Run x2d on the flt frame
                stistools.x2d.x2d(input='oc0102060_dfr.fits',output=im[:-4]+'x2d.fits')
            else:
                stistools.x2d.x2d(input='oc0102050_flc.fits',output=im[:-4]+'x2d.fits')

            h = pyfits.open(im[:-4]+'x2d.fits', mode='update')
            #Replace all of the bad pixels in the image by -666 based on the DQ array
            #save them to a new file
            #Throw away bad reference file pixels and saturated pixels. None of the other error codes 
            #were in the first file so I haven't included them here, but we might want to
            d = h[3].data
            badpix = logical_and(bitwise_and(d,256) == 256,bitwise_and(d,512) == 512)
            h[1].data[badpix] = -10000
            d = h[6].data
            badpix = logical_and(bitwise_and(d,256) == 256,bitwise_and(d,512) == 512)
            h[4].data[badpix] = -10000
            h.flush()

            # Trim the images
            for i in range(1,7):
                h[i].data = h[i].data[100:-100, 100:-100]
                h[i].header['CRPIX1'] -= 100
                h[i].header['CRPIX2'] -= 100

            h.flush()
            h.close()

            # Combine the images
            iraf.unlearn(iraf.imcombine)
            iraf.imcombine(input=im[:-4]+'x2d[1],'+im[:-4]+'x2d[4]', output=outbase+'_com.fits',
                            reject='crreject')

            hdu = pyfits.open(outbase +'_com.fits')
            mask = hdu[0].data == 0.0
            hdu.close()
            fits.writeto(outbase+'_mask.fits', mask.astype('i'), clobber=True)

            iraf.unlearn(iraf.fixpix)
            iraf.fixpix(outbase+'_com.fits', outbase+'_mask.fits')

            iraf.unlearn(iraf.apall)
            iraf.apall(input=outbase+'_com',output='13dh_'+outbase, review='no',
                       nsum = -500, b_order = 1,
                       b_function='legendre',b_niterate=30, b_naverage = -21,
                       nfind=1,t_order=3,background='median',weights='variance',
                       skybox=100 )
            iraf.splot(outbase+'[SCI]')


