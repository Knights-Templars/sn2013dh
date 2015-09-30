__author__ = 'cmccully'

from astropy.io import ascii
import numpy as np
from matplotlib import pyplot

def read_data(filename):
    names = ['filter', 'mjd', 'mag', 'magerr', 'maglim3sig',
             'satlim0.98', 'countrate', 'rateerr', 'aperture', 'frametime', 'exptime',
             'telapse']

    data = ascii.read(filename, names=names)
    for column in data.columns[1:]:
        data[column][data[column] =='NULL'] = -99.9
        print(column)

        for i in data[column]: print(i)
        data[column].dtype = np.float
    return data


def plot():
    dh_data = read_data('SN2013dh_uvotB14.1.dat')
    fe_data = read_data('SN2011fe_uvotB14.1.dat')

    # Calculate the UVW2 - UVW1 colors
    fe_uvw1 = fe_data['filter'] == 'UVW1'
    fe_uvw2 = fe_data['filter'] == 'UVW2'

    fe_uv_mjds = []
    fe_uv_colors = []
    for i, mjd in enumerate(fe_data[fe_uvw2]['mjd']):
        datediff = np.abs(mjd - fe_data[fe_uvw1])
        if np.min(datediff) < 1:
            fe_uvw1_mag = fe_data[fe_uvw1][np.argmin(datediff)]
            fe_uv_colors.append(fe_data[fe_uvw2]['mag'][i] - fe_uvw1_mag)
            fe_uv_mjds.append(mjd)

    pyplot.plot(fe_uv_mjds, fe_uv_colors)
    pyplot.ylim(-2, 2)
    pyplot.show()