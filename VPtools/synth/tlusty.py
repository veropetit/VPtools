# Module tlusty
# Tools to manipulate tlusty models from the OSTAR2002 and BSTARS2006 databases

import numpy as np
#import specpolFlow as pol

def read_tlusty(codename):
    '''
    Function to read in a Tlusty synthetic spectrum from the OSTAR2002 or BSTAR2006 grids. 
    The returned continuum flux array is interpolated over the wavelength grid of the flux.

    :param codename: the path+codename of the model to open (so without the '.7.gz' or '.17.gz')
    :param rtype: 3 arrays with wavelength, flux, and continnum flux (interpolated)

    '''

    file_flux = '{}.7.gz'.format(codename)
    file_cont = '{}.17.gz'.format(codename)

    wl_flux, flux,  = np.genfromtxt(file_flux, unpack=True)
    wl_cont, cont = np.genfromtxt(file_cont, unpack=True)

    cont_interp = np.interp(wl_flux, wl_cont, cont)

    return wl_flux, flux, cont_interp

# def read_tlusty_norm(codename):
#     '''
#     Function to read in a Tlusty synthetic spectrum and return the normalized spectrum
#     in a specpolFlow Spectrum object

#     :param codename: the path+codename of the model to open (so without the '.7.gz' or '.17.gz')
#     :rtype: pol.Spectrum the normalized flux spectrum in a specpolFlow Spectrum object 

#     '''

#     wl, flux, cont = read_tlusty(codename)

#     Z = np.zeros(len(wl)) # arrays of zero to put in the specpolFLow object

#     return pol.Spectrum(wl, flux/cont, 
#                         Z, Z, Z, Z, header='codename')