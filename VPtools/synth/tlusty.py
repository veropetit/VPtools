# Module tlusty
# Tools to manipulate tlusty models from the OSTAR2002 and BSTARS2006 databases

import numpy as np
from . import synth as synthLib

def read(codename):
    '''
    Function to read in a Tlusty synthetic spectrum from the OSTAR2002 or BSTAR2006 grids. 
    Return the calibrated flux spectrum and continuum spectrum on their original wavelength grids.

    :param codename: the path+codename of the model to open (so without the '.7.gz' or '.17.gz')
    :param rtype: 2 Synth objects

    '''

    file_flux = '{}.7.gz'.format(codename)
    file_cont = '{}.17.gz'.format(codename)

    wl_flux, flux,  = np.genfromtxt(file_flux, unpack=True)
    wl_cont, cont = np.genfromtxt(file_cont, unpack=True)

    #cont_interp = np.interp(wl_flux, wl_cont, cont)

    return synthLib.Synth(wl_flux, flux), synthLib.Synth(wl_cont, cont)

def read_norm(codename):
    '''
    Function to read in a Tlusty synthetic spectrum and return the normalized spectrum
    in a Synth object.

    :param codename: the path+codename of the model to open (so without the '.7.gz' or '.17.gz')
    :rtype: Synth the normalized flux spectrum in a specpolFlow Spectrum object 

    '''

    flux, cont = read(codename)

    return synthLib.Synth(flux.wl, flux.specI/cont.interp(flux.wl).specI)