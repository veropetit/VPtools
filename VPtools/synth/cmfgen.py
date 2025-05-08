# module cmfgen

from . import synth as synthLib
import pandas as pd
import numpy as np

galOB_path_name = 'ResearchGroupResources-Files/ReferenceMaterial/CMFGEN/CMFGEN-OB/{}_{}.spec'

def get_table_galOB():
    '''
    Function to return a pandas dataframe with the Galactic CMFGEN models that are available
    '''
    url = "https://drive.google.com/file/d/1-7Z1DfOgQbB2hYhyNEF1mEE3vg1p87Yd/view?usp=sharing"
    url='https://drive.google.com/uc?id=' + url.split('/')[-2]
    table = pd.read_csv(url)
    return table

def read(file):
    '''
    Function to read in a cmfgen synthetic spectrum into a 
    Synth object. The function returns two spectra: 
    the flux spectra and the spectra normalized to the continuum. 
    '''
    # The files are simply two columns
    wl, flux, norm = np.genfromtxt(file, unpack=True)

    return synthLib.Synth(wl, flux), synthLib.Synth(wl, norm)