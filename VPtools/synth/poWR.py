# module poWR

from . import synth as synthLib
import pandas as pd
import numpy as np


gal_path_name = 'ResearchGroupResources-Files/ReferenceMaterial/poWR/griddl-gal-ob-vd3-line_{}/gal-ob-vd3_{}_line.txt'

def get_table_gal():
    '''
    Function to return a pandas dataframe with the Galactic poWR models that are available
    '''
    url = "https://drive.google.com/file/d/18i8-6ijtejHWp8k75sx1_TLAiGdp09XF/view?usp=sharing"
    url='https://drive.google.com/uc?id=' + url.split('/')[-2]
    table = pd.read_csv(url)
    return table


def read(file):
    '''
    Function to read in a poWR synthetic spectrum into a 
    Synth object
    '''
    # The files are simply two columns
    wl, flux = np.genfromtxt(file, unpack=True)

    return synthLib.Synth(wl, flux)