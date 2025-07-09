# Module synth3
# Tools to manipulate synth3 models

import numpy as np
import pandas as pd
from . import synth as synthLib

synth3_path_name='ResearchGroupResources-Files/ReferenceMaterial/synth3/0vsiniSpec/T{0}G{1}'

def get_table_synth3():
    '''
    Function to return a pandas dataframe with the synth3 models that are available
    '''
    url = "https://drive.google.com/file/d/1vAIA1c4jWq8cPAyFoYEyZ2e75NcYv7WI/view?usp=sharing"
    url='https://drive.google.com/uc?id=' + url.split('/')[-2]
    table = pd.read_csv(url)
    return table

def read(codename):
    '''
    Function to read in a synth3 synthetic spectrum into a Synth object.

    :param codename: the path+modelcode (i.e T07000G40) of the model to open
    :param rtype: 2 Synth objects
    '''
    # The files are simply two columns
    normwl = np.genfromtxt('{}.s'.format(codename), skip_header=2)[:,0]*10
    normflux = np.genfromtxt('{}.s'.format(codename), skip_header=2)[:,1]
    conwl = np.genfromtxt('{}_contin.s'.format(codename), skip_header=2)[:,0]*10
    conflux = np.genfromtxt('{}_contin.s'.format(codename), skip_header=2)[:,1]

    return synthLib.Synth(normwl, normflux), synthLib.Synth(conwl, conflux)