# module synth
# Create an object to store and manipulate synthetic models
# Based on the specpolFlow Spectrum class

import numpy as np
import copy 
import specpolFlow as pol
import VPtools as vp
# Only for debugging
import matplotlib.pyplot as plt

const_c = 299792.458  #speed of light in km/s


def tlusty_binary(path_google,model1,model2,vsini1,vsini2,R2R1,vrad1,vrad2):
    '''
    Combines two spectra into a joint binary model spectrum

    :param path_google: the google path of shared drives
    :param model1: a TLUSTY path name for the 1st model
    :param model2: a TLUSTY path name for the 2nd model
    :param vsini1: the vsini of star 1
    :param vsini2: the vsini of star 2
    :param R2R1: the radius ratio of stars (R2/R1)
    :param vrad1: the radial velocity of star 1
    :param vrad2: the radial velocity of star 2

    :param rtype: Synth object
    '''
    synth1, synth1_con = vp.tlusty.read(path_google+model1) #read in the first star
    synth2, synth2_con = vp.tlusty.read(path_google+model2)
    return modelSB(synth1,synth1_con,synth2,synth2_con,vsini1,vsini2,R2R1,vrad1,vrad2)

def cmfgen_binary(path_google,model1,model2,vsini1,vsini2,R2R1,vrad1,vrad2):
    '''
    Combines two spectra into a joint binary model spectrum

    :param path_google: the google path of shared drives
    :param model1: a CMFGEN path name for the 1st model
    :param model2: a CMFGEN path name for the 2nd model
    :param vsini1: the vsini of star 1
    :param vsini2: the vsini of star 2
    :param R2R1: the radius ratio of stars (R2/R1)
    :param vrad1: the radial velocity of star 1
    :param vrad2: the radial velocity of star 2

    :param rtype: Synth object
    '''
    synth1, synth_norm1 = vp.cmfgen.read(path_google+model1)
    synth2, synth_norm2 = vp.cmfgen.read(path_google+model2)

    synth1_con=vp.Synth(synth1.wl,synth1.specI/synth_norm1.specI)
    synth2_con=vp.Synth(synth2.wl,synth2.specI/synth_norm2.specI)
    return modelSB(synth1,synth1_con,synth2,synth2_con,vsini1,vsini2,R2R1,vrad1,vrad2)

def synth3_binary(path_google,model1,model2,vsini1,vsini2,R2R1,vrad1,vrad2):
    '''
    Combines two spectra into a joint binary model spectrum

    :param path_google: the google path of shared drives
    :param model1: a synth3 path name for the 1st model
    :param model2: a synth3 path name for the 2nd model
    :param vsini1: the vsini of star 1
    :param vsini2: the vsini of star 2
    :param R2R1: the radius ratio of stars (R2/R1)
    :param vrad1: the radial velocity of star 1
    :param vrad2: the radial velocity of star 2

    :param rtype: Synth object
    '''
    synth_norm1, synth1_con = vp.synth3.read(path_google+model1)
    synth_norm2, synth2_con = vp.synth3.read(path_google+model2)

    synth1=vp.Synth(synth_norm1.wl,synth1_con.specI*synth_norm1.specI)
    synth2=vp.Synth(synth_norm2.wl,synth2_con.specI*synth_norm2.specI)
    return modelSB(synth1,synth1_con,synth2,synth2_con,vsini1,vsini2,R2R1,vrad1,vrad2)

def modelSB(synth1,synth1_con,synth2,synth2_con,vsini1,vsini2,R2R1,vrad1,vrad2):
    '''
    Combines two spectra into a joint binary model spectrum

    :param synth1: a synth object of the unnormalized spectra for star 1
    :param synth1_con: a synth object of the continuum for star 2
    :param synth2: a synth object of the unnormalized spectra for star 1
    :param synth2_con: a synth object of the continuum for star 2
    :param vsini1: the vsini of star 1
    :param vsini2: the vsini of star 2
    :param R2R1: the radius ratio of stars (R2/R1)
    :param vrad1: the radial velocity of star 1
    :param vrad2: the radial velocity of star 2

    :param rtype: Synth object
    '''

    #shift the spectra to the vrad of star1
    shifted1=synth1.doppler_shift(vrad1)
    shifted1c=synth1_con.doppler_shift(vrad1)

    wl=synth1.wl #defines the output wavelength range

    #interpolates the vrad shifted spectra to the output wavelength range
    interp1=shifted1.interp(wl)
    interp1_con=shifted1c.interp(wl)

    #rotationally broadens the spectra
    broad1=interp1.rot_int_cmj(vsini=vsini1)
    broad1_con=interp1_con.rot_int_cmj(vsini=vsini1)

    #repeat for star 2
    shifted2=synth2.doppler_shift(vrad2)
    shifted2c=synth2_con.doppler_shift(vrad2)

    interp2=shifted2.interp(wl)
    interp2_con=shifted2c.interp(wl)

    broad2=interp2.rot_int_cmj(vsini=vsini2)
    broad2_con=interp2_con.rot_int_cmj(vsini=vsini2)    

    #combines the two synth spectra
    joint=((broad1.specI)+broad2.specI*R2R1**2)/((broad1_con.specI)+broad2_con.specI*R2R1**2)
    
    return(vp.synth.Synth(wl,joint))

class Synth():
    '''
    Contains a synthetic spectra

    Contains arrays: 
    
    * wl - wavelengths
    * specI - spectrum
    '''

    def __init__(self, wl, specI):
        
        self.wl = wl
        self.specI=specI

    def __getitem__(self, key):
        """
        Returns a Spectrum object with only the values at the specified index(s)

        :param key: the index or slice being checked
        :rtype: Spectrum
        """
        wl_s = self.wl[key]
        specI_s = self.specI[key]
 
        slice_spec = Synth(wl_s, specI_s)
        return slice_spec

    def __setitem__(self, key, newval):
        """
        Sets all values of the Spectrum at the specified location equal
        to the input Spectrum's values.

        :param key: the index or slice being overwritten
        :param newval: Spectrum whose values are to replace the overwritten ones
        """
        if not(isinstance(newval, Synth)):
            raise TypeError()
        else:
            self.wl[key] = newval.wl
            self.specI[key] = newval.specI

    def __len__(self):
        return len(self.wl)

    def interp(self, wl):
        '''
        Class function to interpolate a Synth spectrum object to a different wavelength grid
        '''
        return Synth(wl, np.interp(wl, self.wl, self.specI))

    def doppler_shift(self, velocity):
        '''
        Doppler shift the spectrum according to an input radial velocity.

        :param velocity: the radial velocity in km/s
        :rtype: Spectrum
        '''
        spec = copy.deepcopy(self) #work on a copy (not self!)
        spec.wl = spec.wl + spec.wl*velocity/const_c
        return spec

    def extract_wl_range(self, down, up):
        '''
        Function to return the portion of the spectrum between two wavelengths
        '''
        ind0, ind1 = np.searchsorted(self.wl,
                                         [down, up])
        return self[ind0:ind1]

    def convolveRot(self, vsini, eps):
        '''
        Perform an exact rotation convolution on the model's datapoints
        '''
    
        # For debug..
        #fig, ax = plt.subplots(1,1)

        #Verify that the spectrum is in order of wavelength
        if np.all(self.wl[1:] <= self.wl[:-1]):
            print('Warning -- spectrum not in order!!')

        # Create an empty array to store the result
        conv_spec = copy.deepcopy(self)

        # Get the velocity array, with the first wl element for v=0

        for i in range(len(self)):

            # range of relevant wavelengths
            delta_lambda = vsini/const_c * self.wl[i]

            # get the portion of the spectrum relevant for convolution
            # for this pixel, i.e. within n*sigma of this pixel
            ind0, ind1 = np.searchsorted(self.wl,
                                         [self.wl[i] - delta_lambda,
                                          self.wl[i] + delta_lambda])
            
            # if there are no elements in the kernel, return the original flux value
            if ind1-ind0 > 1:

                #Pad to make sure to get the whole of the kernel
                if ind0 > 0:
                    ind0=ind0-1
                if ind1 < len(self):
                    ind1 = ind1 + 1
                specT = self[ind0:ind1]

                # Get the velocity array
                v = (specT.wl - self.wl[i]) / self.wl[i] * const_c
                # Compute the kernel
                G = np.zeros(len(v))
                G[1:-1] = 2*(1-eps) * (1-(v[1:-1]/vsini)**2)**0.5 + 0.5*np.pi*eps*(1-(v[1:-1]/vsini)**2)
                # explicitely set the padding regions to zero

                # This normalizes the kernel in velocity space
                # So when doing the trapeze integral below, I do it in velocity space
                # G = G / ( np.pi * vsini * (1-eps/3) )
                # This is the exact normalization. 
                # But if the wavelength grid is not fine enough for the vsini
                # then there will be some error from the numerical integration. 
                # So here we normalize with the numerical trapeze integral
                G = G / np.sum((G[:-1] + G[1:])*0.5*(v[1:] - v[:-1]))

                #ax.plot(v, G)

                # compute the convolution 
                prod = G*(specT.specI)
                # Make a trapeze integration
                conv_spec.specI[i] = np.sum((prod[:-1] + prod[1:])*0.5*(v[1:] - v[:-1]))

        return conv_spec
    

    def rot_int_cmj(self, vsini, eps=0.6, nr=10, ntheta=100, dif = 0.0):
        '''
        From RNAAS: Carvalho & Johns-Krull (2023)
        A routine to quickly rotationally broaden a spectrum in linear time.

        VP: modified to work on a Synth object

        INPUTS:
        self - a Synth object
        
        vsini (km/s) - projected rotational velocity
        
        OUTPUT:
        ns - a rotationally broadened Synth object on the wavelength scale w

        OPTIONAL INPUTS:
        eps (default = 0.6) - the coefficient of the limb darkening law
        
        nr (default = 10) - the number of radial bins on the projected disk
        
        ntheta (default = 100) - the number of azimuthal bins in the largest radial annulus
                                note: the number of bins at each r is int(r*ntheta) where r < 1
        
        dif (default = 0) - the differential rotation coefficient, applied according to the law
        Omeg(th)/Omeg(eq) = (1 - dif/2 - (dif/2) cos(2 th)). Dif = .675 nicely reproduces the law 
        proposed by Smith, 1994, A&A, Vol. 287, p. 523-534, to unify WTTS and CTTS. Dif = .23 is 
        similar to observed solar differential rotation. Note: the th in the above expression is 
        the stellar co-latitude, not the same as the integration variable used below. This is a 
        disk integration routine.

        '''
        w = self.wl
        s = self.specI
        ns = self.specI*0.0

        tarea = 0.0
        dr = 1./nr
        for j in range(0, nr):
            r = dr/2.0 + j*dr
            area = ((r + dr/2.0)**2 - (r - dr/2.0)**2)/int(ntheta*r) * (1.0 - eps + eps*np.cos(np.arcsin(r)))
            for k in range(0,int(ntheta*r)):
                th = np.pi/int(ntheta*r) + k * 2.0*np.pi/int(ntheta*r)
                if dif != 0:
                    vl = vsini * r * np.sin(th) * (1.0 - dif/2.0 - dif/2.0*np.cos(2.0*np.arccos(r*np.cos(th))))
                    ns += area * np.interp(w + w*vl/2.9979e5, w, s)
                    tarea += area
                else:
                    vl = r * vsini * np.sin(th)
                    ns += area * np.interp(w + w*vl/2.9979e5, w, s)
                    tarea += area
            
        return Synth(w, ns/tarea)

    def to_spf(self):
        '''
        Returns a synth object into a specpolFlow Spectrum object (with error bars and Stokes V, N1, N2 set to zero)

        This allows for the use of the specpolFlow Spectrum functionalities (like performing LSD, line analysis, etc)
        '''
        # Create an array of zeros of the right length for V, N1, N2, errors
        Z = np.zeros(len(self))
        return(pol.Spectrum(self.wl, self.specI, Z, Z, Z, Z ))