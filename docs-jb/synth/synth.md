# Using synthetic spectra

In this section, we cover tools to manipulate and generate synthetic spectra. 

* python class `synth` that stores a synthetic spectrum and provides class functions for manipulation (such as rotational covulution and re-gridding)
* functions to read in synthetic spectra from various pre-computed grids
* function to display table and obtains paths/filenames of synthetic spectra hosted on our Shared Google Drive (requires access)

## What are synthetic spectra?

Synthetic spectra are computed by first solving the equations of stellar structure that have been tailored for the atmospheres of a star (see notes from PHYS 633). This provides at `atmosphere model` that describes the physical parameters at each layers of the atmosphere. Then a more detailed calculation of the emerging spectrum (so the flux at $\tau=0$) is done (including all spectra lines). 

There are some important distinctions in the type of calculations: 

* LTE versus non-LTE: whether the occupation of the atoms' energy levels are computed using thermodynamical equilibrium (i.e. the Saha and Boltzmann equations, which assumes that collision is the only mechanism responsible for the electrons moving around the energy levels -- not good for hot stars)
* Spherical versus plane-parallel (aka flat): whether the curvature of the atmosphere is taken into account -- this is also relevent for codes that can model the spectral lines that originates in the stellar winds. 

