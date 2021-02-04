#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')

import numpy as np
import h5py
import mpl_scatter_density
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u

###### plot style
matplotlib.rc('xtick', labelsize=25, top=True, direction='in')
matplotlib.rc('ytick', labelsize=25, right=True, direction='in')
matplotlib.rc('axes', linewidth=1, labelsize=25)
matplotlib.rc('xtick.major', size=5)
matplotlib.rc('xtick.minor', size=5)
matplotlib.rc('ytick.major', size=5)
matplotlib.rc('ytick.minor', size=5)
plt.rc('text', usetex=True)
plt.rcParams['figure.figsize'] = [8, 8]
######

fname ='/data/des81.b/data/mariaeli/y3_buzz/Buzzard-3_v1.9.8_Y3a_mastercat/Buzzard-3_v1.9.8_Y3a_mastercat.h5'
f     = h5py.File(fname, 'r')

select_idx = f['index/select'][:] #this selects all the source galaxies for science analysis
ra  = f['catalog/metacal/unsheared/ra'][...][select_idx]
dec = f['catalog/metacal/unsheared/dec'][...][select_idx]
f.close()

#Getting the sky coordinates and changing for the referecial -180 to 180
c = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree))
ra[ra>180] -=360     # put in the range [-180, 180]

#Plotting sky distribution using https://pypi.org/project/mpl-scatter-density/
fig = plt.figure(figsize=(20, 10))
ax = fig.add_subplot(1,1,1, projection='scatter_density')
ax.scatter_density(ra,dec)
ax.grid(True)
ax.invert_xaxis()
plt.xlabel('$RA$')
plt.ylabel('$DEC$')
plt.title('Buzzard v1.9.8, N$_{src}$='+str(len(ra)), fontsize=25 )
fig.savefig('sky_sources_v1.9.8.png')
