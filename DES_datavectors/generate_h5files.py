import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
from datetime import datetime
import six

nvec = 12
ndim = 15
rad_vectors = np.zeros(nvec*ndim )#+ ndim)
data_vectors = np.zeros(nvec*ndim )#+ ndim)
#data_vectors_err = np.zeros(nvec*ndim + ndim)
#cov_vectors = np.zeros([nvec*ndim+ndim, nvec*ndim+ndim])
cov_vectors = np.zeros([nvec*ndim, nvec*ndim])

#NC data-vector
ncs = np.genfromtxt('data_files/nc.txt', delimiter=',')
#ncs_err = np.genfromtxt('data_files/nc_err.txt')
nc_covmat = np.genfromtxt('data_files/Cov_ij_bestfit_DESY1_105.txt')

#WL data-vector
for l in range(4):
    for z in range(3):
        #data_vectors[z+l*3] = ncs[z, l]
        #cov_vectors[z+l*3, z+l*3] = ncs_err[z, l]**2
        lam = l+3
        dat = np.genfromtxt('data_files/wl_data_files/full-unblind-v2-mcal-zmix_y1subtr_l%i_z%i_profile.dat'%(lam, z));
        cov = np.genfromtxt('data_files/SAC_files/SAC_z%i_l%i.txt'%(z, lam))
        start = (z + l * 3) * ndim #+ ndim
        stop  = (z + l * 3) * ndim + ndim #+ ndim
        print(dat.shape, cov.shape, start, stop,  np.shape(rad_vectors[start:stop]))
        rad_vectors[start:stop]=dat[:, 0]
        data_vectors[start:stop]=dat[:, 1]
        #data_vectors_err[start:stop]=dat[:, 2]
        cov_vectors[start:stop, start:stop]=cov

now = datetime.now()
timestamp = now.strftime("%d/%m/%Y %H:%M:%S")

# create the HDF5 datavector file
fileName = 'y1clusters_datavector.h5'
f = h5.File(fileName, "w")

## create the y1 group
dvector = f.create_group(u'y1')
# give the HDF5 root some more attributes
f.attrs[u'file_name']        = fileName
f.attrs[u'file_time']        = timestamp
f.attrs[u'creator']          = u'dasilvap@umich.edu'
f.attrs[u'HDF5_Version']     = six.u(h5.version.hdf5_version)
f.attrs[u'h5py_version']     = six.u(h5.version.version)

# create the dvector group
dvdata = dvector.create_group(u'dvector')
dvdata.attrs[u'deltasig'] = u'deltasig'
dvdata.attrs[u'radius'] = u'radius'
dvdata.attrs[u'cov_wl'] = u'covmat_wl'
dvdata.attrs[u'ncounts'] = u'ncounts'
dvdata.attrs[u'cov_nc'] = u'covmat_nc'

print(ncs)

# create dataset
ds = dvdata.create_dataset(u'radius', data=rad_vectors)
ds.attrs[u'R units'] = u'Mpc [physical]'

ds = dvdata.create_dataset(u'deltasig', data=data_vectors)
ds.attrs[u'DS units'] = u'Msun/pc^2'

ds = dvdata.create_dataset(u'covmat_wl', data=cov_vectors)
ds.attrs[u'COV_WL units'] = u'Msun/pc^2'

ds = dvdata.create_dataset(u'ncounts', data=ncs)
ds.attrs[u'COV_NC units'] = u'counts'

ds = dvdata.create_dataset(u'covmat_nc', data=nc_covmat)
ds.attrs[u'NC units'] = u'counts'

f.close()   # be CERTAIN to close the file

print("wrote file:", fileName)
