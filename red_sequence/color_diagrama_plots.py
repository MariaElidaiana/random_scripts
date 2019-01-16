from __future__ import division
import os
import numpy as np
from astropy.table import Table, vstack
import matplotlib.pyplot as plt
import healpy as hp
import pandas as pd
from time import clock
import time

redpath = '/data/des50.b/data/redmapper/v2_2_1/'   #new run DES Y3, change the files below
redclusters = Table.read(redpath + 'y3_gold_2.2.1_wide_sofcol_run2_redmapper_v6.4.22+2_lgt20_vl02_catalog.fit')
redmembers  = Table.read(redpath + 'y3_gold_2.2.1_wide_sofcol_run2_redmapper_v6.4.22+2_lgt20_vl02_catalog_members.fit')

# clusters
c_id = redclusters['MEM_MATCH_ID']
ra = redclusters['RA']
dec = redclusters['DEC']
z = redclusters['Z_LAMBDA']
lambda_ = redclusters['LAMBDA_CHISQ']
mag_g = redclusters['MODEL_MAG'][:,0]
mag_r = redclusters['MODEL_MAG'][:,1]
mag_i = redclusters['MODEL_MAG'][:,2]
mag_z = redclusters['MODEL_MAG'][:,3]
gr = mag_g - mag_r
ri = mag_r - mag_i
iz = mag_i - mag_z

#members
m_id = redmembers['MEM_MATCH_ID']
m_mag_g = redmembers['MODEL_MAG'][:,0]
m_mag_r	= redmembers['MODEL_MAG'][:,1]
m_mag_i	= redmembers['MODEL_MAG'][:,2]
m_mag_z	= redmembers['MODEL_MAG'][:,3] 
m_magerr_g = redmembers['MODEL_MAGERR'][:,0]
m_magerr_r = redmembers['MODEL_MAGERR'][:,1]
m_magerr_i = redmembers['MODEL_MAGERR'][:,2]
m_magerr_z = redmembers['MODEL_MAGERR'][:,3]
m_gr = m_mag_g - m_mag_r
m_ri = m_mag_r - m_mag_i
m_iz = m_mag_i - m_mag_z
m_gr_err = np.sqrt(m_magerr_g**2 + m_magerr_r**2)
m_ri_err = np.sqrt(m_magerr_r**2 + m_magerr_i**2)
m_iz_err = np.sqrt(m_magerr_i**2 + m_magerr_z**2)


def colordiagram(mags, colors, zs, labelx, labely, labelcmap, limx, limy, outname):
    
    fig, axs = plt.subplots(1, 3, figsize=(24,10))
    for i in range(len(mags)):
        if 'lambda' in outname:
            im=axs[i].hexbin(mags[i], colors[i], C=np.log10(zs), vmin=round( np.log10(min(zs)) ,2) , vmax=round( np.log10(max(zs)) ,2), cmap=plt.cm.jet)
            labelcmap='$\log{\lambda}$'
        else: 
            im=axs[i].hexbin(mags[i], colors[i], C=zs, vmin=round(min(zs),2) , vmax=round(max(zs),2), cmap=plt.cm.jet)
        axs[i].set_xlabel(labelx[i])
        axs[i].set_ylabel(labely[i])
        #axs[i].set_ylim(limy)
        #axs[i].set_xlim((13,23))
        axs[i].set_xlim(limx)
    cb = plt.colorbar(im, ax=axs.ravel().tolist(), orientation="horizontal", shrink=0.3)
    cb.set_label(labelcmap, labelpad=+0.5)
    plt.savefig(outname, bbox_inches='tight')
    plt.close()

def colorz(zs, colors, richness, labely, limy, fgsize, outname, plottype):

    fig, axs = plt.subplots(1, 3, figsize=fgsize)
    #axs = axs.ravel()
    for i in range(len(colors)):
        if plottype=='scatter':
            axs[i].plot(zs, colors[i],'r.', alpha=0.3)
            axs[i].set_xlabel('z')
            axs[i].set_xlim(0,1)
            axs[i].set_ylabel(labely[i])
        elif plottype=='hex_density':
            im=axs[i].hexbin(zs, colors[i], cmap=plt.cm.jet)
            axs[i].cmap=plt.cm.jet
            axs[i].set_xlabel('z')
            axs[i].set_xlim(0,1) 
            axs[i].set_ylabel(labely[i])
            #axs[i].set_ylim(limy)
        elif plottype=='hex_lambda':
            im=axs[i].hexbin(zs, colors[i], C=np.log10(richness), vmin=round( np.log10(min(richness)) ,2) , vmax=round( np.log10(max(richness)) ,2), cmap=plt.cm.jet)
            axs[i].set_xlabel('z')
            axs[i].set_xlim(0,1)
            axs[i].set_ylabel(labely[i])
            if 'stdcolorz' in outname:
                axs[i].axhline(np.mean(colors[i]) , color='k', linestyle='-')
                axs[i].text(0.8, np.mean(colors[i]), '$\mu$= '+str(round(np.mean(colors[i]),3)), fontsize=12, horizontalalignment='left', verticalalignment='top', fontweight='bold')
                axs[i].axhline(np.mean(colors[i]) + 3*np.std(colors[i]) , color='k', linestyle='--')
                axs[i].text(0.8, np.mean(colors[i]) + 3*np.std(colors[i]), '$3\sigma$= '+str(round(3*np.std(colors[i]),3)), fontsize=12, horizontalalignment='left', verticalalignment='bottom', fontweight='bold')
                #axs[i].axhline(np.mean(colors[i]) - np.std(colors[i]) , color='r', linestyle='--')
                idxs=(colors[i] > 3*np.std(colors[i]) )
                N_out=len(np.array(colors)[i][idxs])
                N=len(colors[i])
                frac = (N_out*100)/N
                axs[i].set_title('$frac_{out}$= '+str(round(frac,2))+'$\%$')
                #axs[i].text(0.5, 0.5, '$frac_{out}$= '+str(round(frac,2))+'$\%$', fontsize=12, horizontalalignment='center', verticalalignment='center', fontweight='bold', transform=axs[i].transAxes)

    if plottype=='hex_lambda':
        cb = plt.colorbar(im, ax=axs.ravel().tolist(), orientation="horizontal", shrink=0.3)
        cb.set_label(r'$\log{\lambda}$', labelpad=+0.5)
    plt.savefig(outname, bbox_inches='tight')
    plt.close()


start = time.time()

print "Clusters plot...."
#Plotting clusters
colordiagram([mag_r, mag_i, mag_z], [gr, ri, iz], z, ['r', 'i', 'z'], ['g - r', 'r - i', 'i - z'], 'z', (13,23), (0,2.5), 'colordiag_cluster_test_z.png')
colordiagram([mag_r, mag_i, mag_z], [gr, ri, iz], lambda_, ['r', 'i', 'z'], ['g - r', 'r - i', 'i - z'], r'$\lambda$', (13,23), (0,2.5), 'colordiag_cluster_test_lambda.png')
#colorz(z, [gr, ri, iz], lambda_, ['g - r', 'r - i', 'i - z'],  (0.25,2), (16,4), 'colorz_cluster_test_scatter.png', 'scatter')
colorz(z, [gr, ri, iz], lambda_, ['g - r', 'r - i', 'i - z'], (0.25,2), (24,10), 'colorz_cluster_test_hex_lambda.png', 'hex_lambda')
colorz(z, [gr, ri, iz], lambda_, ['g - r', 'r - i', 'i - z'], (0.25,2), (16,4), 'colorz_cluster_test_hex_density.png', 'hex_density')


def groupby(a, b):
    # Get argsort indices, to be used to sort a and b in the next steps
    sidx = b.argsort(kind='mergesort')
    a_sorted = a[sidx]
    b_sorted = b[sidx]

    # Get the group limit indices (start, stop of groups)
    cut_idx = np.flatnonzero(np.r_[True,b_sorted[1:] != b_sorted[:-1],True])

    # Split input array with those start, stop ones
    out = [a_sorted[i:j] for i,j in zip(cut_idx[:-1],cut_idx[1:])]
    return out


def getslope(clus_id, mem_id, mem_mag, mem_color, mem_ecolor ):

    #slope, yint, mean_color, std_color, mean_mag = [],[],[],[],[]
    slope, yint = [],[]

    band = groupby(mem_mag, mem_id)
    color = groupby(mem_color, mem_id)
    ecolor = groupby(mem_ecolor, mem_id)

    mean_mag = [np.mean(i) for i in band] 
    mean_color = [np.mean(i) for i in color]
    std_color = [np.std(i) for i in color]

    for (i, j, k) in zip(band, color, ecolor):
        rfit,info = np.polynomial.polynomial.polyfit(i,j,deg=1,full=True,w=1./k) #i=band, j=color, k=ecolor
        slope.append(rfit[1])
        yint.append(rfit[0]) 
    return np.array(slope), mean_mag, mean_color, std_color

gr_getslope = getslope(c_id, m_id, m_mag_r, m_gr, m_gr_err)
gr_slope    = gr_getslope[0]
magr_mean   = gr_getslope[1]
gr_mean     = gr_getslope[2]
gr_std      = gr_getslope[3]

ri_getslope = getslope(c_id, m_id, m_mag_i, m_ri, m_ri_err)
ri_slope    = ri_getslope[0]
magi_mean   = ri_getslope[1]
ri_mean     = ri_getslope[2]
ri_std      = ri_getslope[3]

iz_getslope = getslope(c_id, m_id, m_mag_z, m_iz, m_iz_err)
iz_slope    = iz_getslope[0]
magz_mean   = iz_getslope[1]
iz_mean     = iz_getslope[2]
iz_std      = iz_getslope[3]


print "Plotting members...."

#Plotting members
colordiagram([magr_mean, magi_mean, magz_mean], [gr_mean, ri_mean, iz_mean], z, ['r', 'i', 'z'], ['g - r', 'r - i', 'i - z'], 'z',(15,25), (0,2.5), 'colordiag_mem_test_z.png')
colordiagram([magr_mean, magi_mean, magz_mean], [gr_mean, ri_mean, iz_mean], lambda_, ['r', 'i', 'z'], ['g - r', 'r - i', 'i - z'], r'$\lambda$',(15,25), (0,2.5), 'colordiag_mem_test_lambda.png')

#colorz(z, [gr_mean, ri_mean, iz_mean], lambda_, ['g - r', 'r - i', 'i - z'], (0.25,2), (16,4), 'colorz_mem_test_scatter.png', 'scatter')
colorz(z, [gr_mean, ri_mean, iz_mean], lambda_, ['g - r', 'r - i', 'i - z'], (0.25,2), (24,10), 'colorz_mem_test_hex_lambda.png', 'hex_lambda')
colorz(z, [gr_mean, ri_mean, iz_mean], lambda_, ['g - r', 'r - i', 'i - z'], (0.25,2), (16,4), 'colorz_mem_test_hex_density.png', 'hex_density')

#colorz(z, [gr_std, ri_std, iz_std], lambda_, [r'$\sigma_{g - r}$', r'$\sigma_{r - i}$', r'$\sigma_{i - z}$'], (0,0.4), (16,4), 'stdcolorz_mem_test_scatter.png', 'scatter')
colorz(z, [gr_std, ri_std, iz_std], lambda_, [r'$\sigma_{g - r}$', r'$\sigma_{r - i}$', r'$\sigma_{i - z}$'], (0,0.4), (24,10), 'stdcolorz_mem_test_hex_lambda.png', 'hex_lambda')
colorz(z, [gr_std, ri_std, iz_std], lambda_, [r'$\sigma_{g - r}$', r'$\sigma_{r - i}$', r'$\sigma_{i - z}$'], (0,0.4), (16,4), 'stdcolorz_mem_test_hex_density.png', 'hex_density')

#colorz(z, [gr_slope, ri_slope, iz_slope], lambda_, ['g-r slope', 'r-i slope', 'i-z slope'], (-0.3,0.3), (16,4), 'slopez_mem_test_scatter.png','scatter')
colorz(z, [gr_slope, ri_slope, iz_slope], lambda_, ['g-r slope', 'r-i slope', 'i-z slope'], (-0.3,0.3), (24,10), 'slopez_mem_test_hex_lambda.png','hex_lambda')
colorz(z, [gr_slope, ri_slope, iz_slope], lambda_, ['g-r slope', 'r-i slope', 'i-z slope'], (-0.3,0.3), (16,4), 'slopez_mem_test_hex_density.png','hex_density')

end = time.time()
print 'Time to do the plots:', end - start, 'seconds'
