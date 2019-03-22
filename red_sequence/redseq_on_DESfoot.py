from __future__ import division
import os
import numpy as np
from astropy.table import Table, vstack
import matplotlib.pyplot as plt
import healpy as hp
import skymapper as skm
import pandas as pd
from time import clock
import time

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

def plot_hpixdist(ra, dec, barvalue, vlim, bartitle, outname):
    # setup figure
    import matplotlib.cm as cm
    #cmap = cm.RdYlBu_r
    cmap = cm.YlOrRd
    
    nside = 1024
    sep=15
    
    fig = plt.figure(figsize=(6.5,6))
    ax = fig.add_subplot(111, aspect='equal')
    
    proj = skm.createConicMap(ax, ra, dec, proj_class=skm.AlbersEqualAreaProjection)
    # add lines and labels for meridians/parallels (separation 5 deg)
    
    meridians = np.arange(-90, 90+sep, sep)
    parallels = np.arange(0, 360+sep, sep*2)
    skm.setMeridianPatches(ax, proj, meridians, linestyle=':', lw=0.5, zorder=2)
    skm.setParallelPatches(ax, proj, parallels, linestyle=':', lw=0.5, zorder=2)
    skm.setMeridianLabels(ax, proj, meridians, loc="left", fmt=skm.pmDegFormatter)
    skm.setParallelLabels(ax, proj, parallels, loc="bottom")#, fmt=skm.hourAngleFormatter)
    
    # convert to map coordinates and plot a marker for each point
    x,y = proj(ra, dec)
    marker = 's'
    markersize = skm.getMarkerSizeToFill(fig, ax, x, y)
    vmin,vmax = min(barvalue), max(barvalue)
    sc = ax.scatter(x,y, c=barvalue, edgecolors='None', marker=marker, s=markersize, cmap=cmap, vmin=vmin, vmax=vmax, rasterized=True, zorder=1)
    
    # add colorbar
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.0)
    cb = fig.colorbar(sc, cax=cax)
    cb.set_label(bartitle)
    ticks = np.linspace(vmin, vmax, 5)
    cb.set_ticks(ticks)
    cb.solids.set_edgecolor("face")
    
    # show (and save) ...
    fig.tight_layout()
    fig.savefig(outname, bbox_inches='tight')

if __name__ == "__main__":

    redpath = '/data/des50.b/data/redmapper/v2_2_1/'   #new run DES Y3
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
    
    '''
    #members
    m_id = redmembers['MEM_MATCH_ID']
    m_mag_g = redmembers['MODEL_MAG'][:,0]
    m_mag_r = redmembers['MODEL_MAG'][:,1]
    m_mag_i = redmembers['MODEL_MAG'][:,2]
    m_mag_z = redmembers['MODEL_MAG'][:,3]
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

    #get red sequences properties
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
    '''

    #Plotting the cluster density on DES footprint
    nside = 64
    sep = 15
    fig = plt.figure(figsize=(8,4))
    ax = fig.add_subplot(111, aspect='equal')
    fig, ax, proj = skm.plotDensity(ra, dec, nside=nside, sep=sep, ax=ax)
    #add DES footprint
    skm.addFootprint('DES', proj, ax, zorder=10, edgecolor='#2222B2', facecolor='None', lw=2)
    plt.savefig('density_map_on_DESfootprint.png', bbox_inches='tight')

    #Plotting the red-sequence properties on DES footprint
    '''
    plot_hpixdist(ra, dec, gr_mean, (0, 2), 'g-r mean', 'gr_mean_over_footprint.png')
    plot_hpixdist(ra, dec, ri_mean, (0, 2), 'r-i mean', 'ri_mean_over_footprint.png')
    plot_hpixdist(ra, dec, iz_mean, (0, 2), 'i-z mean', 'iz_mean_over_footprint.png')

    plot_hpixdist(ra, dec, gr_std, (0, 0.4), 'g-r std', 'gr_std_over_footprint.png')
    plot_hpixdist(ra, dec, ri_std, (0, 0.4), 'r-i std', 'ri_std_over_footprint.png')
    plot_hpixdist(ra, dec, iz_std, (0, 0.4), 'i-z std', 'iz_std_over_footprint.png')

    plot_hpixdist(ra, dec, gr_slope, (-0.3, 0.3), 'g-r slope', 'gr_slope_over_footprint.png')
    plot_hpixdist(ra, dec, ri_slope, (-0.3, 0.3), 'r-i slope', 'ri_slope_over_footprint.png')
    plot_hpixdist(ra, dec, iz_slope, (-0.3, 0.3), 'i-z slope', 'iz_slope_over_footprint.png')
    '''
    print "Done."


