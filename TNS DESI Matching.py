#!/usr/bin/env python
# coding: utf-8

# In[2]:


### Match TNS to DESI

A little blurb about the notebook.


# In[3]:


import os
import numpy as np
import matplotlib.pyplot as plt

from astropy import units as u
from astropy.table import Table, vstack, hstack
from astropy.coordinates import SkyCoord, match_coordinates_sky

from astropy.io.fits import getdata
from matplotlib.projections import get_projection_names
from astropy.coordinates import ICRS

from desiutil.plots import prepare_data, init_sky, plot_grid_map, plot_healpix_map, plot_sky_circles, plot_sky_binned
get_ipython().run_line_magic('matplotlib', 'inline')


# #### Read TNS table (generated elsewhere).

# In[5]:


def make_skyplot(tns, zcat=None, png=None):
    """Write me."""
    
    
    fig = plt.figure(figsize=(10, 8))
    ax = init_sky(projection = 'aitoff', ra_center = 90, galactic_plane_color = 'Red')
    ax.scatter(ax.projection_ra(tns['RA']), ax.projection_dec(tns['DEC']), 
               marker='.', color = 'Blue', s=3, alpha = .35)
            
    if zcat:
        
        MIa = np.where(zcat_match['CLASS'] == 'SN Ia')
        MII = np.where(zcat_match['CLASS'] == 'SN II')
        MMisc = np.where(zcat_match['CLASS'] == 'MISC')
        
        
        ax.scatter(ax.projection_ra(zcat['TARGET_RA'][MIa]), ax.projection_dec(zcat['TARGET_DEC'][MIa]), 
               marker='s', color = 'green', s=100, facecolor='none', label='Type Ia Supernovae')
        ax.scatter(ax.projection_ra(zcat['TARGET_RA'][MII]), ax.projection_dec(zcat['TARGET_DEC'][MII]), 
               marker='s', color = 'orange', s=100, facecolor='none', label='Type II Supernovae')
        ax.scatter(ax.projection_ra(zcat['TARGET_RA'][MMisc]), ax.projection_dec(zcat['TARGET_DEC'][MMisc]), 
               marker='s', color = 'red', s=100, facecolor='none', label='Miscellaneous')
        
        plt.legend()
        plt.show()
    
    if png:
        print('Writing {}'.format(png))
        plt.savefig(png)


# In[6]:


def read_tns(makeplot=False, png=None):
    from astropy.time import Time

    tnsfile = 'tns_all_thru_20220615.csv'
    print('Reading {}'.format(tnsfile))
    tns = Table.read(tnsfile, format='csv')
    [tns.rename_column(col, col.upper()) for col in tns.colnames]
    
    tns_coord = SkyCoord(tns['RA'], tns['DEC'], frame='icrs', unit=['hour', 'degree'])
    tns.rename_column('RA', 'RA_HMS')
    tns.rename_column('DEC', 'DEC_DMS')
    tns['RA'] = tns_coord.ra.value
    tns['DEC'] = tns_coord.dec.value
                      
    # Convert discovery date to MJD.
    mjd = Time(tns['DISCOVERY DATE (UT)'], format='iso', scale='utc').mjd
    tns['DISCOVERY DATE (MJD)'] = mjd
    
    tns['CLASS'] = np.zeros(len(tns), dtype='U5')
    
    indx_I = np.array(['SN Ia' in typ for typ in tns['OBJ. TYPE']])
    tns['CLASS'][indx_I] = 'SN Ia'
    
    indx_II = np.array(['SN II' in typ for typ in tns['OBJ. TYPE']])
    tns['CLASS'][indx_II] = 'SN II'
    
    indx_misc = np.array(['SN II' not in typ and 'SN Ia' not in typ for typ in tns['OBJ. TYPE']])
    tns['CLASS'][indx_misc] = 'MISC'
    
    if makeplot:
         make_skyplot(tns, png=png)
    
    return tns


# In[7]:


get_ipython().run_line_magic('time', "tns = read_tns(makeplot=True, png='tns-allsky.png')")
tns


# In[8]:


def read_desi_zcat(specprod='fuji'):
    """Write me."""
    
    import fitsio
    from astropy.table import Table, vstack
    from glob import glob
    from desitarget.targets import main_cmx_or_sv
    
    specprod_dir = f'/global/cfs/cdirs/desi/spectro/redux/{specprod}/zcatalog'
    zcatfiles = glob(os.path.join(specprod_dir, 'zpix-*.fits'))

    zcat1 = []
    for zcatfile in zcatfiles:
        print('Reading {}'.format(zcatfile))
        zcat1 = Table(fitsio.read(zcatfile))
        zcat1['ISBGS'] = np.zeros(len(zcat1), bool)
        zcat1['ISELG'] = np.zeros(len(zcat1), bool)
        zcat1['ISLRG'] = np.zeros(len(zcat1), bool)
        zcat1['ISQSO'] = np.zeros(len(zcat1), bool)

        # determine each type of object
        for iobj in np.arange(len(zcat1)):
            targetcols, targetmasks, survey = main_cmx_or_sv(zcat1[iobj])
            zcat1['ISBGS'][iobj] = zcat1[iobj]['{}_DESI_TARGET'.format(survey.upper())] & targetmasks[0].BGS_ANY != 0 
            zcat1['ISELG'][iobj] = zcat1[iobj]['{}_DESI_TARGET'.format(survey.upper())] & targetmasks[0].ELG != 0 
            zcat1['ISLRG'][iobj] = zcat1[iobj]['{}_DESI_TARGET'.format(survey.upper())] & targetmasks[0].LRG != 0 
            zcat1['ISQSO'][iobj] = zcat1[iobj]['{}_DESI_TARGET'.format(survey.upper())] & targetmasks[0].QSO != 0        
        zcat.append(zcat1)
        
    zcat = vstack(zcat)

    return zcat

# In[9]:


get_ipython().run_line_magic('time', 'zcat= read_desi_zcat()')
zcat


# In[11]:


def match_tns_desi(tns, zcat, rad=1.0, makeplot=False):
    """Write me."""
    
    tns_coord = SkyCoord(tns['RA'], tns['DEC'], frame='icrs', unit='degree')
    desi_coord = SkyCoord(zcat['TARGET_RA'], zcat['TARGET_DEC'], frame='icrs', unit='degree')
    
    minsep = rad * u.arcsec
    idx, d2d, _ = desi_coord.match_to_catalog_3d(tns_coord)  
    indx_zcat = np.where(d2d < minsep)[0]
    indx_tns = idx[indx_zcat]
    print('{} matches'.format(len(indx_zcat)))
    
    tns_match = tns[indx_tns]
    zcat_match = zcat[indx_zcat]
    
    zcat_match['CLASS'] = tns_match['CLASS']
    
    if makeplot:
        pngfile = 'tns-desi-matches.png'
        fig = plt.figure(figsize=(10, 8))
        plt.scatter(tns_coord.ra.value[indx_tns], tns_coord.dec.value[indx_tns], marker='s', s=100, color='blue', label='TNS')
        plt.scatter(desi_coord.ra.value[indx_zcat], desi_coord.dec.value[indx_zcat], marker='x', s=120, color='orange', label='DESI')
        
        
        plt.xlabel('RA [deg]')
        plt.ylabel('DEC [deg]')
        
        plt.title('DEC vs RA for TNS & DESI Matches')
        plt.legend()
        plt.show()
        
        
        print('Writing {}'.format(pngfile))
        plt.savefig(pngfile)
    
    return zcat_match, tns_match


# In[15]:


get_ipython().run_line_magic('time', 'zcat_match, tns_match = match_tns_desi(tns, zcat, rad=5.0, makeplot=True)')


# In[16]:


make_skyplot(tns, zcat_match, png='tns-allsky-desi-matches.png')


# In[14]:


set(zcat_match['CLASS'])

