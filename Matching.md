### Match TNS to DESI

A little blurb about the notebook.


```python
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
%matplotlib inline
```

#### Read TNS table (generated elsewhere).


```python
def make_skyplot(tns, zcat=None):
    """Write me."""
    
    pngfile = 'tns-allsky.png'
    fig = plt.figure(figsize=(10, 8))
    ax = init_sky(projection = 'aitoff', ra_center = 90, galactic_plane_color = 'Red')
    ax.scatter(ax.projection_ra(tns_coord.ra.degree), ax.projection_dec(tns_coord.dec.degree), 
               marker='.', color = 'Black', s=5)
    print('Writing {}'.format(pngfile))
    plt.savefig(pngfile)   
```


```python
def read_tns(makeplot=False):
    from astropy.time import Time

    tnsfile = 'tns_all_thru_20220615.csv'
    print('Reading {}'.format(tnsfile))
    tns = Table.read(tnsfile, format='csv')
    [tns.rename_column(col, col.upper()) for col in tns.colnames]
    
    tns_coord = SkyCoord(tns['RA'], tns['DEC'], frame='icrs', unit=['hour', 'degree'])
    tns.rename_column('RA', 'RA_HMS')
    tns.rename_column('DEC', 'DEC_DMS')
    tns['RA'] = tns_coord.ra.value
    tns['DEC'] = tns_coord.ra.value
                      
    # Convert discovery date to MJD.
    mjd = Time(tns['DISCOVERY DATE (UT)'], format='iso', scale='utc').mjd
    tns['DISCOVERY DATE (MJD)'] = mjd
    
    #if makeplot:
    
    return tns
```


```python
%time tns = read_tns(makeplot=True)
tns
```

```python
def read_desi_zcat(specprod='fuji'):
    """Write me."""
    import fitsio
    from astropy.table import Table, vstack
    from glob import glob
    
    specprod_dir = f'/global/cfs/cdirs/desi/spectro/redux/{specprod}/zcatalog'
    zcatfiles = glob(os.path.join(specprod_dir, 'zpix-*.fits'))

    zcat = []
    for zcatfile in zcatfiles:
        print('Reading {}'.format(zcatfile))
        zcat1 = Table(fitsio.read(zcatfile))
        zcat.append(zcat1)
        
    zcat = vstack(zcat)

    return zcat   
```


```python
%time zcat = read_desi_zcat()
zcat
```
