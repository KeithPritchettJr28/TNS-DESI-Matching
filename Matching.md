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

    Reading tns_all_thru_20220615.csv
    CPU times: user 2.3 s, sys: 30.5 ms, total: 2.34 s
    Wall time: 2.33 s





<div><i>Table length=5766</i>
<table id="table23453862905264" class="table-striped table-bordered table-condensed">
<thead><tr><th>COL0</th><th>ID</th><th>NAME</th><th>RA_HMS</th><th>DEC_DMS</th><th>OBJ. TYPE</th><th>REDSHIFT</th><th>HOST NAME</th><th>HOST REDSHIFT</th><th>REPORTING GROUP/S</th><th>DISCOVERY DATA SOURCE/S</th><th>CLASSIFYING GROUP/S</th><th>ASSOCIATED GROUP/S</th><th>DISC. INTERNAL NAME</th><th>DISC. INSTRUMENT/S</th><th>CLASS. INSTRUMENT/S</th><th>TNS AT</th><th>PUBLIC</th><th>END PROP. PERIOD</th><th>DISCOVERY MAG/FLUX</th><th>DISCOVERY FILTER</th><th>DISCOVERY DATE (UT)</th><th>SENDER</th><th>REMARKS</th><th>DISCOVERY BIBCODE</th><th>CLASSIFICATION BIBCODES</th><th>EXT. CATALOG/S</th><th>RA</th><th>DEC</th><th>DISCOVERY DATE (MJD)</th></tr></thead>
<thead><tr><th>int64</th><th>int64</th><th>str11</th><th>str12</th><th>str12</th><th>str17</th><th>float64</th><th>str30</th><th>float64</th><th>str58</th><th>str51</th><th>str39</th><th>str51</th><th>str29</th><th>str154</th><th>str67</th><th>int64</th><th>int64</th><th>int64</th><th>float64</th><th>str12</th><th>str23</th><th>str18</th><th>str134</th><th>str19</th><th>str166</th><th>int64</th><th>float64</th><th>float64</th><th>float64</th></tr></thead>
<tr><td>0</td><td>110525</td><td>SN 2022mox</td><td>13:34:44.847</td><td>+33:53:23.10</td><td>SN Ia</td><td>0.023199</td><td>WISEA J133444.74+335325.5</td><td>0.023199</td><td>ALeRCE, ATLAS</td><td>ZTF, ATLAS</td><td>ZTF</td><td>ALeRCE; ATLAS</td><td>ZTF22aanzset</td><td>ATLAS-HKO - ATLAS-02, ATLAS-MLO - ATLAS-01, P48 - ZTF-Cam</td><td>P60 - SEDM</td><td>1</td><td>1</td><td>--</td><td>18.0214</td><td>g-ZTF</td><td>2022-06-12 07:16:40.999</td><td>ALeRCE</td><td>--</td><td>2022TNSTR1643....1F</td><td>2022TNSCR1667....1S</td><td>--</td><td>203.68686249999996</td><td>203.68686249999996</td><td>59742.30325230324</td></tr>
<tr><td>1</td><td>110451</td><td>SN 2022mme</td><td>18:33:36.132</td><td>+41:29:21.84</td><td>SN Ia</td><td>0.08</td><td>--</td><td>--</td><td>ATLAS</td><td>ATLAS</td><td>TCD</td><td>ATLAS</td><td>ATLAS22qom</td><td>ATLAS-HKO - ATLAS-02, ATLAS-MLO - ATLAS-01</td><td>LT - SPRAT</td><td>1</td><td>1</td><td>--</td><td>19.012</td><td>orange-ATLAS</td><td>2022-06-09 11:04:24.960</td><td>ATLAS_Bot1</td><td>--</td><td>2022TNSTR1621....1T</td><td>2022TNSCR1653....1T</td><td>--</td><td>278.40055</td><td>278.40055</td><td>59739.4614</td></tr>
<tr><td>2</td><td>110447</td><td>SN 2022mma</td><td>14:39:01.507</td><td>+15:59:11.82</td><td>SN IIn</td><td>0.038</td><td>SDSS J143901.93+155923.1</td><td>0.037474</td><td>SGLF, ATLAS</td><td>ZTF, ATLAS</td><td>Global SN Project</td><td>ATLAS; SGLF</td><td>ZTF22aanwibf</td><td>ATLAS-HKO - ATLAS-02, ATLAS-MLO - ATLAS-01, P48 - ZTF-Cam</td><td>FTS - FLOYDS-S</td><td>1</td><td>1</td><td>--</td><td>19.431</td><td>r-ZTF</td><td>2022-06-10 05:05:38.000</td><td>Perez-Fournon</td><td>--</td><td>2022TNSTR1626....1P</td><td>--</td><td>--</td><td>219.75627916666664</td><td>219.75627916666664</td><td>59740.21224537037</td></tr>
<tr><td>3</td><td>110375</td><td>SN 2022mji</td><td>09:42:54.053</td><td>+31:51:03.60</td><td>SN II</td><td>0.004</td><td>--</td><td>--</td><td>ATLAS</td><td>ATLAS</td><td>ZTF</td><td>ATLAS</td><td>ATLAS22qnq</td><td>ATLAS-HKO - ATLAS-02, ATLAS-MLO - ATLAS-01</td><td>P60 - SEDM</td><td>1</td><td>1</td><td>--</td><td>17.361</td><td>orange-ATLAS</td><td>2022-06-09 06:21:10.080</td><td>ATLAS_Bot1</td><td>--</td><td>2022TNSTR1607....1T</td><td>2022TNSCR1666....1S</td><td>--</td><td>145.7252208333333</td><td>145.7252208333333</td><td>59739.2647</td></tr>
<tr><td>4</td><td>110327</td><td>SN 2022mhn</td><td>23:26:47.865</td><td>+25:41:13.21</td><td>SN Ia</td><td>0.022279</td><td>UGC 12599</td><td>0.022279</td><td>ALeRCE</td><td>ZTF</td><td>ZTF</td><td>ALeRCE</td><td>ZTF22aanppbi</td><td>P48 - ZTF-Cam</td><td>P60 - SEDM</td><td>1</td><td>1</td><td>--</td><td>17.0282</td><td>g-ZTF</td><td>2022-06-08 09:54:25.004</td><td>ALeRCE</td><td>--</td><td>2022TNSTR1593....1M</td><td>2022TNSCR1641....1S</td><td>--</td><td>351.6994375</td><td>351.6994375</td><td>59738.41278939815</td></tr>
<tr><td>5</td><td>110320</td><td>SN 2022mhg</td><td>16:15:02.056</td><td>+50:19:52.78</td><td>SN Ia</td><td>0.041</td><td>--</td><td>--</td><td>ATLAS, Pan-STARRS</td><td>ATLAS, Pan-STARRS</td><td>TCD</td><td>ATLAS; Pan-STARRS</td><td>ATLAS22qln</td><td>ATLAS-HKO - ATLAS-02, ATLAS-MLO - ATLAS-01, PS2 - GPC2</td><td>LT - SPRAT</td><td>1</td><td>1</td><td>--</td><td>19.156</td><td>orange-ATLAS</td><td>2022-06-07 09:50:23.136</td><td>ATLAS_Bot1</td><td>--</td><td>2022TNSTR1595....1T</td><td>2022TNSCR1653....1T</td><td>--</td><td>243.75856666666664</td><td>243.75856666666664</td><td>59737.40999</td></tr>
<tr><td>6</td><td>110318</td><td>SN 2022mhe</td><td>23:28:06.446</td><td>+08:45:44.71</td><td>SN Ia</td><td>0.02874</td><td>--</td><td>--</td><td>ATLAS, ZTF</td><td>ATLAS, ZTF</td><td>SCAT</td><td>ATLAS; ZTF</td><td>ATLAS22qif</td><td>ATLAS-HKO - ATLAS-02, P48 - ZTF-Cam</td><td>UH88 - SNIFS</td><td>1</td><td>1</td><td>--</td><td>17.994</td><td>cyan-ATLAS</td><td>2022-06-07 13:45:14.976</td><td>ATLAS_Bot1</td><td>--</td><td>2022TNSTR1595....1T</td><td>2022TNSCR1602....1A</td><td>--</td><td>352.0268583333333</td><td>352.0268583333333</td><td>59737.57309</td></tr>
<tr><td>7</td><td>110305</td><td>SN 2022mgr</td><td>18:32:37.027</td><td>+20:36:09.41</td><td>SN IIn</td><td>0.068</td><td>--</td><td>--</td><td>ALeRCE</td><td>ZTF</td><td>SGLF</td><td>ALeRCE</td><td>ZTF22aanjyae</td><td>P48 - ZTF-Cam</td><td>LT - SPRAT</td><td>1</td><td>1</td><td>--</td><td>20.317</td><td>r-ZTF</td><td>2022-06-07 07:28:29.004</td><td>ALeRCE</td><td>--</td><td>2022TNSTR1579....1M</td><td>2022TNSCR1652....1P</td><td>--</td><td>278.1542791666667</td><td>278.1542791666667</td><td>59737.31144680556</td></tr>
<tr><td>8</td><td>110165</td><td>SN 2022mbh</td><td>10:04:21.290</td><td>+54:38:40.67</td><td>SN Ia</td><td>0.05</td><td>--</td><td>--</td><td>ATLAS, GaiaAlerts</td><td>ATLAS, GaiaAlerts</td><td>SCAT</td><td>ATLAS; GaiaAlerts</td><td>ATLAS22qdm</td><td>ATLAS-HKO - ATLAS-02, ATLAS-MLO - ATLAS-01, Gaia - Gaia-photometric</td><td>UH88 - SNIFS</td><td>1</td><td>1</td><td>--</td><td>18.421</td><td>orange-ATLAS</td><td>2022-06-05 06:04:26.976</td><td>ATLAS_Bot1</td><td>--</td><td>2022TNSTR1558....1T</td><td>2022TNSCR1603....1D</td><td>--</td><td>151.08870833333333</td><td>151.08870833333333</td><td>59735.25309</td></tr>
<tr><td>9</td><td>110139</td><td>SN 2022mai</td><td>13:18:04.296</td><td>+66:11:40.06</td><td>SN Ia</td><td>0.07749</td><td>--</td><td>--</td><td>ZTF, ATLAS</td><td>ZTF, ATLAS</td><td>ZTF</td><td>ATLAS; ZTF</td><td>ZTF22aalzyzl</td><td>ATLAS-HKO - ATLAS-02, ATLAS-MLO - ATLAS-01, P48 - ZTF-Cam</td><td>P60 - SEDM</td><td>1</td><td>1</td><td>--</td><td>20.42</td><td>g-ZTF</td><td>2022-05-30 09:48:57.600</td><td>ZTF_Bot1</td><td>--</td><td>2022TNSTR1563....1F</td><td>2022TNSCR1641....1S</td><td>--</td><td>199.5179</td><td>199.5179</td><td>59729.409</td></tr>
<tr><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
<tr><td>256</td><td>43437</td><td>SN 2014ge</td><td>12:04:51.500</td><td>+26:59:46.60</td><td>SN Ib</td><td>0.00189</td><td>--</td><td>--</td><td>MASTER, iPTF</td><td>MASTER, iPTF</td><td>Padova-Asiago</td><td>iPTF; MASTER</td><td>MASTER OT J120451.50+265946.6</td><td>Other - Other, P48 - CFH12k</td><td>Ekar - AFOSC</td><td>1</td><td>1</td><td>--</td><td>13.9</td><td>Clear-</td><td>2014-10-28 20:59:20.000</td><td>MASTER</td><td>--</td><td>2019TNSTR1539....1G</td><td>2021TNSCR2418....1P</td><td>--</td><td>181.2145833333333</td><td>181.2145833333333</td><td>56958.87453703704</td></tr>
<tr><td>257</td><td>6469</td><td>SN 2014dr</td><td>03:22:46.390</td><td>+00:09:20.70</td><td>SN II</td><td>--</td><td>UGC 2705</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>PS1 - GPC1</td><td>--</td><td>0</td><td>1</td><td>--</td><td>16.6</td><td>--</td><td>2014-10-14 00:00:00.000</td><td>--</td><td>[Discoverer=Howerton, Drake et al. (Catalina Real-time Transient Survey); DiscoveryRef=CBET 4009; PositionRef=CBET 4009; Orig Type=II]</td><td>2021TNSTR.486....1H</td><td>--</td><td>--</td><td>50.69329166666666</td><td>50.69329166666666</td><td>56944.0</td></tr>
<tr><td>258</td><td>6436</td><td>SN 2014cl</td><td>02:16:09.101</td><td>-11:56:02.62</td><td>SN IIb</td><td>--</td><td>IC 217</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>PS1 - GPC1</td><td>--</td><td>0</td><td>1</td><td>--</td><td>17.5</td><td>--</td><td>2014-06-16 00:00:00.000</td><td>--</td><td>[Discoverer=Gonzalez, Pignata et al. (CHASE); DiscoveryRef=CBET 3950; PositionRef=CBET 3950; Orig Type=IIb]</td><td>2021TNSTR.635....1G</td><td>--</td><td>--</td><td>34.03792083333333</td><td>34.03792083333333</td><td>56824.0</td></tr>
<tr><td>259</td><td>6403</td><td>SN 2014bf</td><td>15:57:35.230</td><td>+18:01:47.60</td><td>SN IIP</td><td>--</td><td>MCG +03-41-7</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>PS1 - GPC1</td><td>--</td><td>0</td><td>1</td><td>--</td><td>17.4</td><td>--</td><td>2014-05-22 00:00:00.000</td><td>--</td><td>[Discoverer=Mo, Li, Li, Wang, Zhang (THU-NAOC Transient Survey); DiscoveryRef=CBET 3889; PositionRef=CBET 3889; Orig Type=IIP]</td><td>2021TNSTR.507....1A</td><td>--</td><td>--</td><td>239.39679166666664</td><td>239.39679166666664</td><td>56799.0</td></tr>
<tr><td>260</td><td>88250</td><td>SN 2013ld</td><td>16:26:42.918</td><td>+13:06:50.08</td><td>SN Ia</td><td>0.02741</td><td>--</td><td>--</td><td>Pan-STARRS</td><td>Pan-STARRS</td><td>ZTF</td><td>Pan-STARRS</td><td>PS21hrr</td><td>PS1 - GPC1</td><td>P60 - SEDM</td><td>1</td><td>1</td><td>--</td><td>19.93</td><td>z-Sloan</td><td>2013-08-18 08:09:36.000</td><td>PS1_Bot1</td><td>--</td><td>2021TNSTR2390....1C</td><td>2021TNSCR2547....1C</td><td>--</td><td>246.67882499999996</td><td>246.67882499999996</td><td>56522.34</td></tr>
<tr><td>261</td><td>6336</td><td>SN 2013hj</td><td>09:12:06.290</td><td>-15:25:45.98</td><td>SN II</td><td>--</td><td>MCG -02-24-3</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>PS1 - GPC1</td><td>--</td><td>0</td><td>1</td><td>--</td><td>14.0</td><td>--</td><td>2013-12-12 00:00:00.000</td><td>--</td><td>[Discoverer=Antezana, Pignata et al. (CHASE); DiscoveryRef=CBET 3757; PositionRef=CBET 3757; Orig Type=II]</td><td>2021TNSTR.535....1A</td><td>--</td><td>--</td><td>138.0262083333333</td><td>138.0262083333333</td><td>56638.0</td></tr>
<tr><td>262</td><td>6309</td><td>SN 2013gj</td><td>02:41:10.001</td><td>-21:01:29.21</td><td>SN II</td><td>--</td><td>Anon.</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>PS1 - GPC1</td><td>--</td><td>0</td><td>1</td><td>--</td><td>17.5</td><td>--</td><td>2013-11-10 00:00:00.000</td><td>--</td><td>[Discoverer=Drake et al. (Catalina Real-time Transient Survey); DiscoveryRef=CBET 3716; PositionRef=CBET 3716; Orig Type=II]</td><td>2021TNSTR.496....1A</td><td>--</td><td>--</td><td>40.29167083333333</td><td>40.29167083333333</td><td>56606.0</td></tr>
<tr><td>263</td><td>6201</td><td>SN 2013cj</td><td>17:04:52.951</td><td>+12:55:10.42</td><td>SN IIn</td><td>--</td><td>UGC 10685</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>PS1 - GPC1</td><td>--</td><td>0</td><td>1</td><td>--</td><td>18.1</td><td>--</td><td>2013-05-09 00:00:00.000</td><td>--</td><td>[Discoverer=Jin, Gao (Xingming Sky Survey); DiscoveryRef=CBET 3520; PositionRef=CBET 3520; Orig Type=IIn]</td><td>2021TNSTR.670....1J</td><td>--</td><td>--</td><td>256.22062916666664</td><td>256.22062916666664</td><td>56421.0</td></tr>
<tr><td>264</td><td>6158</td><td>SN 2013at</td><td>18:50:11.290</td><td>+30:20:59.50</td><td>SN IIP</td><td>--</td><td>Anon.</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>PS1 - GPC1</td><td>--</td><td>0</td><td>1</td><td>--</td><td>16.9</td><td>--</td><td>2013-03-19 00:00:00.000</td><td>--</td><td>[Discoverer=Shurpakov (MASTER); DiscoveryRef=CBET 3448; PositionRef=CBET 3448; Orig Type=IIP]</td><td>2021TNSTR.748....1S</td><td>--</td><td>--</td><td>282.5470416666666</td><td>282.5470416666666</td><td>56370.0</td></tr>
<tr><td>265</td><td>85276</td><td>SN 2011km</td><td>08:09:12.866</td><td>+46:18:48.81</td><td>SN Ia</td><td>0.0459</td><td>SDSS J080912.92+461846.8</td><td>0.046309</td><td>--</td><td>iPTF</td><td>None</td><td>iPTF</td><td>PTF11kx</td><td>P48 - CFH12k</td><td>Lick-3m - KAST</td><td>1</td><td>1</td><td>--</td><td>19.91</td><td>R-Cousins</td><td>2011-01-16 08:47:32.640</td><td>Nugent</td><td>--</td><td>2021TNSTR1888....1N</td><td>2021TNSCR1897....1N</td><td>--</td><td>122.30360833333333</td><td>122.30360833333333</td><td>55577.36635</td></tr>
</table></div>




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

    Reading /global/cfs/cdirs/desi/spectro/redux/fuji/zcatalog/zpix-sv2-backup.fits
    Reading /global/cfs/cdirs/desi/spectro/redux/fuji/zcatalog/zpix-sv3-bright.fits
    Reading /global/cfs/cdirs/desi/spectro/redux/fuji/zcatalog/zpix-sv3-dark.fits
    Reading /global/cfs/cdirs/desi/spectro/redux/fuji/zcatalog/zpix-sv3-backup.fits
    Reading /global/cfs/cdirs/desi/spectro/redux/fuji/zcatalog/zpix-sv1-bright.fits
    Reading /global/cfs/cdirs/desi/spectro/redux/fuji/zcatalog/zpix-cmx-other.fits
    Reading /global/cfs/cdirs/desi/spectro/redux/fuji/zcatalog/zpix-sv2-dark.fits
    Reading /global/cfs/cdirs/desi/spectro/redux/fuji/zcatalog/zpix-sv1-dark.fits
    Reading /global/cfs/cdirs/desi/spectro/redux/fuji/zcatalog/zpix-sv1-other.fits
    Reading /global/cfs/cdirs/desi/spectro/redux/fuji/zcatalog/zpix-special-dark.fits
    Reading /global/cfs/cdirs/desi/spectro/redux/fuji/zcatalog/zpix-sv1-backup.fits
    Reading /global/cfs/cdirs/desi/spectro/redux/fuji/zcatalog/zpix-sv2-bright.fits
    CPU times: user 16.9 s, sys: 2.88 s, total: 19.8 s
    Wall time: 21.1 s





<div><i>Table length=2847435</i>
<table id="table23453863548768" class="table-striped table-bordered table-condensed">
<thead><tr><th>TARGETID</th><th>HEALPIX</th><th>SPGRPVAL</th><th>Z</th><th>ZERR</th><th>ZWARN</th><th>CHI2</th><th>COEFF [10]</th><th>NPIXELS</th><th>SPECTYPE</th><th>SUBTYPE</th><th>NCOEFF</th><th>DELTACHI2</th><th>COADD_FIBERSTATUS</th><th>TARGET_RA</th><th>TARGET_DEC</th><th>PMRA</th><th>PMDEC</th><th>REF_EPOCH</th><th>FA_TARGET</th><th>FA_TYPE</th><th>OBJTYPE</th><th>SUBPRIORITY</th><th>OBSCONDITIONS</th><th>RELEASE</th><th>BRICKNAME</th><th>BRICKID</th><th>BRICK_OBJID</th><th>MORPHTYPE</th><th>EBV</th><th>FLUX_G</th><th>FLUX_R</th><th>FLUX_Z</th><th>FLUX_W1</th><th>FLUX_W2</th><th>FLUX_IVAR_G</th><th>FLUX_IVAR_R</th><th>FLUX_IVAR_Z</th><th>FLUX_IVAR_W1</th><th>FLUX_IVAR_W2</th><th>FIBERFLUX_G</th><th>FIBERFLUX_R</th><th>FIBERFLUX_Z</th><th>FIBERTOTFLUX_G</th><th>FIBERTOTFLUX_R</th><th>FIBERTOTFLUX_Z</th><th>MASKBITS</th><th>SERSIC</th><th>SHAPE_R</th><th>SHAPE_E1</th><th>SHAPE_E2</th><th>REF_ID</th><th>REF_CAT</th><th>GAIA_PHOT_G_MEAN_MAG</th><th>GAIA_PHOT_BP_MEAN_MAG</th><th>GAIA_PHOT_RP_MEAN_MAG</th><th>PARALLAX</th><th>PHOTSYS</th><th>PRIORITY_INIT</th><th>NUMOBS_INIT</th><th>SV2_DESI_TARGET</th><th>SV2_BGS_TARGET</th><th>SV2_MWS_TARGET</th><th>SV2_SCND_TARGET</th><th>DESI_TARGET</th><th>BGS_TARGET</th><th>MWS_TARGET</th><th>PLATE_RA</th><th>PLATE_DEC</th><th>COADD_NUMEXP</th><th>COADD_EXPTIME</th><th>COADD_NUMNIGHT</th><th>COADD_NUMTILE</th><th>MEAN_DELTA_X</th><th>RMS_DELTA_X</th><th>MEAN_DELTA_Y</th><th>RMS_DELTA_Y</th><th>MEAN_FIBER_RA</th><th>STD_FIBER_RA</th><th>MEAN_FIBER_DEC</th><th>STD_FIBER_DEC</th><th>MEAN_PSF_TO_FIBER_SPECFLUX</th><th>TSNR2_GPBDARK_B</th><th>TSNR2_ELG_B</th><th>TSNR2_GPBBRIGHT_B</th><th>TSNR2_LYA_B</th><th>TSNR2_BGS_B</th><th>TSNR2_GPBBACKUP_B</th><th>TSNR2_QSO_B</th><th>TSNR2_LRG_B</th><th>TSNR2_GPBDARK_R</th><th>TSNR2_ELG_R</th><th>TSNR2_GPBBRIGHT_R</th><th>TSNR2_LYA_R</th><th>TSNR2_BGS_R</th><th>TSNR2_GPBBACKUP_R</th><th>TSNR2_QSO_R</th><th>TSNR2_LRG_R</th><th>TSNR2_GPBDARK_Z</th><th>TSNR2_ELG_Z</th><th>TSNR2_GPBBRIGHT_Z</th><th>TSNR2_LYA_Z</th><th>TSNR2_BGS_Z</th><th>TSNR2_GPBBACKUP_Z</th><th>TSNR2_QSO_Z</th><th>TSNR2_LRG_Z</th><th>TSNR2_GPBDARK</th><th>TSNR2_ELG</th><th>TSNR2_GPBBRIGHT</th><th>TSNR2_LYA</th><th>TSNR2_BGS</th><th>TSNR2_GPBBACKUP</th><th>TSNR2_QSO</th><th>TSNR2_LRG</th><th>ZCAT_NSPEC</th><th>ZCAT_PRIMARY</th><th>SV3_DESI_TARGET</th><th>SV3_BGS_TARGET</th><th>SV3_MWS_TARGET</th><th>SV3_SCND_TARGET</th><th>SV1_DESI_TARGET</th><th>SV1_BGS_TARGET</th><th>SV1_MWS_TARGET</th><th>SV1_SCND_TARGET</th><th>CMX_TARGET</th><th>SCND_TARGET</th></tr></thead>
<thead><tr><th>int64</th><th>int32</th><th>int32</th><th>float64</th><th>float64</th><th>int64</th><th>float64</th><th>float64</th><th>int64</th><th>str6</th><th>str20</th><th>int64</th><th>float64</th><th>int32</th><th>float64</th><th>float64</th><th>float32</th><th>float32</th><th>float32</th><th>int64</th><th>uint8</th><th>str3</th><th>float64</th><th>int32</th><th>int16</th><th>str8</th><th>int32</th><th>int32</th><th>str4</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>int16</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>int64</th><th>str2</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>str1</th><th>int64</th><th>int64</th><th>int64</th><th>int64</th><th>int64</th><th>int64</th><th>int64</th><th>int64</th><th>int64</th><th>float64</th><th>float64</th><th>int16</th><th>float32</th><th>int16</th><th>int16</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float64</th><th>float32</th><th>float64</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>int16</th><th>bool</th><th>int64</th><th>int64</th><th>int64</th><th>int64</th><th>int64</th><th>int64</th><th>int64</th><th>int64</th><th>int64</th><th>int64</th></tr></thead>
<tr><td>616093992420902061</td><td>5611</td><td>5611</td><td>1.5324390612456673</td><td>0.00015800737090060045</td><td>1</td><td>7091.521788910031</td><td>0.00040601428665243313 .. 0.0</td><td>7902</td><td>QSO</td><td></td><td>4</td><td>11.633388578891754</td><td>0</td><td>169.94283602332484</td><td>49.623484693601476</td><td>0.0</td><td>0.0</td><td>0.0</td><td>4294967296</td><td>4</td><td>SKY</td><td>0.9433153925937842</td><td>63</td><td>9011</td><td>1698p495</td><td>582458</td><td>1197</td><td></td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.008927048</td><td>0.05416071</td><td>-0.06267472</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0</td><td></td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td></td><td>-1</td><td>-1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>4294967296</td><td>0</td><td>0</td><td>169.94283602332484</td><td>49.623484693601476</td><td>1</td><td>139.7048</td><td>1</td><td>1</td><td>-0.002</td><td>0.002</td><td>-0.013</td><td>0.013</td><td>169.94281964488044</td><td>0.0</td><td>49.623533378359404</td><td>0.0</td><td>0.7933969</td><td>2.3602135</td><td>0.0008920009</td><td>0.4870128</td><td>0.43644428</td><td>5.9276433</td><td>4.303284</td><td>0.022763047</td><td>0.011190208</td><td>196.85054</td><td>0.48135555</td><td>39.82439</td><td>0.0005672656</td><td>43.711002</td><td>326.72336</td><td>0.13496664</td><td>0.67481124</td><td>3.6294017e-07</td><td>1.3604566</td><td>7.418995e-08</td><td>0.0</td><td>83.72673</td><td>6.325105e-07</td><td>0.30303684</td><td>0.71666944</td><td>199.21075</td><td>1.8427042</td><td>40.3114</td><td>0.43701154</td><td>133.36537</td><td>331.02664</td><td>0.46076652</td><td>1.4026709</td><td>1</td><td>True</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td></tr>
<tr><td>616093996346769480</td><td>5611</td><td>5611</td><td>1.4118471585727892</td><td>0.000230938417354251</td><td>5</td><td>7214.78203972429</td><td>0.00020190209113513124 .. 0.0</td><td>7915</td><td>QSO</td><td></td><td>4</td><td>0.12016404420137405</td><td>0</td><td>169.96167901158464</td><td>49.63702509483948</td><td>0.0</td><td>0.0</td><td>0.0</td><td>4294967296</td><td>4</td><td>SKY</td><td>0.9869011941878603</td><td>63</td><td>9011</td><td>1697p497</td><td>583394</td><td>72</td><td></td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>-0.019832069</td><td>-0.036477424</td><td>-0.059049565</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0</td><td></td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td></td><td>-1</td><td>-1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>4294967296</td><td>0</td><td>0</td><td>169.96167901158464</td><td>49.63702509483948</td><td>1</td><td>139.7048</td><td>1</td><td>1</td><td>-0.003</td><td>0.003</td><td>-0.012</td><td>0.012</td><td>169.96165710461366</td><td>0.0</td><td>49.63707021066727</td><td>0.0</td><td>0.7932764</td><td>3.1076524</td><td>0.0011778958</td><td>0.63929856</td><td>0.5835937</td><td>7.794383</td><td>5.6322613</td><td>0.030148352</td><td>0.014715493</td><td>244.28203</td><td>0.59473354</td><td>49.202515</td><td>0.0007010528</td><td>54.01597</td><td>400.76266</td><td>0.16749799</td><td>0.834453</td><td>4.270615e-07</td><td>1.6360322</td><td>8.705572e-08</td><td>0.0</td><td>98.86646</td><td>7.4100444e-07</td><td>0.36287627</td><td>0.8527024</td><td>247.38968</td><td>2.2319436</td><td>49.841812</td><td>0.5842948</td><td>160.67682</td><td>406.39493</td><td>0.5605226</td><td>1.7018709</td><td>1</td><td>True</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td></tr>
<tr><td>2305843015094114261</td><td>5611</td><td>5611</td><td>2.9628891901018706e-05</td><td>3.812523334597231e-06</td><td>0</td><td>12705.285074676825</td><td>-13216.656936150071 .. 0.0</td><td>7906</td><td>STAR</td><td>M</td><td>5</td><td>17785.598721649993</td><td>0</td><td>169.9297144636733</td><td>49.64493920163109</td><td>-42.193623</td><td>-0.12222115</td><td>2015.5</td><td>2305843009213693952</td><td>1</td><td>TGT</td><td>0.07751337399126856</td><td>568</td><td>-1</td><td></td><td>583394</td><td>-1</td><td>GPSF</td><td>0.015491215</td><td>-99.0</td><td>-99.0</td><td>-99.0</td><td>0.0</td><td>0.0</td><td>-99.0</td><td>-99.0</td><td>-99.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>789818514070421632</td><td>G2</td><td>16.907398</td><td>18.213533</td><td>15.78182</td><td>3.0894492</td><td>G</td><td>8</td><td>1</td><td>2305843009213693952</td><td>0</td><td>2305843009213693952</td><td>0</td><td>0</td><td>0</td><td>0</td><td>169.9297144636733</td><td>49.64493920163109</td><td>1</td><td>139.7048</td><td>1</td><td>1</td><td>-0.002</td><td>0.002</td><td>-0.008</td><td>0.008</td><td>169.92969982719904</td><td>0.0</td><td>49.644969289957515</td><td>0.0</td><td>0.7928637</td><td>2.8348725</td><td>0.001036787</td><td>0.5812703</td><td>0.49597213</td><td>6.871747</td><td>5.174559</td><td>0.02668175</td><td>0.013191072</td><td>227.9924</td><td>0.561436</td><td>45.770233</td><td>0.0006481258</td><td>50.205357</td><td>376.67896</td><td>0.1583804</td><td>0.78410405</td><td>4.1011822e-07</td><td>1.5572447</td><td>8.331422e-08</td><td>0.0</td><td>93.981895</td><td>7.161341e-07</td><td>0.34805506</td><td>0.81127036</td><td>230.82727</td><td>2.1197174</td><td>46.351505</td><td>0.49662024</td><td>151.05899</td><td>381.85352</td><td>0.53311723</td><td>1.6085654</td><td>1</td><td>True</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td></tr>
<tr><td>2305843015094114259</td><td>5611</td><td>5611</td><td>0.00018049483015363774</td><td>8.785155971487304e-06</td><td>0</td><td>8350.513439494453</td><td>-3969.747333324784 .. 0.0</td><td>7912</td><td>STAR</td><td>M</td><td>5</td><td>3266.9647149562024</td><td>0</td><td>169.9634740025819</td><td>49.67011395238441</td><td>-41.692318</td><td>-8.550553</td><td>2015.5</td><td>2305843009213693952</td><td>1</td><td>TGT</td><td>0.9067837874882217</td><td>568</td><td>-1</td><td></td><td>583394</td><td>-1</td><td>GPSF</td><td>0.015475582</td><td>-99.0</td><td>-99.0</td><td>-99.0</td><td>0.0</td><td>0.0</td><td>-99.0</td><td>-99.0</td><td>-99.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>789818720228851456</td><td>G2</td><td>18.391855</td><td>19.584547</td><td>17.203827</td><td>2.7993507</td><td>G</td><td>7</td><td>1</td><td>2305843009213693952</td><td>0</td><td>4611686018427387904</td><td>0</td><td>0</td><td>0</td><td>0</td><td>169.9634740025819</td><td>49.67011395238441</td><td>1</td><td>139.7048</td><td>1</td><td>1</td><td>-0.003</td><td>0.003</td><td>-0.009</td><td>0.009</td><td>169.96345313198785</td><td>0.0</td><td>49.67014797702819</td><td>0.0</td><td>0.79290897</td><td>2.927817</td><td>0.0010755257</td><td>0.60011137</td><td>0.5184494</td><td>7.117638</td><td>5.340162</td><td>0.027695453</td><td>0.013643058</td><td>233.7924</td><td>0.57696813</td><td>46.91244</td><td>0.0006638387</td><td>51.43387</td><td>385.79974</td><td>0.16258451</td><td>0.80516547</td><td>4.2276986e-07</td><td>1.6103836</td><td>8.583723e-08</td><td>0.0</td><td>96.68688</td><td>7.3704257e-07</td><td>0.3593661</td><td>0.8372437</td><td>236.72021</td><td>2.1884272</td><td>47.512554</td><td>0.51911324</td><td>155.23839</td><td>391.1399</td><td>0.549646</td><td>1.6560522</td><td>1</td><td>True</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td></tr>
<tr><td>2305843015098305819</td><td>5612</td><td>5612</td><td>-9.153748044020623e-05</td><td>1.1456058252671455e-05</td><td>0</td><td>8229.175911848926</td><td>9221.186777179724 .. 0.0</td><td>7921</td><td>STAR</td><td>G</td><td>5</td><td>1198.8152144064788</td><td>0</td><td>171.4521877421757</td><td>49.503026247748664</td><td>0.96577805</td><td>-0.43101287</td><td>2015.5</td><td>2305843009213693952</td><td>1</td><td>TGT</td><td>0.9743120206490967</td><td>568</td><td>-1</td><td></td><td>582462</td><td>-1</td><td>GPSF</td><td>0.012459293</td><td>-99.0</td><td>-99.0</td><td>-99.0</td><td>0.0</td><td>0.0</td><td>-99.0</td><td>-99.0</td><td>-99.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>789955884305450496</td><td>G2</td><td>17.215637</td><td>17.48034</td><td>16.778713</td><td>0.3028206</td><td>G</td><td>8</td><td>1</td><td>2305843009213693952</td><td>0</td><td>2305843009213693952</td><td>0</td><td>0</td><td>0</td><td>0</td><td>171.4521877421757</td><td>49.503026247748664</td><td>1</td><td>139.7048</td><td>1</td><td>1</td><td>-0.001</td><td>0.001</td><td>-0.003</td><td>0.003</td><td>171.4521816389004</td><td>0.0</td><td>49.503037438193644</td><td>0.0</td><td>0.7923077</td><td>2.4682126</td><td>0.0013917657</td><td>0.50132346</td><td>1.0653914</td><td>8.483982</td><td>4.503935</td><td>0.03693317</td><td>0.013542183</td><td>225.87901</td><td>0.48856232</td><td>45.19003</td><td>0.0008841664</td><td>52.436287</td><td>383.72012</td><td>0.15118614</td><td>0.70890987</td><td>3.0388176e-07</td><td>1.8087848</td><td>6.144743e-08</td><td>0.0</td><td>91.27699</td><td>5.4308833e-07</td><td>0.38225442</td><td>0.85123986</td><td>228.34723</td><td>2.298739</td><td>45.691353</td><td>1.0662756</td><td>152.19727</td><td>388.22406</td><td>0.5703737</td><td>1.5736918</td><td>1</td><td>True</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td></tr>
<tr><td>616093992441873040</td><td>5612</td><td>5612</td><td>1.4994798482456368</td><td>0.00013146039148774144</td><td>2053</td><td>6725.964683115482</td><td>-476.0183967132606 .. -47.7889182734013</td><td>7909</td><td>GALAXY</td><td></td><td>10</td><td>0.09657600522041321</td><td>1024</td><td>171.6413175641367</td><td>49.50413530583062</td><td>0.0</td><td>0.0</td><td>0.0</td><td>4294967296</td><td>4</td><td>SKY</td><td>0.9874304395982679</td><td>63</td><td>9011</td><td>1717p495</td><td>582463</td><td>656</td><td></td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0070314202</td><td>-0.06831814</td><td>-0.12136616</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0</td><td></td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td></td><td>-1</td><td>-1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>4294967296</td><td>0</td><td>0</td><td>171.6413175641367</td><td>49.50413530583062</td><td>1</td><td>139.7048</td><td>1</td><td>1</td><td>-0.035</td><td>0.035</td><td>0.003</td><td>0.003</td><td>171.64110310924636</td><td>0.0</td><td>49.50412355623945</td><td>0.0</td><td>0.7967332</td><td>1.7721186</td><td>0.0010117285</td><td>0.38068527</td><td>0.69767004</td><td>6.5149264</td><td>2.9403958</td><td>0.02398039</td><td>0.01008215</td><td>165.37267</td><td>0.3555587</td><td>35.095394</td><td>0.00057925354</td><td>40.734966</td><td>258.75363</td><td>0.09861642</td><td>0.5317616</td><td>2.1865927e-07</td><td>1.3056929</td><td>4.6799464e-08</td><td>0.0</td><td>70.18367</td><td>3.5647903e-07</td><td>0.24707851</td><td>0.6336545</td><td>167.14479</td><td>1.6622634</td><td>35.476078</td><td>0.6982493</td><td>117.43356</td><td>261.69403</td><td>0.36967534</td><td>1.1754982</td><td>1</td><td>True</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td></tr>
<tr><td>2305843015098305930</td><td>5612</td><td>5612</td><td>-4.7067493605894956e-05</td><td>2.020237445900416e-06</td><td>0</td><td>14205.571482520158</td><td>39633.630649387924 .. 0.0</td><td>7916</td><td>STAR</td><td>K</td><td>5</td><td>21296.415092121097</td><td>0</td><td>171.29642389914463</td><td>49.49892932262815</td><td>-10.557034</td><td>-20.238108</td><td>2015.5</td><td>2305843009213693952</td><td>1</td><td>TGT</td><td>0.1996474071002604</td><td>568</td><td>-1</td><td></td><td>582462</td><td>-1</td><td>GPSF</td><td>0.013379046</td><td>-99.0</td><td>-99.0</td><td>-99.0</td><td>0.0</td><td>0.0</td><td>-99.0</td><td>-99.0</td><td>-99.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>789956571500220672</td><td>G2</td><td>15.547872</td><td>16.147432</td><td>14.798199</td><td>1.5068805</td><td>G</td><td>9</td><td>1</td><td>2305843009213693952</td><td>0</td><td>1152921504606846976</td><td>0</td><td>0</td><td>0</td><td>0</td><td>171.29642389914463</td><td>49.49892932262815</td><td>1</td><td>139.7048</td><td>1</td><td>1</td><td>-0.005</td><td>0.005</td><td>-0.003</td><td>0.003</td><td>171.29639317893327</td><td>0.0</td><td>49.4989405338948</td><td>0.0</td><td>0.79243493</td><td>2.5167806</td><td>0.0014249525</td><td>0.51165843</td><td>1.0903715</td><td>8.683438</td><td>4.5804987</td><td>0.03773854</td><td>0.013861531</td><td>228.69937</td><td>0.49574354</td><td>45.794926</td><td>0.00088724535</td><td>53.077404</td><td>387.44916</td><td>0.15285365</td><td>0.719283</td><td>3.0826266e-07</td><td>1.8386524</td><td>6.238669e-08</td><td>0.0</td><td>92.93129</td><td>5.493213e-07</td><td>0.38786665</td><td>0.86591643</td><td>231.21616</td><td>2.335821</td><td>46.306583</td><td>1.0912588</td><td>154.69214</td><td>392.02966</td><td>0.57845885</td><td>1.599061</td><td>1</td><td>True</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td></tr>
<tr><td>2305843015098307288</td><td>5612</td><td>5612</td><td>-1.3664007772099054e-05</td><td>3.9143383944587955e-06</td><td>0</td><td>12181.071298493323</td><td>-12435.3300485166 .. 0.0</td><td>7895</td><td>STAR</td><td>M</td><td>5</td><td>17209.169397668662</td><td>0</td><td>171.68736941679126</td><td>49.54470455872235</td><td>-69.38101</td><td>-57.846035</td><td>2015.5</td><td>2305843009213693952</td><td>1</td><td>TGT</td><td>0.68223051364576</td><td>568</td><td>-1</td><td></td><td>582463</td><td>-1</td><td>GPSF</td><td>0.010762735</td><td>-99.0</td><td>-99.0</td><td>-99.0</td><td>0.0</td><td>0.0</td><td>-99.0</td><td>-99.0</td><td>-99.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>789953994519832576</td><td>G2</td><td>17.116034</td><td>18.450846</td><td>15.916414</td><td>4.117403</td><td>G</td><td>8</td><td>1</td><td>2305843009213693952</td><td>0</td><td>2305843009213693952</td><td>0</td><td>0</td><td>0</td><td>0</td><td>171.68736941679126</td><td>49.54470455872235</td><td>1</td><td>139.7048</td><td>1</td><td>1</td><td>-0.002</td><td>0.002</td><td>0.004</td><td>0.004</td><td>171.68735687684952</td><td>0.0</td><td>49.54468951604813</td><td>0.0</td><td>0.7923018</td><td>2.7017956</td><td>0.0015122459</td><td>0.5481989</td><td>1.1526091</td><td>9.215066</td><td>4.920465</td><td>0.040109936</td><td>0.014792004</td><td>245.78578</td><td>0.5312204</td><td>49.087536</td><td>0.00095886947</td><td>56.978313</td><td>415.4498</td><td>0.16458322</td><td>0.77099156</td><td>3.2435977e-07</td><td>1.9311537</td><td>6.5515614e-08</td><td>0.0</td><td>96.9108</td><td>5.7839753e-07</td><td>0.40789756</td><td>0.90728706</td><td>248.48758</td><td>2.4638863</td><td>49.635735</td><td>1.153568</td><td>163.10419</td><td>420.37027</td><td>0.61259073</td><td>1.6930707</td><td>1</td><td>True</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td></tr>
<tr><td>2305843015098307317</td><td>5612</td><td>5612</td><td>8.187076272145688e-06</td><td>2.663950196736789e-06</td><td>0</td><td>12131.354403508387</td><td>48041.18003633886 .. 0.0</td><td>7914</td><td>STAR</td><td>G</td><td>5</td><td>15919.053546092564</td><td>0</td><td>171.4626471249043</td><td>49.571728763263934</td><td>-6.487234</td><td>-8.979567</td><td>2015.5</td><td>2305843009213693952</td><td>1</td><td>TGT</td><td>0.10845467407798826</td><td>568</td><td>-1</td><td></td><td>582462</td><td>-1</td><td>GPSF</td><td>0.012409228</td><td>-99.0</td><td>-99.0</td><td>-99.0</td><td>0.0</td><td>0.0</td><td>-99.0</td><td>-99.0</td><td>-99.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>789957224335242112</td><td>G2</td><td>15.332517</td><td>15.647106</td><td>14.842437</td><td>0.5882234</td><td>G</td><td>9</td><td>1</td><td>2305843009213693952</td><td>0</td><td>1152921504606846976</td><td>0</td><td>0</td><td>0</td><td>0</td><td>171.4626471249043</td><td>49.571728763263934</td><td>1</td><td>139.7048</td><td>1</td><td>1</td><td>-0.001</td><td>0.001</td><td>-0.004</td><td>0.004</td><td>171.4626410100553</td><td>0.0</td><td>49.57174383066378</td><td>0.0</td><td>0.7922408</td><td>2.6429763</td><td>0.0014779883</td><td>0.5357533</td><td>1.126929</td><td>8.998085</td><td>4.822091</td><td>0.039301917</td><td>0.014445578</td><td>238.74109</td><td>0.51726633</td><td>47.650562</td><td>0.00092558406</td><td>55.271835</td><td>404.8719</td><td>0.16029207</td><td>0.7498533</td><td>3.1255968e-07</td><td>1.8509881</td><td>6.309495e-08</td><td>0.0</td><td>93.48485</td><td>5.593211e-07</td><td>0.3925457</td><td>0.8718289</td><td>241.38406</td><td>2.3697324</td><td>48.186317</td><td>1.1278546</td><td>157.75476</td><td>409.69397</td><td>0.5921397</td><td>1.6361278</td><td>1</td><td>True</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td></tr>
<tr><td>2305843015098307300</td><td>5612</td><td>5612</td><td>-9.071880405177143e-05</td><td>2.2039071072745236e-06</td><td>0</td><td>17519.42825152598</td><td>97226.82042910179 .. 0.0</td><td>7898</td><td>STAR</td><td>K</td><td>5</td><td>26668.015562488148</td><td>0</td><td>171.522010239778</td><td>49.539644541926656</td><td>11.656263</td><td>-17.613537</td><td>2015.5</td><td>2305843009213693952</td><td>1</td><td>TGT</td><td>0.2606906966336341</td><td>568</td><td>-1</td><td></td><td>582462</td><td>-1</td><td>GPSF</td><td>0.012510003</td><td>-99.0</td><td>-99.0</td><td>-99.0</td><td>0.0</td><td>0.0</td><td>-99.0</td><td>-99.0</td><td>-99.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>789957430493671168</td><td>G2</td><td>14.528907</td><td>14.91449</td><td>13.972879</td><td>1.9625176</td><td>G</td><td>9</td><td>1</td><td>2305843009213693952</td><td>0</td><td>1152921504606846976</td><td>0</td><td>0</td><td>0</td><td>0</td><td>171.522010239778</td><td>49.539644541926656</td><td>1</td><td>139.7048</td><td>1</td><td>1</td><td>-0.006</td><td>0.006</td><td>-0.005</td><td>0.005</td><td>171.52197351199567</td><td>0.0</td><td>49.53966324619804</td><td>0.0</td><td>0.792537</td><td>2.5464804</td><td>0.0013881064</td><td>0.5177623</td><td>1.0218446</td><td>8.553396</td><td>4.629491</td><td>0.03663849</td><td>0.013815446</td><td>228.41553</td><td>0.49797174</td><td>45.754147</td><td>0.00086462416</td><td>52.960907</td><td>386.95078</td><td>0.15297045</td><td>0.7216756</td><td>3.0635482e-07</td><td>1.8440832</td><td>6.201617e-08</td><td>0.0</td><td>92.53374</td><td>5.4565237e-07</td><td>0.38729942</td><td>0.8651746</td><td>230.962</td><td>2.343443</td><td>46.271908</td><td>1.0227093</td><td>154.04803</td><td>391.58026</td><td>0.57690835</td><td>1.6006656</td><td>1</td><td>True</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td></tr>
<tr><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
<tr><td>39633365241956461</td><td>6013</td><td>6013</td><td>0.25004622998506426</td><td>3.4640081372709474e-05</td><td>0</td><td>7147.220601081848</td><td>219.29282171260482 .. 0.8472593965517198</td><td>7904</td><td>GALAXY</td><td></td><td>10</td><td>251.51401315256953</td><td>0</td><td>177.91967886982204</td><td>58.327432477875576</td><td>0.0</td><td>0.0</td><td>2015.5</td><td>1152921504606846976</td><td>1</td><td>TGT</td><td>0.5207872915855238</td><td>516</td><td>9011</td><td>1778p582</td><td>612290</td><td>2157</td><td>DEV</td><td>0.025377158</td><td>6.872843</td><td>19.698236</td><td>37.61195</td><td>40.676037</td><td>28.184603</td><td>152.8685</td><td>44.26062</td><td>22.551344</td><td>2.63561</td><td>0.80983424</td><td>2.3678026</td><td>6.7863526</td><td>12.95791</td><td>2.3678026</td><td>6.7863526</td><td>12.95791</td><td>0</td><td>4.0</td><td>1.1563034</td><td>-0.21848726</td><td>-0.055595778</td><td>0</td><td></td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>N</td><td>2100</td><td>1</td><td>1152921504606846976</td><td>514</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>177.91967886982204</td><td>58.327432477875576</td><td>1</td><td>420.9098</td><td>1</td><td>1</td><td>-0.009</td><td>0.009</td><td>0.0</td><td>0.0</td><td>177.9196113146371</td><td>0.0</td><td>58.327431869341</td><td>0.0</td><td>0.73560476</td><td>17.465498</td><td>0.01072765</td><td>3.5135489</td><td>6.828758</td><td>64.236115</td><td>29.126053</td><td>0.28108925</td><td>0.124423124</td><td>1236.3129</td><td>2.8179736</td><td>244.6197</td><td>0.0062433775</td><td>275.36774</td><td>1915.3945</td><td>0.8270845</td><td>4.0092807</td><td>3.5262087e-06</td><td>15.447745</td><td>6.982991e-07</td><td>0.0</td><td>726.4123</td><td>5.48854e-06</td><td>3.1494138</td><td>7.0912127</td><td>1253.7783</td><td>18.276447</td><td>248.13326</td><td>6.835001</td><td>1066.0161</td><td>1944.5206</td><td>4.2575874</td><td>11.224916</td><td>1</td><td>True</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td></tr>
<tr><td>39633365241957327</td><td>6013</td><td>6013</td><td>0.29923634505142904</td><td>6.411668820672214e-05</td><td>0</td><td>7104.332029595971</td><td>92.0788520797761 .. 8.064678423271925</td><td>7864</td><td>GALAXY</td><td></td><td>10</td><td>122.14909310638905</td><td>0</td><td>178.01296014251585</td><td>58.233671759041414</td><td>0.0</td><td>0.0</td><td>2015.5</td><td>1152921504606846976</td><td>1</td><td>TGT</td><td>0.8505689818070435</td><td>516</td><td>9011</td><td>1778p582</td><td>612290</td><td>3023</td><td>SER</td><td>0.02326719</td><td>8.744032</td><td>27.030731</td><td>53.964527</td><td>86.07819</td><td>70.80018</td><td>93.34585</td><td>26.105515</td><td>12.9111395</td><td>1.9885625</td><td>0.6636225</td><td>1.4386163</td><td>4.4472446</td><td>8.878541</td><td>1.4711175</td><td>4.5068827</td><td>8.901539</td><td>0</td><td>1.3526341</td><td>2.1857274</td><td>-0.17959079</td><td>-0.2231453</td><td>0</td><td></td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>N</td><td>2100</td><td>1</td><td>1152921504606846976</td><td>514</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>178.01296014251585</td><td>58.233671759041414</td><td>1</td><td>420.9098</td><td>1</td><td>1</td><td>-0.007</td><td>0.007</td><td>0.0</td><td>0.0</td><td>178.01290741677596</td><td>0.0</td><td>58.23367140098279</td><td>0.0</td><td>0.7249042</td><td>16.565802</td><td>0.010134473</td><td>3.322152</td><td>6.519705</td><td>60.646935</td><td>27.922485</td><td>0.26719707</td><td>0.11709482</td><td>1185.9607</td><td>2.7286441</td><td>234.0317</td><td>0.0059893443</td><td>262.1917</td><td>1860.4282</td><td>0.8018252</td><td>3.8626204</td><td>3.5502044e-06</td><td>14.995379</td><td>7.0051715e-07</td><td>0.0</td><td>708.8078</td><td>5.572211e-06</td><td>3.0878954</td><td>6.9067483</td><td>1202.5265</td><td>17.734158</td><td>237.35385</td><td>6.5256944</td><td>1031.6465</td><td>1888.3507</td><td>4.1569176</td><td>10.886463</td><td>1</td><td>True</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td></tr>
<tr><td>616094117545378425</td><td>6013</td><td>6013</td><td>1.633433845547278</td><td>8.350731274008447e-05</td><td>5</td><td>6911.096681892872</td><td>-534.3286311138368 .. -33.988643782019814</td><td>7900</td><td>GALAXY</td><td></td><td>10</td><td>5.256517753005028</td><td>0</td><td>178.06789837672608</td><td>58.258914069362476</td><td>0.0</td><td>0.0</td><td>0.0</td><td>4294967296</td><td>4</td><td>SKY</td><td>0.9622348355735991</td><td>63</td><td>9011</td><td>1778p582</td><td>612290</td><td>633</td><td></td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>-0.0045087235</td><td>-0.001621763</td><td>0.0025148385</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0</td><td></td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td></td><td>-1</td><td>-1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>4294967296</td><td>0</td><td>0</td><td>178.06789837672608</td><td>58.258914069362476</td><td>1</td><td>420.9098</td><td>1</td><td>1</td><td>-0.004</td><td>0.004</td><td>0.001</td><td>0.001</td><td>178.0678682304532</td><td>0.0</td><td>58.25891020375297</td><td>0.0</td><td>0.7924418</td><td>18.922798</td><td>0.012129488</td><td>3.790143</td><td>8.103597</td><td>71.858376</td><td>31.89588</td><td>0.32096797</td><td>0.13709714</td><td>1338.0668</td><td>3.0361586</td><td>263.57407</td><td>0.0067711184</td><td>295.5335</td><td>2094.2188</td><td>0.89963865</td><td>4.3121686</td><td>3.8650005e-06</td><td>16.556255</td><td>7.614365e-07</td><td>0.0</td><td>777.57196</td><td>6.0584516e-06</td><td>3.4091833</td><td>7.5968184</td><td>1356.9896</td><td>19.604544</td><td>267.3642</td><td>8.110368</td><td>1144.9639</td><td>2126.1147</td><td>4.62979</td><td>12.046084</td><td>1</td><td>True</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td></tr>
<tr><td>39633365241956284</td><td>6013</td><td>6013</td><td>0.33866003723057986</td><td>1.0743872643550098e-05</td><td>0</td><td>9336.443044304848</td><td>-864.7354373685928 .. -143.17542540017487</td><td>7882</td><td>GALAXY</td><td></td><td>10</td><td>4807.291874483693</td><td>0</td><td>177.90170652335553</td><td>58.297648143990926</td><td>0.0</td><td>0.0</td><td>2015.5</td><td>1152921504606846976</td><td>1</td><td>TGT</td><td>0.4603449391815422</td><td>519</td><td>9011</td><td>1778p582</td><td>612290</td><td>1980</td><td>DEV</td><td>0.024074892</td><td>2.58818</td><td>9.684639</td><td>20.48293</td><td>31.595804</td><td>31.162659</td><td>259.60184</td><td>79.954605</td><td>50.896637</td><td>3.2299082</td><td>0.90416026</td><td>1.2404065</td><td>4.641443</td><td>9.816612</td><td>1.2404065</td><td>4.641443</td><td>9.816612</td><td>0</td><td>4.0</td><td>0.69541645</td><td>0.4006547</td><td>-0.0467836</td><td>0</td><td></td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>N</td><td>2000</td><td>1</td><td>1152921504606846976</td><td>257</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>177.90170652335553</td><td>58.297648143990926</td><td>1</td><td>420.9098</td><td>1</td><td>1</td><td>-0.007</td><td>0.007</td><td>0.006</td><td>0.006</td><td>177.90165372431255</td><td>0.0</td><td>58.29762564645565</td><td>0.0</td><td>0.75474757</td><td>15.905952</td><td>0.00994939</td><td>3.2011817</td><td>6.3660626</td><td>59.46956</td><td>26.685738</td><td>0.260907</td><td>0.11512282</td><td>1146.6636</td><td>2.6249971</td><td>227.06905</td><td>0.0057732663</td><td>253.86412</td><td>1789.8021</td><td>0.7685811</td><td>3.7282414</td><td>3.38807e-06</td><td>14.983832</td><td>6.708386e-07</td><td>0.0</td><td>701.6087</td><td>5.290269e-06</td><td>3.052361</td><td>6.857854</td><td>1162.5696</td><td>17.618778</td><td>230.27023</td><td>6.3718357</td><td>1014.9424</td><td>1816.4879</td><td>4.081849</td><td>10.701218</td><td>1</td><td>True</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td></tr>
<tr><td>39633365241957363</td><td>6013</td><td>6013</td><td>0.5143961903252277</td><td>9.132162758078593e-05</td><td>0</td><td>7301.54689759016</td><td>518.7270214329949 .. 18.007455886267635</td><td>7896</td><td>GALAXY</td><td></td><td>10</td><td>125.40195485204458</td><td>0</td><td>178.01647090733843</td><td>58.33858556091349</td><td>0.0</td><td>0.0</td><td>2015.5</td><td>1152921504606846976</td><td>1</td><td>TGT</td><td>0.2524009032661023</td><td>519</td><td>9011</td><td>1778p582</td><td>612290</td><td>3059</td><td>DEV</td><td>0.029326499</td><td>3.3817282</td><td>11.615948</td><td>27.954557</td><td>57.627365</td><td>38.777466</td><td>106.47783</td><td>40.3256</td><td>31.856411</td><td>2.5536206</td><td>0.85161835</td><td>1.3209742</td><td>4.537434</td><td>10.919639</td><td>1.3224176</td><td>4.5405917</td><td>10.922338</td><td>0</td><td>4.0</td><td>0.78460646</td><td>0.024017273</td><td>0.028956374</td><td>0</td><td></td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>N</td><td>2000</td><td>1</td><td>1152921504606846976</td><td>257</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>178.01647090733843</td><td>58.33858556091349</td><td>1</td><td>420.9098</td><td>1</td><td>1</td><td>-0.01</td><td>0.01</td><td>-0.002</td><td>0.002</td><td>178.0163958504652</td><td>0.0</td><td>58.338592327019455</td><td>0.0</td><td>0.7485425</td><td>14.94997</td><td>0.009348906</td><td>3.0156715</td><td>5.9390492</td><td>56.061104</td><td>25.080227</td><td>0.24443887</td><td>0.10841598</td><td>1075.0309</td><td>2.4634438</td><td>213.49242</td><td>0.0053826403</td><td>238.93777</td><td>1681.7593</td><td>0.7187806</td><td>3.501831</td><td>3.1931522e-06</td><td>14.017997</td><td>6.341206e-07</td><td>0.0</td><td>657.38654</td><td>4.9992004e-06</td><td>2.8453622</td><td>6.4176908</td><td>1089.9808</td><td>16.49079</td><td>216.50809</td><td>5.944432</td><td>952.3854</td><td>1706.8395</td><td>3.8085816</td><td>10.027938</td><td>1</td><td>True</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td></tr>
<tr><td>39633365241957469</td><td>6013</td><td>6013</td><td>0.31613469941419087</td><td>7.631133947109325e-05</td><td>0</td><td>7454.033276565373</td><td>398.9613053665811 .. -0.353706600985458</td><td>7896</td><td>GALAXY</td><td></td><td>10</td><td>195.02904924750328</td><td>0</td><td>178.02761946091448</td><td>58.31644769067127</td><td>0.0</td><td>0.0</td><td>2015.5</td><td>1152921504606846976</td><td>1</td><td>TGT</td><td>0.7149079054822364</td><td>516</td><td>9011</td><td>1778p582</td><td>612290</td><td>3165</td><td>DEV</td><td>0.028268065</td><td>6.2711134</td><td>25.704601</td><td>52.268925</td><td>66.94721</td><td>46.34028</td><td>80.11901</td><td>25.910759</td><td>21.82172</td><td>2.1570017</td><td>0.7435695</td><td>2.0604255</td><td>8.4454565</td><td>17.173382</td><td>2.0834978</td><td>8.485226</td><td>17.227318</td><td>0</td><td>4.0</td><td>1.0648221</td><td>0.037470546</td><td>-0.022615695</td><td>0</td><td></td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>N</td><td>2100</td><td>1</td><td>1152921504606846976</td><td>514</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>178.02761946091448</td><td>58.31644769067127</td><td>1</td><td>420.9098</td><td>1</td><td>1</td><td>-0.01</td><td>0.01</td><td>-0.003</td><td>0.003</td><td>178.02754437185502</td><td>0.0</td><td>58.316458164069964</td><td>0.0</td><td>0.7372002</td><td>15.657929</td><td>0.0098034805</td><td>3.1557755</td><td>6.2844357</td><td>58.633327</td><td>26.227812</td><td>0.25651154</td><td>0.113439605</td><td>1139.59</td><td>2.614507</td><td>225.98457</td><td>0.0056386557</td><td>252.38261</td><td>1775.6364</td><td>0.7618699</td><td>3.7126472</td><td>3.4037691e-06</td><td>14.701818</td><td>6.7477487e-07</td><td>0.0</td><td>693.4621</td><td>5.3015715e-06</td><td>2.9958122</td><td>6.7502675</td><td>1155.2479</td><td>17.32613</td><td>229.14035</td><td>6.2900743</td><td>1004.478</td><td>1801.8641</td><td>4.0141935</td><td>10.576354</td><td>1</td><td>True</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td></tr>
<tr><td>39633365241956347</td><td>6013</td><td>6013</td><td>0.35569677281126383</td><td>8.098014249121414e-05</td><td>0</td><td>7188.658670205623</td><td>415.68008955268874 .. -18.271361476317722</td><td>7886</td><td>GALAXY</td><td></td><td>10</td><td>255.08269446715713</td><td>0</td><td>177.9082840707307</td><td>58.2389189365421</td><td>0.0</td><td>0.0</td><td>2015.5</td><td>1152921504606847233</td><td>1</td><td>TGT</td><td>0.4534677191459664</td><td>517</td><td>9011</td><td>1778p582</td><td>612290</td><td>2043</td><td>SER</td><td>0.021861132</td><td>6.2471185</td><td>31.017164</td><td>68.59942</td><td>112.355606</td><td>75.30025</td><td>158.29771</td><td>28.124582</td><td>19.698257</td><td>1.6565682</td><td>0.62351745</td><td>1.5766715</td><td>7.828229</td><td>17.313381</td><td>1.5766715</td><td>7.828229</td><td>17.313381</td><td>0</td><td>3.2551847</td><td>1.58709</td><td>0.11802119</td><td>-0.015740044</td><td>0</td><td></td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>N</td><td>2100</td><td>1</td><td>1152921504606847233</td><td>514</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>177.9082840707307</td><td>58.2389189365421</td><td>1</td><td>420.9098</td><td>1</td><td>1</td><td>-0.013</td><td>0.013</td><td>0.001</td><td>0.001</td><td>177.9081862419601</td><td>0.0</td><td>58.23891436015811</td><td>0.0</td><td>0.7199077</td><td>16.071577</td><td>0.0103326095</td><td>3.2418313</td><td>6.6936836</td><td>61.55877</td><td>26.86742</td><td>0.26989946</td><td>0.11850463</td><td>1153.9398</td><td>2.6344311</td><td>229.13242</td><td>0.0056826314</td><td>256.33847</td><td>1798.2148</td><td>0.76810306</td><td>3.7480443</td><td>3.4683767e-06</td><td>15.012391</td><td>6.883628e-07</td><td>0.0</td><td>704.49774</td><td>5.398496e-06</td><td>3.045844</td><td>6.88144</td><td>1170.0114</td><td>17.657154</td><td>232.37425</td><td>6.699366</td><td>1022.395</td><td>1825.0823</td><td>4.0838466</td><td>10.747989</td><td>1</td><td>True</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td></tr>
<tr><td>39633365241956215</td><td>6013</td><td>6013</td><td>0.14662384155326322</td><td>2.3908997454831108e-05</td><td>0</td><td>8596.399122469127</td><td>1868.6202392422365 .. -59.50317164496653</td><td>7887</td><td>GALAXY</td><td></td><td>10</td><td>2590.241582274437</td><td>0</td><td>177.89490141739614</td><td>58.35852792661476</td><td>0.4301085</td><td>0.76327</td><td>2015.5</td><td>1152921504606846976</td><td>1</td><td>TGT</td><td>0.3891468072281197</td><td>516</td><td>9011</td><td>1778p582</td><td>612290</td><td>1911</td><td>SER</td><td>0.024757447</td><td>93.01637</td><td>189.10666</td><td>305.42184</td><td>195.54918</td><td>124.573006</td><td>36.425156</td><td>10.2032795</td><td>7.571413</td><td>0.9629097</td><td>0.42319807</td><td>20.41846</td><td>41.511692</td><td>67.044586</td><td>21.181093</td><td>43.684822</td><td>70.91316</td><td>0</td><td>2.6966248</td><td>1.6853192</td><td>0.055675875</td><td>-0.041010804</td><td>846272113961956352</td><td>G2</td><td>18.834333</td><td>18.055773</td><td>16.869219</td><td>1.4527224</td><td>N</td><td>2100</td><td>1</td><td>1152921504606846976</td><td>514</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>177.89490141739614</td><td>58.35852792661476</td><td>1</td><td>420.9098</td><td>1</td><td>1</td><td>-0.011</td><td>0.011</td><td>-0.001</td><td>0.001</td><td>177.89481899342462</td><td>0.0</td><td>58.35853076330482</td><td>0.0</td><td>0.71744066</td><td>14.758839</td><td>0.009321137</td><td>2.9826393</td><td>5.82738</td><td>56.02078</td><td>24.704695</td><td>0.24292165</td><td>0.10888249</td><td>1097.4716</td><td>2.5046775</td><td>218.10712</td><td>0.005422326</td><td>243.30623</td><td>1704.586</td><td>0.7283551</td><td>3.5634873</td><td>3.2480973e-06</td><td>13.850452</td><td>6.4576835e-07</td><td>0.0</td><td>656.0079</td><td>5.05841e-06</td><td>2.8165429</td><td>6.3857346</td><td>1112.2303</td><td>16.36445</td><td>221.08975</td><td>5.8328023</td><td>955.33484</td><td>1729.2908</td><td>3.7878196</td><td>10.0581045</td><td>1</td><td>True</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td></tr>
<tr><td>39633365241957523</td><td>6013</td><td>6013</td><td>0.3382530298537338</td><td>2.0956507875619487e-05</td><td>0</td><td>7136.956998258829</td><td>399.2937157077019 .. 0.25989794965738083</td><td>7892</td><td>GALAXY</td><td></td><td>10</td><td>422.62999925017357</td><td>0</td><td>178.03421923603548</td><td>58.23109857422396</td><td>0.0</td><td>0.0</td><td>2015.5</td><td>1152921504606846976</td><td>1</td><td>TGT</td><td>0.7296236131406513</td><td>519</td><td>9011</td><td>1778p582</td><td>612290</td><td>3219</td><td>REX</td><td>0.023558747</td><td>3.292489</td><td>10.391624</td><td>20.506294</td><td>28.92958</td><td>21.6674</td><td>275.367</td><td>88.74141</td><td>55.492905</td><td>3.6641133</td><td>1.0646143</td><td>1.6489635</td><td>5.2043934</td><td>10.2700815</td><td>1.6501979</td><td>5.206145</td><td>10.271273</td><td>0</td><td>1.0</td><td>0.5178733</td><td>0.0</td><td>0.0</td><td>0</td><td></td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>N</td><td>2100</td><td>1</td><td>1152921504606846976</td><td>265</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>178.03421923603548</td><td>58.23109857422396</td><td>1</td><td>420.9098</td><td>1</td><td>1</td><td>-0.01</td><td>0.01</td><td>-0.001</td><td>0.001</td><td>178.03414391761183</td><td>0.0</td><td>58.231101803872725</td><td>0.0</td><td>0.7655765</td><td>15.598676</td><td>0.009909298</td><td>3.1361609</td><td>6.411245</td><td>59.002144</td><td>26.280376</td><td>0.26061594</td><td>0.1141834</td><td>1151.0795</td><td>2.6326618</td><td>227.72981</td><td>0.005706269</td><td>254.44676</td><td>1804.6569</td><td>0.77194405</td><td>3.7357605</td><td>3.5133014e-06</td><td>14.846387</td><td>6.9460725e-07</td><td>0.0</td><td>699.8964</td><td>5.4973552e-06</td><td>3.043437</td><td>6.8309984</td><td>1166.6781</td><td>17.488958</td><td>230.86597</td><td>6.416951</td><td>1013.34534</td><td>1830.9373</td><td>4.075997</td><td>10.680943</td><td>1</td><td>True</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td></tr>
<tr><td>39633365241956510</td><td>6013</td><td>6013</td><td>9.435878078810976e-05</td><td>4.368909273144145e-06</td><td>0</td><td>12531.322434358151</td><td>19554.60669899227 .. 0.0</td><td>7913</td><td>STAR</td><td>G</td><td>5</td><td>15585.968943333211</td><td>0</td><td>177.92464209948585</td><td>58.213943163967556</td><td>-8.031009</td><td>-12.605257</td><td>2015.5</td><td>2305843052163366912</td><td>3</td><td>TGT</td><td>0.3117450453320081</td><td>519</td><td>9011</td><td>1778p582</td><td>612290</td><td>2206</td><td>PSF</td><td>0.021217449</td><td>234.12079</td><td>306.48447</td><td>318.3136</td><td>74.15675</td><td>39.963947</td><td>23.813667</td><td>10.026296</td><td>14.550704</td><td>2.5543613</td><td>0.93732435</td><td>181.45976</td><td>237.54662</td><td>246.715</td><td>181.45976</td><td>237.54662</td><td>246.715</td><td>16</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>846254865374029312</td><td>G2</td><td>16.3488</td><td>16.6024</td><td>15.932671</td><td>0.21672922</td><td>N</td><td>1500</td><td>1</td><td>2305843052163366912</td><td>0</td><td>768</td><td>0</td><td>0</td><td>0</td><td>0</td><td>177.92464209948585</td><td>58.213943163967556</td><td>1</td><td>420.9098</td><td>1</td><td>1</td><td>-0.009</td><td>0.009</td><td>0.001</td><td>0.001</td><td>177.92457431086444</td><td>0.0</td><td>58.213938870323204</td><td>0.0</td><td>0.789</td><td>15.601844</td><td>0.009918077</td><td>3.134498</td><td>6.4335527</td><td>59.01516</td><td>26.355299</td><td>0.2611379</td><td>0.11411495</td><td>1154.7517</td><td>2.6235275</td><td>228.14279</td><td>0.0058516376</td><td>254.91153</td><td>1810.1</td><td>0.7726246</td><td>3.7255268</td><td>3.3842687e-06</td><td>14.50318</td><td>6.688177e-07</td><td>0.0</td><td>680.041</td><td>5.3167305e-06</td><td>2.9706128</td><td>6.6481895</td><td>1170.3535</td><td>17.136625</td><td>231.27728</td><td>6.4394045</td><td>993.9677</td><td>1836.4553</td><td>4.0043755</td><td>10.487831</td><td>1</td><td>True</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td></tr>
</table></div>




```python
def match_tns_desi(tns, zcat, rad=1.0, makeplot=False):
    """Write me."""
    tns_coord = SkyCoord(tns['RA'], tns['DEC'], frame='icrs', unit=['hour', 'degree'])
    desi_coord = SkyCoord(zcat['TARGET_RA'], zcat['TARGET_DEC'], frame='icrs', unit='degree')
    
    minsep = rad * u.arcsec
    idx, d2d, _ = desi_coord.match_to_catalog_3d(tns_coord)  
    indx_desi = np.where(d2d < minsep)[0]
    indx_tns = idx[indx_desi]
    print('{} matches'.format(len(indx_desi)))
    
    if makeplot:
        pngfile = 'tns-desi-matches.png'
        fig = plt.figure(figsize=(10, 8))
        plt.scatter(tns_coord.ra.value[indx_tns], tns_coord.dec.value[indx_tns], marker='s', s=100, color='k')
        plt.scatter(desi_coord.ra.value[indx_desi], desi_coord.dec.value[indx_desi], marker='x', s=120, color='orange')
        
        print('Writing {}'.format(pngfile))
        plt.savefig(pngfile)
    
    return indx_desi, indx_tns
```


```python
%time indx_desi, indx_tns = match_tns_desi(tns, zcat, rad=1.0, makeplot=True)
```


    ---------------------------------------------------------------------------

    ValueError                                Traceback (most recent call last)

    <timed exec> in <module>


    /tmp/ipykernel_51062/325338605.py in match_tns_desi(tns, zcat, rad, makeplot)
          1 def match_tns_desi(tns, zcat, rad=1.0, makeplot=False):
          2     """Write me."""
    ----> 3     tns_coord = SkyCoord(tns['RA'], tns['DEC'], frame='icrs', unit=['hour', 'degree'])
          4     desi_coord = SkyCoord(zcat['TARGET_RA'], zcat['TARGET_DEC'], frame='icrs', unit='degree')
          5 


    /global/common/software/desi/cori/desiconda/20211217-2.0.0/conda/lib/python3.9/site-packages/astropy/coordinates/sky_coordinate.py in __init__(self, copy, *args, **kwargs)
        329             # creating the internal self._sky_coord_frame object
        330             args = list(args)  # Make it mutable
    --> 331             skycoord_kwargs, components, info = _parse_coordinate_data(
        332                 frame_cls(**frame_kwargs), args, kwargs)
        333 


    /global/common/software/desi/cori/desiconda/20211217-2.0.0/conda/lib/python3.9/site-packages/astropy/coordinates/sky_coordinate_parsers.py in _parse_coordinate_data(frame, args, kwargs)
        294                                                                   repr_attr_names, units):
        295                 attr_class = frame.representation_type.attr_classes[repr_attr_name]
    --> 296                 _components[frame_attr_name] = attr_class(arg, unit=unit)
        297 
        298         else:


    /global/common/software/desi/cori/desiconda/20211217-2.0.0/conda/lib/python3.9/site-packages/astropy/coordinates/angles.py in __new__(cls, angle, unit, **kwargs)
        562             raise TypeError("A Latitude angle cannot be created from a Longitude angle")
        563         self = super().__new__(cls, angle, unit=unit, **kwargs)
    --> 564         self._validate_angles()
        565         return self
        566 


    /global/common/software/desi/cori/desiconda/20211217-2.0.0/conda/lib/python3.9/site-packages/astropy/coordinates/angles.py in _validate_angles(self, angles)
        583                               np.any(angles.value > upper))
        584         if invalid_angles:
    --> 585             raise ValueError('Latitude angle(s) must be within -90 deg <= angle <= 90 deg, '
        586                              'got {}'.format(angles.to(u.degree)))
        587 


    ValueError: Latitude angle(s) must be within -90 deg <= angle <= 90 deg, got [203.6868625  278.40055    219.75627917 ... 256.22062917 282.54704167
     122.30360833] deg



```python
indx_I = np.array(['SN I' in typ for typ in tns['OBJ. TYPE']])
set(tns[indx_I]['OBJ. TYPE'])
```


```python
indx_II = np.array(['SN II' in typ for typ in tns['OBJ. TYPE']])
set(tns[indx_II]['OBJ. TYPE'])
```




    {'SN II', 'SN II-pec', 'SN IIL', 'SN IIP', 'SN IIb', 'SN IIn', 'SN IIn-pec'}




```python
def read_matches(makeplot=False):
    from astropy.time import Time

    tnsfile = 'tns_all_thru_20220615.csv'
    print('Reading {}'.format(tnsfile))
    tns = Table.read(tnsfile, format='csv')

    # Convert discovery date to MJD.
    mjd = Time(tns['Discovery Date (UT)'], format='iso', scale='utc').mjd
    tns['Discovery Date (MJD)'] = mjd
    
    if makeplot:
        pngfile = 'tns-desi-matchsky.png'
        fig = plt.figure(figsize=(10, 8))
        ax = init_sky(projection = 'aitoff', ra_center = 90, galactic_plane_color = 'Red')
        ax.scatter(ax.projection_ra(tns_coord.ra.value[indx_tns]), ax.projection_dec(tns_coord.dec.value[indx_tns]), 
                   marker='.', color = 'Black', s=120)
        ax.scatter(ax.projection_ra(desi_coord.ra.value[indx_desi]), ax.projection_dec(desi_coord.dec.value[indx_desi]), 
                   marker='x', color = 'orange', s=25)
        print('Writing {}'.format(pngfile))
        plt.savefig(pngfile)
    
    return tns
```


```python
tns_desi_match = read_matches(makeplot=True)
```

    Reading tns_all_thru_20220615.csv



    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    /tmp/ipykernel_51062/2723765960.py in <module>
    ----> 1 tns_desi_match = read_matches(makeplot=True)
    

    /tmp/ipykernel_51062/287048112.py in read_matches(makeplot)
         14         fig = plt.figure(figsize=(10, 8))
         15         ax = init_sky(projection = 'aitoff', ra_center = 90, galactic_plane_color = 'Red')
    ---> 16         ax.scatter(ax.projection_ra(tns_coord.ra.value[indx_tns]), ax.projection_dec(tns_coord.dec.value[indx_tns]), 
         17                    marker='.', color = 'Black', s=120)
         18         ax.scatter(ax.projection_ra(desi_coord.ra.value[indx_desi]), ax.projection_dec(desi_coord.dec.value[indx_desi]), 


    NameError: name 'tns_coord' is not defined



    <Figure size 720x576 with 0 Axes>



    
![png](output_13_3.png)
    



```python
def read_matches2(makeplot=False):
    from astropy.time import Time

    tnsfile = 'tns_all_thru_20220615.csv'
    print('Reading {}'.format(tnsfile))
    tns = Table.read(tnsfile, format='csv')

    # Convert discovery date to MJD.
    mjd = Time(tns['Discovery Date (UT)'], format='iso', scale='utc').mjd
    tns['Discovery Date (MJD)'] = mjd
    
    if makeplot:
        pngfile = 'tns-desi-matchsky2.png'
        fig = plt.figure(figsize=(10, 8))
        ax = init_sky(projection = 'aitoff', ra_center = 90, galactic_plane_color = 'Red')
        ax.scatter(ax.projection_ra(tns_coord.ra.degree), ax.projection_dec(tns_coord.dec.degree), 
                   marker='.', color = 'Black', s=5)
        ax.scatter(ax.projection_ra(desi_coord.ra.value[indx_desi]), ax.projection_dec(desi_coord.dec.value[indx_desi]), 
                   marker='x', color = 'orange', s=60)
        print('Writing {}'.format(pngfile))
        plt.savefig(pngfile)
    
    return tns
```


```python
tns_desi_match2 = read_matches2(makeplot=True)
```


```python
def read_matches3(makeplot=False):
    from astropy.time import Time

    tnsfile = 'tns_all_thru_20220615.csv'
    print('Reading {}'.format(tnsfile))
    tns = Table.read(tnsfile, format='csv')

    # Convert discovery date to MJD.
    mjd = Time(tns['Discovery Date (UT)'], format='iso', scale='utc').mjd
    tns['Discovery Date (MJD)'] = mjd
    
    if makeplot:
        pngfile = 'tns-desi-matchsky3.png'
        fig = plt.figure(figsize=(10, 8))
        ax = init_sky(projection = 'aitoff', ra_center = 90, galactic_plane_color = 'Red')
        ax.scatter(ax.projection_ra(tns_coord.ra.degree), ax.projection_dec(tns_coord.dec.degree), 
                   marker='.', color = 'Black', s=15)
        ax.scatter(ax.projection_ra(tns_coord.ra.value[indx_tns]), ax.projection_dec(tns_coord.dec.value[indx_tns]), 
                   marker='.', color = 'Blue', s=120)
        ax.scatter(ax.projection_ra(desi_coord.ra.value[indx_desi]), ax.projection_dec(desi_coord.dec.value[indx_desi]), 
                   marker='x', color = 'Orange', s=25)
        print('Writing {}'.format(pngfile))
        plt.savefig(pngfile)
    
    return tns
```


```python
tns_desi_match3 = read_matches3(makeplot=True)
```


```python
combined = hstack([tns_all[indx_tns],zcat[indx_desi]])
combined
```


```python
tns_all[indx_tns]['Redshift'].data
```


```python
zcat[indx_desi]['Z'].data
```


```python
import matplotlib.pyplot as plt

plt.figure(figsize=(10,10))

plt.scatter(tns_coord.ra.value[indx_tns],tns_coord.dec.value[indx_tns], marker= '.', color = 'red')
plt.scatter(desi_coord.ra.value[indx_desi], desi_coord.dec.value[indx_desi], marker= 'x', color = 'black')
```


```python
# idx -- indices of TNS matches in zcat
idx, d2d, d3d = match_coordinates_sky(tns_coord, desi_coord)
zcat[idx]
```


```python
stop
```


```python
ax = init_sky(projection = 'aitoff', ra_center = 90, galactic_plane_color = 'Red')
p = ax.scatter(ax.projection_ra(matches.ra.degree), ax.projection_dec(matches.dec.degree), marker='.', color = 'Black')
```


```python

```


```python
stop
```


```python
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import SphericalCircle
from datetime import datetime, timedelta
from scipy.ndimage import gaussian_filter1d
from glob import glob
import psycopg2
from tqdm.notebook import tqdm_notebook
import matplotlib as mpl
```


```python
from desispec.io import read_spectra
from desispec.coaddition import coadd_cameras
from desispec.interpolation import resample_flux
from desispec.resolution import Resolution

import redrock.templates
```


```python
rrtemplates = dict()
for filename in redrock.templates.find_templates():
    t = redrock.templates.Template(filename)
    rrtemplates[(t.template_type, t.sub_type)] = t
```


```python
mpl.rc('font', size=14)
mpl.rc('figure', max_open_warning = 0)
```


```python
tns_all = Table.read('tns_all_thru_20220615.csv', format='csv')

# Convert discovery date to MJD.
mjd = Time(tns_all['Discovery Date (UT)'], format='iso', scale='utc').mjd
tns_all['Discovery Date (MJD)'] = mjd

# Convert coordinates to decimal degree format. Do not insert into the TNS table.
coords = SkyCoord(tns_all['RA'], tns_all['DEC'], frame='icrs', unit=['hour', 'degree'])
```


```python
tns_all
```


```python

```
