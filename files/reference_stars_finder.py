#given an input catalogue, generates the list of reference targets and appends
#to a new file appended with 'reference_targets' """
import numpy as np
from astropy.io import fits, ascii
import astropy.units as u
from astropy.coordinates import SkyCoord
from scipy.spatial import cKDTree
import os
import tqdm
refmag = 'NRS_F110W'
fullcat_path = 'Final_input_catalog_10nov2022.csv'
fullcat = ascii.read(fullcat_path)
brightest_in_3p2 = []
deltamag = []
data = np.dstack([fullcat['RA_final'], fullcat['DEC_final']])[0]
brighter = ( (fullcat[refmag] < 20.0) & (fullcat[refmag] > 19.3) )
skycoords = SkyCoord(ra=fullcat['RA_final'][brighter]*u.deg, dec=fullcat['DEC_final'][brighter]*u.deg)
for i in tqdm.tqdm(range(len(skycoords))):
    #make a mask to remove the star being checked
    checker = np.ones(len(skycoords), dtype=bool)
    checker[i] = False
    #make the list of stars to check (removing the target)
    checkcoords = skycoords[checker]
    seps = skycoords[i].separation(checkcoords)
    #get all stars in the list within 3.2 arcsec
    nearby = seps < 3.2*u.arcsec
    within2 = seps < 2*u.arcsec
    #check if the target object is the brightest
    magnitudes = fullcat[refmag][brighter][checker][nearby].data
    magnitudes2 = fullcat[refmag][brighter][checker][within2].data
    #if brightest in 3.2 arcsec, add to list
    if (magnitudes > fullcat[refmag][brighter][i]).all():
        brightest_in_3p2.append(fullcat['ID'][brighter][i])
        if len(magnitudes2) > 0:
            #compute magnitude difference between reference star and (next) brightest in 3.2 (for testing)
            deltamags = magnitudes2 - fullcat[refmag][brighter][i]
            deltamag.append(np.min(deltamags))
        else:
            deltamag.append(np.inf)
#make a list of all stars with no neighbours brighter than reference mag = 23.3 in 0.3 arcsec:
mask = fullcat[refmag] < 24.
data = np.dstack([fullcat['RA_final'], fullcat['DEC_final']])[0]
tree = cKDTree(data)
bright_in_0p3 = tree.query_ball_point(data, 0.3/3600)
print(bright_in_0p3)
close_neighbour = np.ones(len(bright_in_0p3), dtype=bool)
for i in range(len(bright_in_0p3)):
    if len(bright_in_0p3[i]) == 1:
        close_neighbour[i] = 0
brightest_in_3p2_mask = np.zeros(len(fullcat), dtype=bool)
brightest_in_3p2_mask[brightest_in_3p2] = 1
#finally, make the list of magnitude differences
all_deltamag = np.ones(len(fullcat))*np.nan
all_deltamag[brightest_in_3p2_mask] = deltamag
refmask = brightest_in_3p2_mask & ~close_neighbour & (fullcat[refmag] < 22.) & (fullcat[refmag] > 19.5) & (all_deltamag > 1.5)
reference = np.zeros(len(fullcat), dtype=bool)
reference[refmask] = 1
rec = np.recarray(len(fullcat), dtype=[('RA_final', float),('DEC_final', float),('ref_mag', float), ('Reference', int)])
rec['RA_final'] = fullcat['RA_final']
rec['DEC_final'] = fullcat['DEC_final']
rec['Reference'] = reference
ascii.write(rec, os.path.splitext(os.path.basename(fullcat_path))[0]+'_reference_targets.dat',  overwrite=True, format='csv')











Message r.p.schiavon



