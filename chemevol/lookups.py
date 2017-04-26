'''
Chemevol - Python package to read in a star formation history file,
input galaxy parameters and run a chemical evolution model to determine the evolution
of gas, metals and dust in galaxies.

Running this script will produce a results data file (filename.dat) with file name given by user

The code is based on Morgan & Edmunds 2003 (MNRAS, 343, 427)
and described in detail in Rowlands et al 2014 (MNRAS, 441, 1040).

If you use this code, please do cite the above papers.

Copyright (C) 2015 Haley Gomez, Edward Gomez and Simon Schofield, Cardiff University and LCOGT
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
'''

'''
This file lists the different tables required for lifetimes of stars of
different initial mass, and metal mass yields from
supernovae or from stellar winds at different metallicities.

It also includes a function to lookup nearest neighbour quantities
'''

from astropy.table import Table
import numpy as np
from bisect import bisect_left

'''
- lifetime: the lifetime of stars of a given mass
            1st col: initial mass of star Msolar

            2nd col: lifetime in Gyrs, metallicity (Z) = 0.001
            3rd col: lifetime in Gyrs, Z = 0.02

            Values from Schaller et al 1992 (A & AS 96 269)
            this table allows us to call lifetime star but also
            mass of star that corresponds to age of system
'''
lifetime =  np.array([(0.8, 15.0, 26.0),
            (0.9, 9.5, 15.0),
            (1.0, 6.3, 10.0),
            (1.5, 1.8, 2.7),
            (2.0, 0.86, 1.1),
            (3.0, 0.29, 0.35),
            (4.0, 0.14, 0.16),
            (5.0, 0.088, 0.094),
            (7.0, 0.045, 0.043),
            (9.0, 0.029, 0.026),
            (12.0, 0.018, 0.016),
            (20.0, 0.0094, 0.0081),
            (40.0, 0.0049, 0.0043),
            (60.0, 0.0037, 0.0034),
            (85.0, 0.0031, 0.0028),
            (120.0, 0.0028, 0.0026)])

t_lifetime = Table(rows=lifetime, names=('mass','lifetime_low_metals','lifetime_high_metals'),meta={'name': 'Lifetime'})

t_lifetime['mass'].unit = 'solMass' # write units to each column
t_lifetime['lifetime_low_metals'].unit = 'Gyr'
t_lifetime['lifetime_high_metals'].unit = 'Gyr'




'''
- mass_yields: ejected yield (all heavy elements) in Msolar
            1st col: initial mass of star Msolar

            2nd col: ejected yield Msolar from supernovae Z=0.001
            3rd col: ejected yield Msolar from winds Z=0.001

            4th col: ejected yield Msolar from supernovae Z=0.004
            5th col: ejected yield Msolar from winds Z=0.004

            6th col: ejected yield Msolar from supernovae Z=0.008
            7th col: ejected yield Msolar from winds Z=0.008

            8th col: ejected yield Msolar from supernovae Z=0.02
            9th col: ejected yield Msolar from winds Z=0.02

            M >= 9 Msolar: mp_Z from Maeder 1992 (A & A 264 105)
            Z = 0.001 (used for Z < 0.008), Z = 0.02 (used for 0.008 =< Z < inf.)
            M < 9 Msolar: mp_Z from van den Hoek &
            Groenewegen 1997 (A & AS 123 305) Z=0.001,0.004,0.008,0.02
'''
mass_yields =np.array([
      (     0.9      ,       0      ,      0.0      ,       0      , 9.72e-06      ,       0      , 6.147e-05      ,      0      ,      0.0       )      ,
      (     1.0      ,       0      ,      0.0      ,       0      , 0.000854      ,       0      ,  0.000112      ,      0      ,  0.00161       )      ,
      (     1.3      ,       0      , 0.004017      ,       0      , 0.002587      ,       0      ,   0.00221      ,      0      , 0.003939       )      ,
      (     1.5      ,       0      , 0.005295      ,       0      , 0.003855      ,       0      ,   0.00459      ,      0      ,  0.00312       )      ,
      (     1.7      ,       0      , 0.006409      ,       0      , 0.006783      ,       0      ,  0.005967      ,      0      , 0.004029       )      ,
      (     2.0      ,       0      ,  0.01058      ,       0      ,  0.01172      ,       0      ,   0.01086      ,      0      ,  0.00788       )      ,
      (     2.5      ,       0      , 0.013975      ,       0      ,  0.01705      ,       0      ,   0.01645      ,      0      ,    0.014       )      ,
      (     3.0      ,       0      ,  0.01605      ,       0      ,  0.02529      ,       0      ,   0.02394      ,      0      ,  0.02076       )      ,
      (     4.0      ,       0      ,  0.02624      ,       0      ,  0.02348      ,       0      ,   0.02176      ,      0      ,  0.02496       )      ,
      (     5.0      ,       0      ,   0.0386      ,       0      ,  0.03535      ,       0      ,   0.03295      ,      0      ,   0.0314       )      ,
      (     7.0      ,       0      ,  0.06727      ,       0      ,  0.06216      ,       0      ,   0.05845      ,      0      ,  0.05418       )      ,
      (     8.0      ,       0      ,  0.07024      ,       0      ,  0.07656      ,       0      ,    0.0816      ,      0      ,  0.06728       )      ,
             (9, 0.27, 0, 0.27, 0, 0.173, 0, 0.173, 0),
             (12, 0.83, 0, 0.83, 0, 0.686, 0, 0.686, 0),
             (15, 1.53, 0, 1.53, 0, 1.32, 0, 1.32, 0),
             (20, 2.93, 0, 2.93, 0, 2.73, 0, 2.73, 0),
             (25, 4.45, 0, 4.45, 0, 4.48, 0, 4.48, 0),
             (40, 9.71, 0, 9.71, 0, 1.61, 6.4, 1.61, 6.4),
             (60, 17.1, 0, 17.1, 0, 1.16, 8.69, 1.16, 8.69),
             (85, 26.7, 0, 26.7, 0, 1.56, 17.75, 1.56, 17.75),
             (120, 41.6, 0, 41.6, 0, 0.72, 9.39, 0.72, 9.39)])
yield_names = ['mass','yields_sn_001','yields_winds_001',
             'yields_sn_004','yields_winds_004',
             'yields_sn_008','yields_winds_008',
             'yields_sn_02','yields_winds_02']

t_yields = Table(rows=mass_yields, names=yield_names,
                                          meta={'name': 'Mass Yields'})

t_yields['mass'].unit = 'solMass'
t_yields['yields_sn_001'].unit = 'solMass'
t_yields['yields_winds_001'].unit = 'solMass'
t_yields['yields_sn_004'].unit = 'solMass'
t_yields['yields_winds_004'].unit = 'solMass'
t_yields['yields_sn_008'].unit = 'solMass'
t_yields['yields_winds_008'].unit = 'solMass'
t_yields['yields_sn_02'].unit = 'solMass'
t_yields['yields_winds_02'].unit = 'solMass'

'''
- oxymass_yields: ejected oxygen yield in Msolar
            1st col: initial mass of star Msolar

            2nd col: ejected yield Msolar from supernovae Z=0.001
            3rd col: ejected yield Msolar from winds Z=0.001

            4th col: ejected yield Msolar from supernovae Z=0.004
            5th col: ejected yield Msolar from winds Z=0.004

            6th col: ejected yield Msolar from supernovae Z=0.008
            7th col: ejected yield Msolar from winds Z=0.008

            8th col: ejected yield Msolar from supernovae Z=0.02
            9th col: ejected yield Msolar from winds Z=0.02

            M >= 9 Msolar: mp_Z from Maeder 1992 (A & A 264 105)
            Z = 0.001 (used for Z < 0.008), Z = 0.02 (used for 0.008 =< Z < inf.)
            M < 9 Msolar: mp_Z from van den Hoek &
            Groenewegen 1997 (A & AS 123 305) Z=0.001,0.004,0.008,0.02
'''
oxymass_yields =np.array([(0.9, 0, -1.773e-06, 0, -6.498e-07, 0, 2.565e-05, 0, -3.483e-05),
                (1, 0, -2.23e-06, 0, 6.36e-05, 0, 5.36e-05, 0, 0.000981),
                (1.3, 0, 0.0003237, 0, 0.0001872, 0, 0.0001807, 0, 0.002431),
                (1.5, 0, 0.000426, 0, 0.0003105, 0, 0.000324, 0, 0.0005565),
                (1.7, 0, 0.0005168, 0, 0.0005253, 0, 0.0004233, 0, 0.0003026),
                (2, 0, 0.0008, 0, 0.000722, 0, 0.000508, 0, 0.0001554),
                (2.5, 0, 0.0010725, 0, 0.00099, 0, 0.0005875, 0, -0.00010475),
                (3, 0, 0.001251, 0, 0.001536, 0, 0.000915, 0, -4.41e-06),
                (4, 0, 0.001724, 0, 0.001104, 0, 0.0002288, 0, -0.000864),
                (5, 0, 0.00206, 0, 0.001285, 0, 0.00033, 0, -0.001455),
                (7, 0, 0.0004403, 0, -0.001778, 0, -0.004424, 0, -0.00896),
                (8, 0, 0.000768, 0, -0.002568, 0, -0.007936, 0, -0.0128),
                (9, 0.004, 0, 0.004, 0, 0, 0, 0, 0),
                (12, 0.15, 0, 0.15, 0, 0.11, 0, 0.11, 0),
                (15, 0.46, 0, 0.46, 0, 0.41, 0, 0.41, 0),
                (20, 1.27, 0, 1.27, 0, 1.27, 0, 1.27, 0),
                (25, 2.40, 0, 2.40, 0, 2.60, -0.03, 2.60, -0.03),
                (40, 6.80, 0, 6.80, 0, 0.62, 1.46, 0.62, 1.46),
                (60, 14.20, 0, 14.20, 0, 0.40, 1.03, 0.40, 1.03),
                (85, 22.60, 0, 22.60, 0, 0.59, 3.37, 0.59, 3.37),
                (120,35.30, 0,35.30, 0, 0.18, -0.13, 0.18, -0.13)])


oxyt_yields = Table(rows=oxymass_yields, names=yield_names,
                                          meta={'name': 'Oxygen Mass Yields'})

oxyt_yields['mass'].unit = 'solMass'
oxyt_yields['yields_sn_001'].unit = 'solMass'
oxyt_yields['yields_winds_001'].unit = 'solMass'
oxyt_yields['yields_sn_004'].unit = 'solMass'
oxyt_yields['yields_winds_004'].unit = 'solMass'
oxyt_yields['yields_sn_008'].unit = 'solMass'
oxyt_yields['yields_winds_008'].unit = 'solMass'
oxyt_yields['yields_sn_02'].unit = 'solMass'
oxyt_yields['yields_winds_02'].unit = 'solMass'


'''
dust_mass_sn: dust mass returned by supernovae from Todini & Ferrara 2001 (MNRAS 325 276)
            1st column: initial mass of star
            2nd column: dust mass returned

            As TF01 only have SN masses for m_i > 12Msun
            (ie no 9Msun entry), we add dust for 8.5 + 9Msun progenitor
            by assuming no dust from 8.5 Msun, and similar dust/metals
            ratio for 12-20Msun for 9Msun
'''
dust_mass_sn = ((8.5,0),
                (9, 0.17),
                (12, 0.2),
                (15, 0.5),
                (20, 0.5),
                (22, 0.8),
                (25, 1.0),
                (30, 1.0),
                (35, 0.6),
                (40, 0.4))

t_dustmass_sn = Table(rows=dust_mass_sn, names=('mass','dustmass'), meta={'name': 'Dust Mass SN'})
t_dustmass_sn['mass'].unit = 'solMass'
t_dustmass_sn['dustmass'].unit = 'solMass'

def find_nearest(lookup,value):
    '''
    Take a 2D array and return pair of nearest neighbour values based on first column
    '''
    # indices = zip(*lookup)
    # pos = bisect_left(indices[0], value)
    # if pos == 0:
    #     return lookup[0]
    # if pos == len(indices[0]):
    #     return lookup[-1]
    # before = indices[0][pos - 1]
    # after = indices[0][pos]
    # if after - value < value - before:
    #    return lookup[pos]
    # else:
    #    return lookup[pos - 1]
    col1 = lookup[:,0]
    idx = (np.abs(col1-value)).argmin()
    return lookup[idx]

def find_nearest_col(lookup,value,colnum):
    '''
    Take a 2D array and return pair of nearest neighbour values based on first column
    '''
    col1 = lookup[:,0]
    idx = (np.abs(col1-value)).argmin()
    return lookup[idx][colnum]

def lookup_fn(lookup, column, value):
    '''
    return nearest neighbour value from function based on column input
    '''
    col1 = lookup[column]
    idx = (np.abs(col1-value)).argmin()
    return lookup[idx]

def lookup_taum(mass, colnum):
    '''
    Looks up nearest neighbour for lifetime of star of mass m from lifetime table
    reads in the relevant colnum for correct metallicity
    '''
    col1 = lifetime[:,0]
    idx = (np.abs(col1-mass)).argmin()
    #print mass, col1, lifetime[idx][colnum]
    return lifetime[idx][colnum]
