'''
Chemevol - Python package to read in a star formation history file,
input galaxy parameters and run chemical evolution to determine the evolution
of gas, metals and dust in galaxies.

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

'''
- lifetime: the lifetime of stars of a given mass
            1st col: initial mass of star Msolar

            2nd col: lifetime in Gyrs, metallicity (Z) < 0.008
            3rd col: lifetime in Gyrs, Z >= 0.008

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

t_lifetime['mass'].unit = 'solMass'
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
mass_yields =np.array([(0.9, 0, 0, 0, 1.08e-5, 0, 6.83e-3, 0, 6.83e-3),
             (1.0, 0, 0, 0, 8.54e-4, 0, 1.12e-4, 0, 1.61e-4),
             (1.3, 0, 3.09e-3, 0, 1.99e-3, 0, 1.70e-30, 0, 1.70e-3),
             (1.5, 0, 3.53e-3, 0, 2.57e-3, 0, 3.06e-3, 0, 3.06e-3),
             (1.7, 0, 3.77e-3, 0, 3.99e-3, 0, 3.51e-3, 0, 3.51e-3),
             (2.0, 0, 5.29e-3, 0, 5.86e-3, 0, 5.43e-3, 0, 5.43e-3),
             (2.5, 0, 5.59e-3, 0, 6.82e-3, 0, 6.58e-3, 0, 6.58e-3),
             (3, 0, 5.35e-3, 0, 8.43e-3, 0, 7.98e-3, 0, 7.98e-3),
             (4, 0, 6.56e-3, 0, 5.87e-3, 0, 5.44e-3, 0, 5.44e-3),
             (5, 0, 7.72e-3, 0, 7.07e-3, 0, 6.59e-3, 0, 6.59e-3),
             (7, 0, 9.61e-3, 0, 8.88e-3, 0, 8.35e-3, 0, 8.35e-3),
             (8, 0, 8.78e-3, 0, 9.57e-3, 0, 1.02e-2, 0, 1.02e-2),
             (9, 0.27, 0, 0.27, 0, 0.173, 0, 0.173, 0),
             (12, 0.83, 0, 0.83, 0, 0.686, 0, 0.686, 0),
             (15, 1.53, 0, 1.53, 0, 1.32, 0, 1.32, 0),
             (20, 2.93, 0, 2.93, 0, 2.73, 0, 2.73, 0),
             (25, 4.45, 0, 4.45, 0, 4.48, 0, 4.48, 0),
             (40, 4.71, 0, 4.71, 0, 1.61, 6.4, 1.61, 6.4),
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
dust_mass_sn: dust mass returned by supernovae
            1st column: initial mass of star
            2nd column: dust mass returned

            From Todinin & Ferrara 2001 (MNRAS 325 276)
            As TF01 only have SN masses for m_i > 12Msun
            (ie no 9Msun entry), we add dust for 9Msun progenitor
            by assuming similar dust/metals
            ratio for 12-20Msun stars
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
    Take a 2D lifetime table and return nearest neighbour based on either 2 or 3 col
    '''
    col1 = lookup[column]
    idx = (np.abs(col1-value)).argmin()
    return lookup[idx]

def lookup_taum(mass, colnum):
#    '''
#    Take a lifetime list
#    '''
    col1 = lifetime[:,0]
    idx = (np.abs(col1-mass)).argmin()
    #print mass, col1, lifetime[idx][colnum]
    return lifetime[idx][colnum]
