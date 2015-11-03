'''
Chemevol - Python package to read in a star formation history file,
input galaxy parameters and run chemical evolution to determine the evolution
of gas, metals and dust in galaxies.

The code is based on Morgan & Edmunds 2003 (MNRAS, 343, 427)
and described in detail in Rowlands et al 2014 (MNRAS, 441, 1040).

If you use this code, please do cite the above papers.

Copyright (C) 2015 Haley Gomez and Edward Gomez, Cardiff University and LCOGT
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
This lists the different tables required for lifetimes of stars of
different initial mass, and metal mass yields from
supernovae or from stellar winds at different metallicities
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
lifetime =  ((0.8, 15.0, 26.0),
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
            (120.0, 0.0028, 0.0026))

t_lifetime = Table(rows=lifetime, names=('mass','lifetime_low','lifetime_high'),meta={'name': 'Lifetime'})

t_lifetime['mass'].unit = 'solMass'
t_lifetime['lifetime_low'].unit = 'Gyr'
t_lifetime['lifetime_high'].unit = 'Gyr'

'''
- mass_yields_001: ejected yield (all heavy elements) in Msolar
            1st col: initial mass of star Msolar
            2nd col: ejected yield Msolar from supernovae
            3rd col: ejected yield Msolar from winds

            M >= 9 Msolar: mp_Z from Maeder 1992 (A & A 264 105)
            Z = 0.001
            M < 9 Msolar: mp_Z from van den Hoek &
            Groenewegen 1997 (A & AS 123 305) Z = 0.001
'''
mass_yields_001 =((0.9, 0, 0.),
             (1., 0, 0.),
             (1.3, 0, 3.09e-3),
             (1.5, 0, 3.53e-3),
             (1.7, 0, 3.77e-3),
             (2., 0, 5.29e-3),
             (2.5, 0, 5.59e-3),
             (3., 0, 5.35e-3),
             (4., 0, 6.56e-3),
             (5., 0, 7.72e-3),
             (7., 0, 9.61e-3),
             (8., 0, 8.78e-3),
             (9., 0.27, 0.),
             (12., 0.83, 0.),
             (15., 1.53, 0.),
             (20., 2.93, 0.),
             (25., 4.45, 0.),
             (40., 4.71, 0.),
             (60., 17.1, 0.),
             (85., 26.7, 0.),
             (120., 41.6, 0.))

t_yields_001 = Table(rows=mass_yields_001, names=('mass','yields_sn','yields_winds'),meta={'name': 'Yields Z=0.001'})

t_yields_001['mass'].unit = 'solMass'
t_yields_001['yields_sn'].unit = 'solMass'
t_yields_001['yields_winds'].unit = 'solMass'


'''
- mass_yields_004: ejected yield (all heavy elements) in Msolar
            1st col: initial mass of star Msolar
            2nd col: ejected yield Msolar from supernovae
            3rd col: ejected yield Msolar from winds

            M >= 9 Msolar: mp_Z from Maeder 1992 (A & A 264 105)
            Z = 0.001 (no other information available)
            M < 9 Msolar: mp_Z from van den Hoek &
            Groenewegen 1997 (A & AS 123 305) Z = 0.004
'''

mass_yields_004 =[(0.9, 0, 1.08e-5),
                 (1., 0, 8.54e-4),
                 (1.3, 0, 1.99e-3),
                 (1.5, 0, 2.57e-3),
                 (1.7, 0, 3.99e-3),
                 (2., 0, 5.86e-3),
                 (2.5, 0, 6.82e-3),
                 (3., 0, 8.43e-3),
                 (4., 0, 5.87e-3),
                 (5., 0, 7.07e-3),
                 (7., 0, 8.88e-3),
                 (8., 0, 9.57e-3),
                 (9., 0.27, 0.),
                 (12., 0.83, 0.),
                 (15., 1.53, 0.),
                 (20., 2.93, 0.),
                 (25., 4.45, 0.),
                 (40., 4.71, 0.),
                 (60., 17.1, 0.),
                 (85., 26.7, 0.),
                 (120., 41.6, 0.)]

t_yields_004 = Table(rows=mass_yields_004, names=('mass','yields_sn','yields_winds'),meta={'name': 'Yields Z=0.004'})

t_yields_004['mass'].unit = 'solMass'
t_yields_004['yields_sn'].unit = 'solMass'
t_yields_004['yields_winds'].unit = 'solMass'

'''
- mass_yields_008: ejected yield (all heavy elements) in Msolar
            1st col: initial mass of star Msolar
            2nd col: ejected yield Msolar from supernovae
            3rd col: ejected yield Msolar from winds

            M >= 9 Msolar: mp_Z from Maeder 1992 (A & A 264 105)
            Z = 0.001 (no other information available)
            M < 9 Msolar: mp_Z from van den Hoek &
            Groenewegen 1997 (A & AS 123 305) Z = 0.008
'''

mass_yields_008 =((0.9, 0, 6.83e-5),
                 (1.0, 0, 1.12e-4),
                 (1.3, 0, 1.70e-3),
                 (1.5, 0, 3.06e-3),
                 (1.7, 0, 3.51e-3),
                 (2., 0, 5.43e-3),
                 (2.5, 0, 6.58e-3),
                 (3., 0, 7.98e-3),
                 (4., 0, 5.44e-3),
                 (5., 0, 6.59e-3),
                 (7., 0, 8.35e-3),
                 (8., 0, 1.02e-2),
                 (9., 0.27, 0.),
                 (12., 0.83, 0.),
                 (15., 1.53, 0.),
                 (20., 2.93, 0.),
                 (25., 4.45, 0.),
                 (40., 4.71, 0.),
                 (60., 17.1, 0.),
                 (85., 26.7, 0.),
                 (120., 41.6, 0.))

t_yields_008 = Table(rows=mass_yields_008, names=('mass','yields_sn','yields_winds'),meta={'name': 'Yields Z=0.008'})

t_yields_008['mass'].unit = 'solMass'
t_yields_008['yields_sn'].unit = 'solMass'
t_yields_008['yields_winds'].unit = 'solMass'

'''
- mass_yields_02: ejected yield (all heavy elements) in Msolar
            1st col: initial mass of star Msolar
            2nd col: ejected yield Msolar from supernovae
            3rd col: ejected yield Msolar from winds

            M >= 9 Msolar: mp_Z from Maeder 1992 (A & A 264 105)
            Z = 0.02 (no other information available)
            M < 9 Msolar: mp_Z from van den Hoek &
            Groenewegen 1997 (A & AS 123 305) Z = 0.019
'''
mass_yields_02 =((0.9, 0, 0.),
             (1., 0, 1.61e-3),
             (1.3, 0, 3.03e-3),
             (1.5, 0, 2.08e-3),
             (1.7, 0, 2.36999e-3),
             (2., 0, 3.94e-3),
             (2.5, 0, 5.5999e-3),
             (3., 0, 6.92e-3),
             (4., 0, 6.24e-3),
             (5., 0, 6.27999e-3),
             (7., 0, 7.73999e-3),
             (8., 0, 8.41e-3),
             (9., 0.173, 0.),
             (12., 0.686, 0.),
             (15., 1.32, 0.),
             (20., 2.73, 0.),
             (25., 4.48, 0.),
             (40., 1.61, 6.4),
             (60., 1.16, 8.69),
             (85., 1.56, 17.75),
             (120., 0.72, 9.39))

t_yields_02 = Table(rows=mass_yields_02, names=('mass','yields_sn','yields_winds'),meta={'name': 'Yields Z=0.019'})

t_yields_02['mass'].unit = 'solMass'
t_yields_02['yields_sn'].unit = 'solMass'
t_yields_02['yields_winds'].unit = 'solMass'

def find_nearest(lookup,value):
    '''
    Take a 2D array and return pair of nearest neighbour values based on first column
    '''
    col1 = lookup[:,0]
    idx = (np.abs(col1-value)).argmin()
    return lookup[idx]
