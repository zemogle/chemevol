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
from astropy.table import Table
'''
- lifetime: the lifetime of stars of a given mass (1st column)
            2nd col: lifetime in Gyrs when metallicity < 0.008
            3rd col: lifetime in Gyrs when metallicity >= 0.008
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

t['mass'].unit = 'solMass'
t['lifetime_low'].unit = 'Gyr'
t['lifetime_high'].unit = 'Gyr'

mass_yields
