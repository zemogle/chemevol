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
to make sure input dictionary is set up correctly
'''

def validate_initial_dict(keysdict, data_dict):
    '''
    Validate the initial data
    '''
    for run,keys in data_dict.items():
        for k in keysdict:
            try:
                dummy = data_dict[run][k]
            except KeyError:
                print("Oops key %r is missing in %r" % (k,run))
            else:

            # check gasmass, gamma, inflows and outflows are numbers:
                if (k == 'gasmass_init') or (k == 'inflows') \
                    or (k == 'outflows') or (k == 'gamma'):
                    try:
                        dummy = int(data_dict[run][k])
                    except ValueError:
                        print("Oops we were expecting a number in %r:%r" % (run,k))

            # check SFH is a string :
                if (k == 'SFH'):
                    if not isinstance(data_dict[run][k], basestring):
                        raise TypeError("Oops %r:%r should be a string" % (run,k))

            # check dust_source options correct
                if (k == 'dust_source'):
                    dummy = data_dict[run][k]
                    if not ((dummy == 'SN') or (dummy == 'Sn') or (dummy == 'sn') or \
                            or (dummy == 'LIMS') or (dummy == 'lims') or (dummy == 'Lims') or \
                           (dummy == 'LIMS+SN') or (dummy == 'lims+sn') or
                           (dummy == 'Lims+sn') or (dummy =="sn+lims") or (dummy =="SN+LIMS") or \
                           (dummy == 'GG') or (dummy == 'gg') or (dummy == 'Gg') or \
                           (dummy == 'ALL') or (dummy == 'all') or (dummy == 'All') ):
                        raise ValueError("Oops double check %r:%r" % (run,k))

            # check IMF_fn source options correct
                if (k == 'IMF_fn'):
                    dummy = data_dict[run][k]
                    if not ((dummy == 'Chab') or (dummy == 'chab') or (dummy == 'c') or  \
                            (dummy == 'TopChab') or (dummy == 'topchab') or (dummy == 'tc') or \
                            (dummy == 'Kroup') or (dummy == 'kroup') or (dummy == 'k') or \
                            (dummy == 'Salp') or (dummy == 'salp') or (dummy == 's')):
                        raise ValueError("Oops check %r in %r" % (k,run))
