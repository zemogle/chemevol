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

********************************************************************************
'''
'''
This sets all the unit tests for the code. To run: py.test tests.py
'''
import pytest
import numpy as np
from astropy.table import Table
from functions import remnant_mass, destruction_timescale, destroy_dust, graingrowth, \
                    grow_timescale, dust_masses_fresh, initial_mass_function_integral, \
                    inflows, ejected_gas_mass, astration, fresh_metals, \
                    ejected_dust_mass, imf_chab, outflows_feldmann, gas_inandout,\
                    metals_inandout, dust_inandout

from lookups import lifetime, mass_yields, dust_mass_sn, find_nearest, \
                    lookup_taum, lookup_fn

dustchoice_all = {'sn' : True, 'lims' : True, 'gg':True}
dustchoice_sn = {'sn' : True, 'lims' : False, 'gg':False}
dustchoice_lims = {'sn' : False, 'lims' : True, 'gg':False}
dustchoice_snlims = {'sn' : True, 'lims' : True, 'gg':False}
dustchoice_gg = {'sn' : False, 'lims' : False, 'gg':True}


class TestFunctions:
    '''
    Test all the function values
    '''

    '''
    Test remnant masses
    '''
    def test_remnant_mass_mid(self):
        mass = remnant_mass(10.)
        assert mass == 1.5

    def test_remnant_mass_low(self):
        mass = remnant_mass(1.)
        assert 0.5515 < mass < 0.5525

    def test_remnant_mass_high(self):
        mass = remnant_mass(40.)
        assert 10.6 < mass <10.65

    '''
    Test astration
    '''
    def test_astration(Self):
        gasmass = 1e10
        sfr = 10
        ast = astration(1,gasmass, sfr)
        assert ast == 1e-9

    '''
    Test grain growth and destruction
    '''
    def test_timescale_destruction(self):
        destroy = destruction_timescale(True,1000.,1.0992e10,6.6e6)
        destroy = destroy*1e-6 #in Myr
        assert  1660 < destroy < 1670

    def test_dust_destruction(self):
        dust_sink = destroy_dust(True,1000.,1.02e10,6.66e6,6.765e08,0.5)[0]
        assert 2.02e8 < dust_sink < 2.29e8

    def test_timescale_destruction_off(self):
        destroy = destruction_timescale(False,1000.,1.0992e10,6.6e6)
        destroy = destroy*1e-6 #in Myr
        assert destroy == 0

    def test_dust_destruction_off(self):
        dust_sink = destroy_dust(False,1000.,1.02e10,6.66e6,6.765e08,0.5)[0]
        assert dust_sink == 0

    def test_timescale_graingrowth(self):
        grow = grow_timescale(True,500.,3.35e9,1.169e9,6.64e-2,(0.671*6.64e-2))
        grow = grow*1e-6 #in Myr
        assert 86.3 < grow < 86.4

    def test_timescale_graingrowth_off(self):
        grow = grow_timescale(False,500.,3.35e9,1.169e9,6.64e-2,(0.671*6.64e-2))
        grow = grow*1e-6 #in Myr
        assert grow == 0

    def test_dust_graingrowth(self):
        dust_ism = graingrowth(True,500,1.02e10,1e9,0.07,6.765e8,0.5)[0]
        assert  3.200e6 < dust_ism < 3.202e6

    def test_dust_graingrowth_off(self):
        dust_ism = graingrowth(False,500,1.02e10,1e9,0.07,6.765e8,0.5)[0]
        assert dust_ism == 0

    '''
    Test inflows and outflows
    '''
    def test_inflow_func(Self):
       gas_inflow = inflows(1.0,1.5)
       assert gas_inflow == 1.5

    def test_gas_inandout_f1(Self):
       val_in,val_out = gas_inandout(False,False,1.5,1,1e10)
       assert val_in == 0 and val_out == 0

    def test_gas_inandout_f2(Self):
       val_in,val_out = gas_inandout(True,False,0,1,1e10)
       assert val_in == 0 and val_out == 0

    def test_gas_inandout_true(Self):
       val_in,val_out = gas_inandout(True,True,1.5,1,1e10)
       assert val_in > 0 and val_out > 0

    def test_metals_inandout_f1(Self):
       val_in,val_out = metals_inandout(False,1.5,1e-3,False,False,1,0.02,1e10)
       assert val_in == 0 and val_out == 0

    def test_metals_inandout_f2(Self):
       val_in,val_out = metals_inandout(True,0,1e-3,False,False,1,0.02,1e10)
       assert val_in == 0 and val_out == 0

    def test_metals_inandout_f3(Self):
       val_in,val_out = metals_inandout(True,1,0,False,False,1,0.02,1e10)
       assert val_in == 0 and val_out == 0

    def test_metals_inandout_f4(Self):
       val_in,val_out = metals_inandout(True,1,1e-4,True,False,1,0.02,1e10)
       assert val_in > 0 and val_out == 0

    def test_metals_inandout_true(Self):
       val_in,val_out = metals_inandout(True,1,1e-4,True,True,1,0.02,1e10)
       assert val_in > 0 and val_out > 0

    def test_dust_inandout_f1(Self):
       val_in,val_out = dust_inandout(False,1.5,1e-4,False,False,1,0.01,1e10)
       assert val_in == 0 and val_out == 0

    def test_dust_inandout_f2(Self):
       val_in,val_out = dust_inandout(True,0,1e-3,False,False,1,0.01,1e10)
       assert val_in == 0 and val_out == 0

    def test_dust_inandout_f3(Self):
       val_in,val_out = dust_inandout(True,1,0,False,False,1,0.01,1e10)
       assert val_in == 0 and val_out == 0

    def test_dust_inandout_f4(Self):
       val_in,val_out = dust_inandout(True,1,1e-4,True,False,1,0.01,1e10)
       assert val_in > 0 and val_out == 0

    def test_dust_inandout_true(Self):
       val_in,val_out = dust_inandout(True,1,1e-4,True,True,1,0.02,1e10)
       assert val_in > 0 and val_out > 0

    def test_outflows_feld_intm(self):
        val = outflows_feldmann(1,6.8e10)
        assert 2 < val < 4

    def test_outflows_feld_off(self):
        val = outflows_feldmann(1,1e7)
        assert val  == 0

    def test_outflows_feld_lowm(self):
        val = outflows_feldmann(1,1e8)
        assert val <= 30

    '''
    Test IMF
    '''
    def test_imf_integral_chab(self):
        unity = initial_mass_function_integral('Chab')
        assert 0.99 < unity < 1.09

    def test_imf_integral_kroup(self):
        unity = initial_mass_function_integral('Kroup')
        assert 0.99 < unity < 1.09

    def test_imf_integral_salp(self):
        unity = initial_mass_function_integral('Salp')
        assert 0.99 < unity < 1.09

    def test_imf_integral_topchab(self):
        unity = initial_mass_function_integral('TopChab')
        assert 0.99 < unity < 1.09

    '''
    Test the metals and dust that are ejected by stars
    '''
    def test_fresh_metals_lowmassa_lowmetals(self):
        mass_yields = fresh_metals(1.,0.001)
        assert mass_yields == 0

    def test_fresh_metals_lowmassb_lowmetals(self):
        mass_yields = fresh_metals(2.,0.001)
        assert mass_yields == 0.00529

    def test_fresh_metals_midmass_lowmetals(self):
        mass_yields = fresh_metals(30.,0.001)
        assert mass_yields == 4.45

    def test_fresh_metals_highmass_lowmetals(self):
        mass_yields = fresh_metals(119,0.001)
        assert mass_yields == 0

    def test_fresh_metals_lowmassa_highmetals(self):
        mass_yields = fresh_metals(1.,0.02)
        assert mass_yields == 1.61e-4

    def test_fresh_metals_lowmassb_highmetals(self):
        mass_yields = fresh_metals(2.,0.02)
        assert mass_yields == 0.00543

    def test_fresh_metals_midmass_highmetals(self):
        mass_yields = fresh_metals(30.,0.02)
        assert mass_yields == 4.48

    def test_fresh_metals_highmass_highmetals(self):
        mass_yields = fresh_metals(119,0.02)
        assert mass_yields == 9.39

    def test_fresh_dust_mass_lowmass_lowmetals(self):
        dust_mass = dust_masses_fresh(dustchoice_all,1.0,1.0,0.001)
        assert dust_mass == 0

    def test_fresh_dust_mass_midmass_lowmetals(self):
        dust_mass = dust_masses_fresh(dustchoice_all,1.0,5.0,0.001)
        assert 0.00346 < dust_mass < 0.00348

    def test_fresh_dust_mass_highmass_lowmetals(self):
        dust_mass = dust_masses_fresh(dustchoice_all,1.0,30.0,0.001)
        assert dust_mass == 1.0

    def test_fresh_dust_mass_highermass_lowmetals(self):
        dust_mass = dust_masses_fresh(dustchoice_all,1.0,40.0,0.001)
        assert dust_mass == 0.4

    def test_fresh_dust_mass_lowmass_highmetals(self):
        dust_mass = dust_masses_fresh(dustchoice_all,1.0,1.0,0.02)
        assert  0.000070 < dust_mass < 0.000073

    def test_fresh_dust_mass_midmass_highmetals(self):
        dust_mass = dust_masses_fresh(dustchoice_all,1.0,2.0,0.02)
        assert 2.442e-3 < dust_mass < 2.446e-3

    def test_fresh_dust_mass_highmass_highmetals(self):
        dust_mass = dust_masses_fresh(dustchoice_all,1.0,30.0,0.02)
        assert dust_mass == 1.0

    def test_fresh_dust_mass_highermass_highmetals(self):
        dust_mass = dust_masses_fresh(dustchoice_all,1.0,40.0,0.02)
        assert dust_mass == 0.4

    def test_fresh_metals_highmass_highmetals(self):
        mass_yields = fresh_metals(119,0.02)
        assert mass_yields == 9.39

    def test_fresh_dust_mass_lowmass_lowmetals_no(self):
        dust_mass = dust_masses_fresh(dustchoice_sn,1,2.0,0.001)
        assert dust_mass == 0

    def test_fresh_dust_mass_midmass_lowmetals_no(self):
        dust_mass = dust_masses_fresh(dustchoice_gg,1,8.0,0.001)
        assert dust_mass == 0

    def test_fresh_dust_mass_highmass_lowmetals_no(self):
        dust_mass = dust_masses_fresh(dustchoice_gg,1,30.0,0.001)
        assert dust_mass == 0

    def test_fresh_dust_mass_highermass_lowmetals_no(self):
        dust_mass = dust_masses_fresh(dustchoice_sn,1,40.0,0.001)
        assert dust_mass == 0.4

    def test_fresh_dust_mass_lowmass_highmetals_no(self):
        dust_mass = dust_masses_fresh(dustchoice_lims,1,1.0,0.02)
        assert  0.000070 < dust_mass < 0.000073

    def test_fresh_dust_mass_midmass_highmetals_no(self):
        dust_mass = dust_masses_fresh(dustchoice_sn,1,2.0,0.02)
        assert dust_mass  == 0

    def test_fresh_dust_mass_highmass_highmetals_no(self):
        dust_mass = dust_masses_fresh(dustchoice_snlims,1,30.0,0.02)
        assert dust_mass == 1.0

    def test_freshdust_sn_only(self):
        dustmass_low = dust_masses_fresh(dustchoice_sn, 1.0,4.9, 0.002)
        dustmass_high = dust_masses_fresh(dustchoice_sn, 1.0,15, 0.002)
        assert dustmass_low == 0 and dustmass_high == 0.5

    def test_freshdust_lims_only_highmetals(self):
        dustmass_low = dust_masses_fresh(dustchoice_lims, 1.0,4.9, 0.01)
        dustmass_high = dust_masses_fresh(dustchoice_lims, 1.0,15, 0.01)
        assert dustmass_low == 0.0029655000000000003 and dustmass_high == 0.0

    def test_ejected_gasmass_outofbounds(self):
        mass = ejected_gas_mass(120.5,10.5,1)
        assert mass == 0.

    def test_ejected_dust_mass(self):
        dustmass_all = ejected_dust_mass(dustchoice_all,1,5.0,10389385569.1, 7.70733489684e-06, 0.000166298678684,imf_chab)
        dustmass_sn = ejected_dust_mass(dustchoice_sn,1,5.0,10389385569.1, 7.70733489684e-06, 0.000166298678684,imf_chab)
        dustmass_lims = ejected_dust_mass(dustchoice_lims,1,5.0,10389385569.1, 7.70733489684e-06, 0.000166298678684,imf_chab)
        dustmass_both = ejected_dust_mass(dustchoice_snlims,1,5.0,10389385569.1, 7.70733489684e-06, 0.000166298678684,imf_chab)
        assert dustmass_all == 214655.06895476999 and dustmass_both == 214655.06895476999 and \
                dustmass_sn == 0 and dustmass_lims == 214655.06895476999

class TestInitials:
    '''
    Tests whether things are turned on or off correctly from init file
    '''
    def test_destruction_turned_off(self):
        dustmass = destroy_dust(False,1000,4e10,0.006,1e5,0.5)[0]
        assert dustmass == 0

    def test_graingrowth_turned_off(Self):
        dustmass = graingrowth(False,500,1.02e10,1e9,0.07,6.765e8,0.5)[0]
        assert dustmass == 0

class TestTables:
    '''
    Tests the table entries in lookups haven't changed!
    '''

    def test_lifetimes(self):
        assert lifetime[1][1] == 9.5
        assert lifetime[15][2] == 0.0026

    def test_yields_return_by_mass(self):
        assert mass_yields[1][1] == 0.
        assert mass_yields[12][1] == 0.27
        assert mass_yields[1][4] == 8.54e-4
        assert mass_yields[20][3] == 41.6
        assert mass_yields[1][6] == 1.12e-4
        assert mass_yields[20][6] == 9.39
        assert mass_yields[0][8] == 6.83e-3
        assert mass_yields[19][8] == 17.75

class Testlookups:
    '''
    Test the lookup functions
    '''
    def test_stellarlifetime_lookup_lowZ(self):
        lifetime_cols = {'low_metals':1, 'high_metals':2}
        taum =  lookup_taum(41.,lifetime_cols['low_metals'])
        assert taum == 0.0048999999999999998

    def test_stellarlifetime_lookup_highZ(self):
        lifetime_cols = {'low_metals':1, 'high_metals':2}
        taum =  lookup_taum(60.,lifetime_cols['high_metals'])
        assert taum == 0.0034

    def test_minimum_mass_lookup_lowZ(self):
        t_lifetime = Table(rows=lifetime, names=('mass','lifetime_low_metals','lifetime_high_metals'),meta={'name': 'Lifetime'})
        minimum_mass =  lookup_fn(t_lifetime,'lifetime_low_metals',0.029)['mass']
        assert minimum_mass == 9.0

    def test_minimum_mass_lookup_highZ(self):
        t_lifetime = Table(rows=lifetime, names=('mass','lifetime_low_metals','lifetime_high_metals'),meta={'name': 'Lifetime'})
        minimum_mass =  lookup_fn(t_lifetime,'lifetime_high_metals',1.75)['mass']
        assert minimum_mass == 2.0
