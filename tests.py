import pytest
import numpy as np
from astropy.table import Table
from functions import remnant_mass, destruction_timescale, destroy_dust, graingrowth, \
                    grow_timescale, dust_masses_fresh, initial_mass_function_integral, \
                    inflows, outflows, ejected_gas_mass, astration, fresh_metals
from lookups import lifetime, mass_yields, dust_mass_sn, find_nearest, \
                    lookup_taum, lookup_fn

class TestFunctions:
    '''
    Test all the function values
    '''

    def test_remnant_mass_mid(self):
        mass = remnant_mass(10.)
        assert mass.value == 1.5

    def test_remnant_mass_low(self):
        mass = remnant_mass(1.)
        assert 0.5515 < mass.value < 0.5525

    def test_remnant_mass_high(self):
        mass = remnant_mass(40.)
        assert 10.6 < mass.value <10.65

    def test_ejected_gasmass_outofbounds(self):
        mass = ejected_gas_mass(120.5,10.5,1)
        assert mass == 0.

    def test_timescale_destruction(self):
        destroy = destruction_timescale(1000.,1.0992e10,6.6e6)
        destroy = destroy*1e-6 #in Myr
        assert  1660 < destroy.value < 1670

    def test_dust_destruction(self):
        dust_sink = destroy_dust(1000.,1.02e10,6.66e6,6.765e08,0.5)[0]
        assert 2.02e8 < dust_sink < 2.29e8

    def test_astration(Self):
        gasmass = 1e10
        sfr = 10
        ast = astration(gasmass, sfr)
        assert ast == 1e-9

    def test_timescale_graingrowth(self):
        grow = grow_timescale(500.,3.35e9,1.169e9,6.64e-2,(0.671*6.64e-2))
        grow = grow*1e-6 #in Myr
        assert 86.3 < grow.value < 86.4

    def test_dust_graingrowth(self):
        dust_ism = graingrowth(500,1.02e10,1e9,0.07,6.765e8,0.5)[0]
        assert  3.200e6 < dust_ism < 3.202e6

    def test_outflow_func(Self):
        gas_outflow = outflows(1.0,1.5)
        assert gas_outflow.value == 1.5

    def test_inflow_func(Self):
       gas_inflow = inflows(1.0,1.5)
       assert gas_inflow.value == 1.5

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
        dust_mass = dust_masses_fresh(1.0,0.001).value
        assert dust_mass == 0

    def test_fresh_dust_mass_midmass_lowmetals(self):
        dust_mass = dust_masses_fresh(5.0,0.001).value
        assert 0.00346 < dust_mass < 0.00348

    def test_fresh_dust_mass_highmass_lowmetals(self):
        dust_mass = dust_masses_fresh(30.0,0.001).value
        assert dust_mass == 1.0

    def test_fresh_dust_mass_highermass_lowmetals(self):
        dust_mass = dust_masses_fresh(40.0,0.001).value
        assert dust_mass == 0.4

    def test_fresh_dust_mass_lowmass_highmetals(self):
        dust_mass = dust_masses_fresh(1.0,0.02).value
        assert  0.000070 < dust_mass < 0.000073

    def test_fresh_dust_mass_midmass_highmetals(self):
        dust_mass = dust_masses_fresh(2.0,0.02).value
        assert 2.442e-3 < dust_mass < 2.446e-3

    def test_fresh_dust_mass_highmass_highmetals(self):
        dust_mass = dust_masses_fresh(30.0,0.02).value
        assert dust_mass == 1.0

    def test_fresh_dust_mass_highermass_highmetals(self):
        dust_mass = dust_masses_fresh(40.0,0.02).value
        assert dust_mass == 0.4

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
    def test_stellarlifetime_lookup(self):
        lifetime_cols = {'low_metals':1, 'high_metals':2}
        taum =  lookup_taum(41.,lifetime_cols['low_metals'])
        assert taum == 0.0048999999999999998
