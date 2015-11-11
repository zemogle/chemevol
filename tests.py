import pytest
import numpy as np
from astropy.table import Table
from functions import remnant_mass, destruction_timescale, \
                    grow_timescale, dust_masses, \
                    inflows, outflows, ejected_gas_mass, astration
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

    def test_destruction(self):
        destroy = destruction_timescale(1000.,6e9,2.1e7, X)
        destroy = destroy/1e6 #in Myr
        assert  285.5 < destroy.value < 285.9

    def test_astration(Self):
        gasmass = 1e10
        sfr = 10
        ast = astration(gasmass, sfr)
        assert ast == 1e-9


    def test_graingrowth(self):
        grow = grow_timescale(500.,3.35e9,1.169e9,6.64e-2,(0.671*6.64e-2))
        grow = grow/1e6 #in Myr
        assert 86.3 < grow.value < 86.4

    def test_outflow_func(Self):
        gas_outflow = outflows(1.0,1.5)
        assert gas_outflow.value == 1.5

    def test_inflow_func(Self):
       gas_inflow = inflows(1.0,1.5)
       assert gas_inflow.value == 1.5

class TestTables:
    '''
    Tests the table entries in lookups haven't changed!
    '''

    def test_lifetimes(self):
        assert lifetime[1][1] == 9.5
        assert lifetime[15][2] == 0.0026

    def test_yields_001(self):
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
