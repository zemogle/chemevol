import pytest
from functions import remnant_mass, destruction_timescale, grow_timescale
from lookups import lifetime, mass_yields_001, mass_yields_004, mass_yields_008, mass_yields_02
from main import initial_galaxy_params

class TestFunctions:

    def test_remnant_mass_mid(self):
        mass = remnant_mass(10.)
        assert mass.value == 1.5

    def test_remnant_mass_low(self):
        mass = remnant_mass(1.)
        assert 0.5515 < mass.value < 0.5525
        
    def test_remnant_mass_high(self):
        mass = remnant_mass(40.)
        assert 10.6 < mass.value <10.65

    def test_destruction(self):
        destroy = destruction_timescale(1000.,6e9,2.1e7)
        destroy = destroy/1e6 #in Myr
        assert  280 < destroy.value < 290
    
    def test_graingrowth(self):
        grow = grow_timescale(500.,3.35e9,1.169e9,6.64e-2,(0.671*6.64e-2))
        grow = grow/1e6 #in Myr
        assert 80 < grow.value < 90 

class TestTables:
    
    def test_lifetimes(self):
        assert lifetime[1][1] == 9.5
        assert lifetime[15][2] == 0.0026        

    def test_yields_001(self):
        assert mass_yields_001[1][1] == 0.
        assert mass_yields_001[12][1] == 0.27   

    def test_yields_004(self):
        assert mass_yields_004[1][2] == 8.54e-4
        assert mass_yields_004[20][1] == 41.6 
        
    def test_yields_008(self):
        assert mass_yields_008[1][2] == 1.12e-4
        assert mass_yields_008[20][2] == 0. 
    
    def test_yields_002(self):
        assert mass_yields_02[1][2] == 1.61e-3
        assert mass_yields_02[19][2] == 17.75

 


            