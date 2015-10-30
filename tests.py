import pytest
from functions import remnant_mass, destruction_timescale, grow_timescale

class TestFunctions:

    def test_remnant_mass_mid(self):
        mass = remnant_mass(10.)
        assert mass.value == 1.5

    def test_remnant_mass_low(self):
        mass = remnant_mass(1.)
        assert mass.value ==  0.552
        
    def test_remnant_mass_high(self):
        mass = remnant_mass(40.)
        assert mass.value ==  10.649999999999999

    def test_destruction(self):
        destroy = destruction_timescale(1000.,6e9,2.1e7)
        destroy = destroy/1e6 #in Myr
        assert destroy.value == 285.7142857142857
    
    def test_graingrowth(self):
        grow = grow_timescale(500.,3.35e9,1.169e9,6.64e-2,(0.671*6.64e-2))
        grow = grow/1e6 #in Myr
        assert grow.value == 86.31618004965111
        

