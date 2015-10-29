'''
Chemevol - Python package to read in a star formation history file, input galaxy parameters and run chemical evolution to determine gas, metals, dust evolution of galaxies.  The code is based on Morgan & Edmunds 2003 (MNRAS, 343, 427) and described in detail in Rowlands et al 2014 (MNRAS, 441, 1040).   

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

def remnant_mass(mass):
    '''
    Calculates the remnant mass of a star of mass m. 
    
    The formulism is based on Ferreras & Silk 2000 (ApJ 532 193)
    which in turn is based on Iben & Tsutukov 1984, Woosley & Weaver 1995.
    '''
    if mass <= 9.0:
        rem_mass = 0.106*mass+0.446
    elif (mass < 25.0) & (mass > 9.0):
        rem_mass = 1.5
    else:
        rem_mass = 0.61*m-13.75
    return rem_mass
    
def grow_timescale(choice,e,G,SFR,Z,D):
    '''
    Calculates the grain growth timescale in years.
    
    Based on Mattsson & Andersen 2012 (MNRAS 423, 38).  
    The grain growth timescale is depends on the metallicity (Z), 
    gas mass (G), dust mass (D) and SFR. e is a free parameter. 
    Typical values for the Milky Way would be e=500, 
    e = 2e4 would result in timescales of 1Myr and below. 
    
    In dust evolution, dMd/dt is proportional to Md/t_grow  
    '''
    SFR_in_years = SFR/1e9 # to convert from per Gyr to per yr
    t_grow = G/(e*Z*SFR_in_years)
    t_grow = t_grow/(1-((D/G)/Z)) #to account for metals already locked up in grains
    return t_grow

def destruction_timescale(m,G,SN_rate):
    '''
    Calculates the dust destruction timescale in years.
    
    Based on Dwek, Galliano & Jones 2004 (ApJ, 662, 927).  
    This depends on the supernova rate and the gas mass (G).  
    m is a free parameter and is the amount of gas cleared by each supernova event.  
    This is often selected to be 100 or 1000Msun, appropriate for SN expanding 
    into galactic densities of 1cm^-3 or 0.1cm^-3 respectively.    
    
    In dust evolution, dMd/dt is proportional to Md/t_destroy
    '''
    SN_rate_in_years = SN_rate/1e9 # to convert from per Gyr to per yr
    t_destroy = G/(m*SN_rate_in_years)
    return t_destroy
