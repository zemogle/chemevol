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

def remnant_mass(m):
    '''
    Calculates the remnant mass of a star of mass m. 
    
    The formulism is based on Ferreras & Silk 2000 (ApJ 532 193)
    which in turn is based on Iben & Tsutukov 1984, Woosley & Weaver 1995.
    '''
    if m <= 9.0:
        rem_mass = 0.106*m + 0.446
    elif (m > 9.0) & (m < 25.0):
        rem_mass = 1.5
    else:
        rem_mass = 0.61*m - 13.75
    return rem_mass
    
def grow_timescale(e,G,SFR,Z,D):
    '''
    Calculates the grain growth timescale in years.
      
    - e: free parameter with MW value set to ~500.  
         e = 2e4 would result in timescales of 1Myr and below.
    - G: gas mass in Msolar
    - D: dust mass in Msolar
    - Z: metallicity mass fraction
    - SFR: star formation rate in units of Msolar/Gyr
    
    Based on Mattsson & Andersen 2012 (MNRAS 423, 38)
    In dust evolution, dMd/dt is proportional to Md/t_grow  
    '''
    
    SFR_in_years = SFR/1e9 # to convert from per Gyr to per yr
    t_grow = G/(e*Z*SFR_in_years)
    t_grow = t_grow/(1-((D/G)/Z)) #to account for metals already locked up in grains
    return t_grow

def destruction_timescale(m,G,SN_rate):
    '''
    Calculates the dust destruction timescale in years.
    
    
    - SN_rate: supernova rate in N/Gyr 
    - G: gas mass in Msolar
    - m: the amount of gas cleared by each supernova event.  
    This is often 100 or 1000 Msolar, appropriate for SNe expanding 
    into galactic densities of 1cm^-3 or 0.1cm^-3 respectively.    
    
    Based on Dwek, Galliano & Jones 2004 (ApJ, 662, 927)
    In dust evolution, dMd/dt is proportional to Md/t_destroy
    '''
    
    SN_rate_in_years = SN_rate/1e9 # to convert from per Gyr to per yr
    t_destroy = G/(m*SN_rate_in_years)
    return t_destroy
    
def initial_mass_function(choice,m):
    '''
    Returns the IMF for a given choice of function and mass range.
    
    - "Chab" selects the Chabrier 2003 IMF (PASP 115 763)  
    - "TopChab" selects a top heavy Chabrier IMF with high mass slope -0.8
    - "Kroup" selects the Kroupa & Weidner 2003 IMF (ApJ 598 1076)
    - "Salp" selects the Salpeter 1955 IMF (ApJ 121 161)
    '''
    
    if choice == "Chab":
        if m <= 0.5:
            imf = np.exp(-1.*(np.log10(m)-np.log10(0.079))**2.)
            imf = (0.85*imf)/((2.*0.69**2.))/m
        else:
            imf = 0.24*(m**-1.3)/m
            
    if choice == "TopChab":
        # If you want to do -0.5 slope, need to change norm factor by
        # 4.72424 
        if m <= 1.0:
            imf = np.exp(-1.*(np.log10(m)-np.log10(0.079))**2.)
            imf = imf*(0.85/2.21896)/((2.*0.69**2.))/m
        else:
            imf = (0.24/2.21896)*(m**-0.8)/m   
            
    if choice == "Kroup":
        if m <= 0.5:
            imf = 0.58*(m**-0.30)/m
        elif (m > 0.5) & (m <= 1.0):
            imf = 0.31*(m**-1.20)/m
        else:
            imf = 0.31*(m**-1.70)/m
            
    if choice == "Salp":
        imf = (0.17/0.990465)*(m**-1.35)/m 
    return imf
