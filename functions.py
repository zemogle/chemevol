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

from astropy import units as u
import numpy as np
import logging
from lookups import find_nearest, dust_mass_sn, t_yields, lookup_fn, mass_yields
from lookups import yield_names as yn

logger = logging.Logger('chem')

def extra_sfh(sfh, gamma):
    '''
    This extrapolates the SFH provided to start at 0.001Gyr
    with 100 extra steps between 0.001Gyr and the first non-zero
    entry in the SFH list.

    Returns a new SFH list made from joining the
    extrapolated SFH in this routine to the original input SFH file
    '''
    #to start integral at t_0 regardless of when SFH file starts
    t_0 = 1e-3 # we want it to start at 1e-3
    tend_sfh = sfh[1][0] # 1st time array after 0
    # work out difference between t_0 and [1] entry in SFH
    dlogt = (np.log10(tend_sfh) - np.log10(t_0))/100
    norm = sfh[1][1]*(1./np.exp(-1.*gamma*tend_sfh))
    sfr_extra = norm * np.exp(-1.*gamma*t_0)
    sfr_new = sfr_extra
    t_new = t_0
    newlist = []
    #create new array between 0.001 Gyr and start of SFH data
    while t_new < tend_sfh:
        t_new = 10.**(np.log10(t_new)+dlogt)
        sfr_new = norm * np.exp(-1.*gamma*t_new)
        newlist.append([t_new,sfr_new])
    # start from [2:] to account for t[0],t[1] repeated entries
    # when new and in old SFHs combined
    final_sfh = newlist + (sfh.tolist()[2:])
    return final_sfh

def astration(gasmass,sfr):
        mdust_astration = (1./gasmass)*sfr
        return mdust_astration

def remnant_mass(m):
    '''
    Calculates the remnant mass of a star of mass m.

    The formulism is based on Ferreras & Silk 2000 (ApJ 532 193)
    which in turn is based on Iben & Tsutukov 1984, Woosley & Weaver 1995.
    Stars with mass above m_bh don't eject material into the ISM
    '''
    if m <= 9.0:
        rem_mass = 0.106*m + 0.446
    elif (m > 9.0) & (m < 25.0):
        rem_mass = 1.5
    else:
        rem_mass = 0.61*m - 13.75

    rem_mass = rem_mass*u.solMass
    return rem_mass

def initial_mass_function(m, choice):
    '''
    Returns the IMF for a given choice of function and mass range.

    - "Chab", "chab" or "c" selects the Chabrier 2003 IMF (PASP 115 763)
    - "TopChab", "topchab", or "tc" selects a top heavy Chabrier IMF with high mass slope -0.8
    - "Kroup", "kroup" or "k" selects the Kroupa & Weidner 2003 IMF (ApJ 598 1076)
    - "Salp", "salp" or "s" selects the Salpeter 1955 IMF (ApJ 121 161)
    '''

    if (choice == "Chab" or choice == "chab" or choice == "c"):
        if m <= 1.0:
            imf = np.exp(-1.*(np.log10(m)-np.log10(0.079))**2.)
            imf = (0.85*imf)/((2.*0.69**2.))/m
        else:
            imf = 0.24*(m**-1.3)/m

    if (choice == "TopChab" or choice == 'topchab' or choice == "tc"):
        # If you want to do -0.5 slope, need to change norm factor by
        # 4.72424
        if m <= 1.0:
            imf = np.exp(-1.*(np.log10(m)-np.log10(0.079))**2.)
            imf = imf*(0.85/2.21896)/((2.*0.69**2.))/m
        else:
            imf = (0.24/2.21896)*(m**-0.8)/m

    if (choice == "Kroup" or choice == "kroup" or choice == "k"):
        if m <= 0.5:
            imf = 0.58*(m**-0.30)/m
        elif (m > 0.5) & (m <= 1.0):
            imf = 0.31*(m**-1.20)/m
        else:
            imf = 0.31*(m**-1.70)/m

    if (choice == "Salp" or choice == "salp" or choice == "s"):
        imf = (0.17/0.990465)*(m**-1.35)/m
    return imf

def initial_mass_function_integral(choice):
    '''
    Calculates the sum of IMF integral from 0.8 to 120Msun.
    '''
    # initialize
    # only works if dm < 0.005!
    dm = 0.002
    mmax = 120.
    imf_norm = 0.
    m = 0.1
    while m < mmax:
        imf_norm += m*initial_mass_function(m, choice)*dm
        m += dm
    return imf_norm

def ejected_gas_mass(m, sfr, choice):
    '''
    Calculate the ejected mass from stars by mass loss/stellar death
    at time t, needs to be integrated from mass corresponding to
    age of system (tau(m)) -- 120 Msolar

    de/dm = (m-m_R(m)) x SFR(t-tau(m)) x phi(m)
    '''
    if m >= 120.0:
        dej = 0.0
    else:
        dej = (m - (remnant_mass(m).value)) * sfr * initial_mass_function(m, choice)
    return dej

def ejected_metal_mass(m, sfr, zdiff, choice):
    '''
    Calculate the ejected metal mass from stars by mass loss/stellar death
    at time t, needs to be integrated from mass corresponding to
    age of system (tau(m)) -- 120 Msolar

    metals for LIMS are from van Hoek
    Massive stars are from Maeder 1992

    de/dm = (m-m_R(m)*Z(t-taum) + mp_z) x SFR(t-taum x phi(m)
    '''
    if m >= 120.0:
        dej = 0.0
    else:
        massyields = find_nearest(mass_yields, m)
        #sum_mass = metals from winds and SN unless m>40
        # where sum_mass = winds only
        if m <= 40.0:
            sum_mass = massyields[yn.index('yields_sn_001')]+massyields[yn.index('yields_winds_001')]
        else:
            sum_mass = massyields[yn.index('yields_winds_001')]
        dej = ((m - (remnant_mass(m).value))*zdiff + sum_mass) * \
                sfr * initial_mass_function(m, choice)
    return dej

def ejected_dust_mass(m, sfr, zdiff, choice):
    '''
    Calculate the ejected dust mass from stars by mass loss/stellar death
    at time t, needs to be integrated from mass corresponding to
    age of system (tau(m)) -- 120 Msolar

    Re-released dust for LIMS d_LIMS are 0.45 x yields from van den Hoek & Groenewegen
    New dust from fresh heavy elements returned in dust_masses function where
    dust from massive stars (in SN only) are from Todini & Ferrara 2001 (TF01)
    DELTA = fraction of new metals in LIMS (0.45) and SN (from TF01)

    de/dm = (m-m_R(m)*Z(t-taum)*d_LIMS + mp_z*DELTA) x SFR(t-taum x phi(m)
    '''

    delta_LIMS = 0.45
    # no dust from stars with m>40Msun.
    if m >= 40.:
        dej = 0.0
    else:
        massyields = find_nearest(mass_yields, m)
        sum_yields = massyields[yn.index('yields_sn_001')]+massyields[yn.index('yields_winds_001')]
        # sum mass gets dust mass ejected from new elements from SN and winds
        sum_mass = dust_masses(m, sum_yields).value
        dej = ((m - (remnant_mass(m).value))*zdiff*delta_LIMS + sum_mass) * \
                sfr * initial_mass_function(m, choice)
    return dej

def dust_masses(m,yields):
    '''
    This function returns the dust mass ejected by a star
    of initial mass m

    For dust re-released ejecta from stars we multiply the yields
    by a dust condensation efficiency which ranges from 0.16-0.45
    in Morgan & Edmunds 2003 (MNRAS, 343, 427)

    For dust formed from newly processed metals we split into
    two categories: winds from LIMS and SN.

    LIMS: we multiply the metal yields by a dust condensation
    efficiency parameter assumed to be 0.45

    see Figure 3a in Rowlands et al 2014 (MNRAS 441, 1040)

    For high mass stars we use the SN yields of
    Todini & Ferrara 2001 (MNRAS 325 276)
    see Figure 3b in Rowlands et al 2014 (MNRAS 441, 1040)

    m - mass of star
    delta - dust condensation efficiency for LIMS =0.45
    yields - metal yields by mass
    '''
    delta_new_LIMS = 0.45
    if (m <= 8.0):
        dustmass = delta_new_LIMS*yields
    elif (m >= 9.0) & (m <= 40.0):
        #find dust mass from TF01 in dust_mass_sn table
        dustmass = find_nearest(np.array(dust_mass_sn),m)[1]
    else:
        dustmass = 0.
    dustmass = dustmass*u.solMass
    return dustmass

def grow_timescale(e,g,sfr,z,d):
    '''
    Calculates the grain growth timescale in years.

    - e: free parameter with MW value set to ~500.
         e = 2e4 would result in timescales of 1Myr and below.
    - G: gas mass in Msolar
    - D: dust mass in Msolar
    - Z: metallicity mass fraction
    - SFR: star formation rate in units of Msolar/Gyr
    - f_c: fraction of gas in cold dense state for grain growth

    Based on Mattsson & Andersen 2012 (MNRAS 423, 38)
    In dust evolution, dMd/dt is proportional to Md/t_grow
    '''
    sfr_in_years = sfr/1e9 # to convert from per Gyr to per yr
    if z <= 0.:
        t_grow = 0
    else:
        t_grow = g/(e*z*sfr_in_years)
        t_grow = t_grow/(1-((d/g)/z)) #to account for metals already locked up in grains
    t_grow = t_grow*u.year
    return t_grow

def graingrowth(e,g,sfr,z,md,f_c):
        gg = 1e-9*grow_timescale(e,g,sfr,z,md).value #in Gyrs
        if gg <= 0:
            mdust_gg = 0.
        else:
            mdust_gg = md*f_c/gg
        return mdust_gg

def destruction_timescale(destruct,g,supernova_rate):
    '''
    Calculates the dust destruction timescale in years.

    - sn_rate: supernova rate in N/Gyr
    - g: gas mass in Msolar
    - destruct: the amount of gas mass cleared by each supernova event

    destruct is often 100 or 1000 Msolar, appropriate for SNe expanding
    into galactic densities of 1cm^-3 or 0.1cm^-3 respectively.

    Based on Dwek, Galliano & Jones 2004 (ApJ, 662, 927)
    In dust evolution, dMd/dt is proportional to Md/t_destroy
    '''
    if supernova_rate <= 0:
        t_destroy = 0.
    else:
        # to convert from per Gyr to per yr
        supernova_rate = supernova_rate/1e9
        t_destroy = g/(destruct*supernova_rate)
    t_destroy = t_destroy*u.year
    return t_destroy

def inflows(sfr,parameter):
    '''
    Define inflow rate, parameterised by N x SFR
    See Rowlands et al 2014 (MNRAS 441 1040)

    -sfr: SFR at time t
    -parameter: inflow parameter defined in dictionary
    '''
    inflow_rate = sfr*parameter
    inflow_rate = inflow_rate*u.solMass/u.Gyr
    return inflow_rate

def outflows(sfr,parameter):
    '''
    Define outflow rate, parameterised by N x SFR
    See Rowlands et al 2014 (MNRAS 441 1040)

    -sfr: SFR at time t
    -parameter: outflow parameter defined in dictionary
    '''
    outflow_rate = sfr*parameter
    outflow_rate = outflow_rate*u.solMass/u.Gyr
    return outflow_rate

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
                    if not ((dummy == 'SN') or (dummy == 'LIMS') or \
                           (dummy == 'LIMS+SN') or (dummy == 'GG') or\
                           (dummy == 'ALL')):
                        raise ValueError("Oops double check %r:%r" % (run,k))

            # check IMF_fn source options correct
                if (k == 'IMF_fn'):
                    dummy = data_dict[run][k]
                    if not ((dummy == 'Chab') or (dummy == 'chab') or (dummy == 'c') or  \
                            (dummy == 'TopChab') or (dummy == 'topchab') or (dummy == 'tc') or \
                            (dummy == 'Kroup') or (dummy == 'kroup') or (dummy == 'k') or \
                            (dummy == 'Salp') or (dummy == 'salp') or (dummy == 's')):
                        raise ValueError("Oops check %r in %r" % (k,run))
