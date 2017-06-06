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
import numpy as np
from numpy import abs, array
import logging
from lookups import find_nearest, dust_mass_sn, t_yields, t_lifetime, \
                    lookup_fn, lookup_taum, mass_yields, oxymass_yields
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
    dlogt = (np.log10(tend_sfh) - np.log10(t_0))/1000
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

def astration(parameter,gasmass,sfr):
        astration_term = (parameter/gasmass)*sfr
        return astration_term

def remnant_mass(m):
    '''
    Calculates the remnant mass of a star of mass m.

    The formulism is based on Ferreras & Silk 2000 (ApJ 532 193)
    which in turn is based on Iben & Tsutukov 1984, Woosley & Weaver 1995.
    This accounts for stars with mass above 40Msun not returning all their material into the ISM
    '''
    if m <= 9.0:
        rem_mass = 0.106*m + 0.446
    elif (m > 9.0) & (m < 25.0):
        rem_mass = 1.5
    else:
        rem_mass = 0.61*m - 13.75

    rem_mass = rem_mass#*u.solMass
    return rem_mass

def imf_chab(m):
    if m <= 1.0:
        imf = np.exp(-1.*(np.log10(m)+1.1023729) * (np.log10(m)+1.1023729))
        imf = (0.85*imf)/0.952199/m
    else:
        imf = 0.24*(m**-1.3)/m
    return imf

def imf_topchab(m):
    # If you want to do -0.5 slope, need to change norm factor by
    # 4.72424
    if m <= 1.0:
        imf = np.exp(-1.*(np.log10(m)-np.log10(0.079))**2.)
        imf = imf*(0.85/2.21896)/((2.*0.69**2.))/m
    else:
        imf = 0.1081587*(m**-0.8)/m
    return imf

def imf_kroup(m):
    if m <= 0.5:
        imf = 0.58*(m**-0.30)/m
    elif (m > 0.5) & (m <= 1.0):
        imf = 0.31*(m**-1.20)/m
    else:
        imf = 0.31*(m**-1.70)/m
    return imf

def imf_salp(m):
    imf = (0.17/0.990465)*(m**-1.35)/m
    return imf

def initial_mass_function(m, choice):
    '''
    Returns the IMF for a given choice of function and mass range.

    - "Chab", "chab" or "c" selects the Chabrier 2003 IMF (PASP 115 763)
    - "TopChab", "topchab", or "tc" selects a top heavy Chabrier IMF with high mass slope -0.8
    - "Kroup", "kroup" or "k" selects the Kroupa & Weidner 2003 IMF (ApJ 598 1076)
    - "Salp", "salp" or "s" selects the Salpeter 1955 IMF (ApJ 121 161)
    '''

    if (choice == "Chab" or choice == "chab" or choice == "c"):
        imf = imf_chab(m)

    elif (choice == "TopChab" or choice == 'topchab' or choice == "tc"):
        imf = imf_topchab(m)

    elif (choice == "Kroup" or choice == "kroup" or choice == "k"):
        imf = imf_kroup(m)

    elif (choice == "Salp" or choice == "salp" or choice == "s"):
        imf = imf_salp(m)
    return imf

def initial_mass_function_integral(choice):
    '''
    Calculates the sum of IMF integral from 0.8 to 120Msun
    This isn't called in the code but there to check if
    normalised IMF works
    '''
    # only works if dm < 0.005!
    dm = 0.002
    mmax = 120.
    imf_norm = 0.
    m = 0.1
    while m < mmax:
        imf_norm += m*initial_mass_function(m, choice)*dm
        m += dm
    return imf_norm

def ejected_gas_mass(m, sfr, imf):
    '''
    Calculate the ejected mass from stars by mass loss/stellar death
    at time t, needs to be integrated from mass corresponding to
    age of system (tau(m)) -- 120 Msolar

    de/dm = (m-m_R(m)) x SFR(t-tau(m)) x phi(m)
    '''
    if m >= 120.0:
        dej = 0.0
    else:
        dej = (m - (remnant_mass(m))) * sfr * imf(m)
    return dej

def fresh_metals(m, metallicity):
    '''
    Function to return the fresh new elements made by stars
    These are metallicity dependent and calls mass_yields table
    in lookups.py

    metals for LIMS are from van Hoek
    Massive stars are from Maeder 1992

    For m > 40, then only winds contribute to ejected metals
    For m <= 40, winds + SNe contribute
    '''
    massyields = find_nearest(mass_yields, m)
    if metallicity <= 0.0025:
        if m <= 40:
            sum_yields = massyields[yn.index('yields_sn_001')]+massyields[yn.index('yields_winds_001')]
        else:
            sum_yields = massyields[yn.index('yields_winds_001')]
    elif metallicity <= 0.006:
        if m <= 40:
            sum_yields = massyields[yn.index('yields_sn_004')]+massyields[yn.index('yields_winds_004')]
        else:
            sum_yields = massyields[yn.index('yields_winds_004')]
    elif metallicity <= 0.01:
        if m <= 40:
            sum_yields = massyields[yn.index('yields_sn_008')]+massyields[yn.index('yields_winds_008')]
        else:
            sum_yields = massyields[yn.index('yields_winds_008')]
    else:
        if m <= 40:
            sum_yields = massyields[yn.index('yields_sn_02')]+massyields[yn.index('yields_winds_02')]
        else:
            sum_yields = massyields[yn.index('yields_winds_02')]
    return sum_yields

def fresh_oxygen(m, metallicity):
    '''
    Function to return the fresh oxygen made by stars
    These are metallicity dependent and calls oxymass_yields table
    in lookups.py

    metals for LIMS are from van Hoek
    Massive stars are from Maeder 1992

    For m > 40, then only winds contribute to ejected metals
    For m <= 40, winds + SNe contribute
    '''
    oxymassyields = find_nearest(oxymass_yields, m)
    if metallicity <= 0.0025:
        if m <= 40:
            sum_yields = oxymassyields[yn.index('yields_sn_001')]+oxymassyields[yn.index('yields_winds_001')]
        else:
            sum_yields = oxymassyields[yn.index('yields_winds_001')]
    elif metallicity <= 0.006:
        if m <= 40:
            sum_yields = oxymassyields[yn.index('yields_sn_004')]+oxymassyields[yn.index('yields_winds_004')]
        else:
            sum_yields = oxymassyields[yn.index('yields_winds_004')]
    elif metallicity <= 0.01:
        if m <= 40:
            sum_yields = oxymaessyields[yn.index('yields_sn_008')]+oxymassyields[yn.index('yields_winds_008')]
        else:
            sum_yields = oxymassyields[yn.index('yields_winds_008')]
    else:
        if m <= 40:
            sum_yields = oxymassyields[yn.index('yields_sn_02')]+oxymassyields[yn.index('yields_winds_02')]
        else:
            sum_yields = oxymassyields[yn.index('yields_winds_02')]
    return sum_yields

def ejected_oxygen_mass(m, sfrdiff, oxydiff, metallicity, imf):
    '''
    Calculate the ejected oxygen mass from stars by mass loss/stellar death
    at time t, needs to be integrated from mass corresponding to
    age of system (tau(m)) -- 120 Msolar.

    It calls function fresh_oxygen to find correct mass of new
    oxygen ejected by stars of mass m (metallicity dependent)

    de (m,t) = (m-m_R(m)*Z(t-taum) + mp(m,Z)) x SFR(t-taum x phi(m)
    '''
    if m >= 120.0:
        dej = 0.0
    else:
        dej = ((m - (remnant_mass(m)))*oxydiff + fresh_oxygen(m, metallicity)) * \
                sfrdiff * imf(m)
    return dej

def ejected_metal_mass(m, sfrdiff, zdiff, metallicity, imf):
    '''
    Calculate the ejected metal mass from stars by mass loss/stellar death
    at time t, needs to be integrated from mass corresponding to
    age of system (tau(m)) -- 120 Msolar.

    It calls function fresh_metals to find correct mass of new
    heavy elements ejected by stars of mass m (metallicity dependent)

    de (m,t) = (m-m_R(m)*Z(t-taum) + mp(m,Z)) x SFR(t-taum x phi(m)
    '''
    if m >= 120.0:
        dej = 0.0
    else:
        dej = ((m - (remnant_mass(m)))*zdiff + fresh_metals(m, metallicity)) * \
                sfrdiff * imf(m)
    return dej

def ejected_dust_mass(choice, delta_lims, reduce_sn, m, sfrdiff, zdiff, metallicity, imf):
    '''
    Calculate the ejected dust mass from stars by mass loss/stellar death
    at time t, needs to be integrated from mass corresponding to
    age of system (tau(m)) -- 120 Msolar

    1st term: dust re-released by stars
              Calculated by ejected gas mass * Z(t-taum) * dust condensation efficiency (delta_lims_rec)
              delta_lims_rec ranges from 0.15-0.4 for stars in Morgan & Edmunds 2003 (MNRAS, 343, 427)
              delta_lims_rec is set to 0.15 in De Vis et al 2017b in press 1705.02340

    2nd term: new dust from fresh heavy elements returned in dust_masses function where
              dust from massive stars (in SN only) are from Todini & Ferrara 2001 (TF01) and dust
              from Van den Hoek & Groenewegen:
              delta_lims = fraction of new metals in LIMS that condense into dust, set by user in inits
              NB delta_lims ranges from 0.15-0.4 for stars in Morgan & Edmunds 2003 (MNRAS, 343, 427)
              NB delta_lims is set to 0.15 in De Vis et al 2017b in press 1705.02340
              md_SN = dust mass SN (from TF01)

    de/dm = (m-m_R(m)*Z(t-taum)*delta_lims_rec + (mp*delta_lims_new)+md_SN) x SFR(t-taum x phi(m)

    In:
    -- choice: array of dust source choices, set by user in inits
               0th element = sn value 1 or 0
               1st element = lims value 1 or 0
               2nd element = grain growth value 1 or 0
    -- delta_lims: efficiency of fresh metals condensing into dust in lims
    -- reduce_sn: factor to reduce SN dust contribution, set by user in inits
    -- m: mass of star
    -- metallicity: metal mass fraction Mz/Mg
    -- sfrdiff: SFR calculated at t-taum (when stars that are dying now were born)
    -- zdiff: metallicity calculated at t-taum (when stars that are dying now were born)
    -- imf choice, set by user in inits

    '''
    # If LIMS is turned on or off
    choice_lims = choice['lims']
    # condensation efficiency of recycled stars in LIMS ONLY for m <= 8Msun
    if m <= 8. and choice_lims:
        delta_LIMS_rec = 0.15
    else:
        delta_LIMS_rec = 0.

    if m > 40.: # no dust from stars with m>40Msun.
        dej = 0.0
        dej_fresh = 0
        dej_recycled = 0
    else: # recycled + fresh dust
        # read in dust mass from freshly formed metals as function m and Z (chi_2 * LIMS yields)
        dej_fresh = dust_masses_fresh(choice, delta_lims, reduce_sn, m, metallicity) * sfrdiff * imf(m)
        # recycled dust = LIMS condensation efficiency
        dej_recycled = ((m - (remnant_mass(m)))*zdiff*delta_LIMS_rec) * sfrdiff * imf(m)

    dej = dej_recycled + dej_fresh
    return dej

def dust_masses_fresh(choice, delta_lims, reduce_sn, m, metallicity):
    '''
    This function returns the dust mass ejected by a star
    of initial mass m made from freshly synthesised elements

    For dust formed from newly processed metals we split into
    two categories: winds from LIMS and SN.

    LIMS: we multiply the fresh metal yields by a dust condensation
    efficiency parameter input by user (delta_lims)

    For high mass stars we use the SN yields of
    Todini & Ferrara 2001 (MNRAS 325 276) in lookups.py

    In:
    -- choice: array of dust source choices, set by user in inits
               0th element = sn value 1 or 0
               1st element = lims value 1 or 0
               2nd element = grain growth value 1 or 0
    -- reduce_sn: factor to reduce SN dust contribution, set by user in inits
    -- m: mass of star
    -- metallicity: metal mass fraction Mz/Mg

    yields - metal yields by mass

    See Figure 3 in Rowlands et al 2014 (MNRAS 441, 1040)
    '''

#    delta_new_LIMS = 0.15
    if (m <= 8.0) and choice['lims']:
        dustmass = delta_lims * fresh_metals(m, metallicity)
    elif (m > 8.0) and (m <= 40.0) and choice['sn']:
        # find dust mass from TF01 in dust_mass_sn table
        # assume massive star winds don't form dust
        dustmass = reduce_sn**-1*find_nearest(np.array(dust_mass_sn),m)[1]
    else:
        dustmass = 0.
    return dustmass

def grow_timescale(on,e,g,sfr,z,d):
    '''
    Calculates the grain growth timescale in years
    Based on Mattsson & Andersen 2012 (MNRAS 423, 38)

    In:
    - e: free parameter with MW value set to ~500-1000
    - G: gas mass in Msolar
    - D: dust mass in Msolar
    - Z: metallicity mass fraction
    - SFR: star formation rate in units of Msolar/Gyr
    - f_c: fraction of gas in cold dense state for grain growth

    '''
    if (on == False or z <= 0 or e <= 0): # set to zero if destroy not turned on
        t_grow = 0
    else:
        t_grow = g/(e*z*sfr)
        t_grow = t_grow/(1-((d/g)/z)) #to account for metals already locked up in grains
    return t_grow #units of Gyrs

def graingrowth(on,e,g,sfr,z,md,f_c):
    '''
    Calculates the grain growth contribution to dust mass, also
    returns grain growth timescale
    Based on Mattsson & Andersen 2012 (MNRAS 423, 38)

    In:
    -- choice: is grain growth turned on or off (True or False)
    -- e: the grain growth epsilon factor given by user
    -- g: gas mass at time t in Msolar
    -- sfr: SFR at time t in Msolar per Gyr
    -- z: metallicity of system (Mz/Mg)
    -- md: dust mass at time t in Msolar
    -- f_c: fraction of gas in cold dense clouds

    e between 500-1000 appropriate for timescales < 1 Gyr.
    In dust evolution, dMd/dt is proportional to Md/t_grow
    '''

    if (on == False or md == 0 or z == 0 or e == 0): #accounts for 1/z in equation
        mdust_gg = 0.
        time_gg = 0.
    else:
        time_gg = grow_timescale(on,e,g,sfr,z,md) # units of Gyrs as SFR = Msun/Gyr
        mdust_gg = md * f_c * (1.-((md/g)/z)) * time_gg**-1  # units of mdust per Gyr
    return mdust_gg, time_gg # units of Gyrs

def destruction_timescale(on,destruct,g,supernova_rate):
    '''
    Calculates the dust destruction timescale in years
    Based on Dwek, Galliano & Jones 2004 (ApJ, 662, 927)

    In
    --  destruct: the amount of gas mass cleared by each supernova event
    --  g: gas mass in Msolar
    -- sn_rate: supernova rate in N/yr

    destruct is often 100 or 1000 Msolar, appropriate for SNe expanding
    into galactic densities of 1cm^-3 or 0.1cm^-3 respectively.
    '''
    if (supernova_rate <= 0 or on == False or destruct == 0): # set to zero if destroy not turned on
        t_destroy = 0.
    else:
        # sn_rate is in units of N per Gyr
        t_destroy = g/(destruct*supernova_rate)  # units are in Gyrs
    return t_destroy # units are in Gyrs

def destroy_dust(on,destruct,gasmass,supernova_rate,md,f_c):
    '''
    Determine how much dust mass is removed by destruction in SN shocks
    Calls destruction_timescale function

    In:
    -- choice: is destruction turned on or off (True or False)
    -- destruct: value of destruction parameter MISM
    -- gasmass: gas mass of system in Msolar at time t
    -- supernova_rate: rate of core-collapse SN at time t (in units of Gyr^-1)
    -- md: dust mass at time t
    -- f_c: fraction of gas in cold dense clouds

    In dust evolution, dMd/dt is proportional to (1-cold fraction) * Md/t_destroy
    '''

    if (on == False or md == 0 or supernova_rate == 0 or destruct == 0):
        mdust_des = 0
        t_des = 0
    else:
        t_des = destruction_timescale(on,destruct,gasmass,supernova_rate)
        mdust_des = md*(1-f_c)*t_des**-1
    #print t_des, mdust_des
    return mdust_des, t_des # in Gyrs

def inflows(sfr,parameter):
    '''
    Define inflow rate, parameterised by N x SFR
    See Rowlands et al 2014 (MNRAS 441 1040)

    In:
    -- sfr: SFR at time t
    -- parameter: inflow parameter defined in input dictionary
    '''
    inflow_rate = sfr*parameter
    return inflow_rate

def gas_inandout(in_on,out_on,in_sfr,sfr,m):
    '''
    Derive the gas lost and gained from inflows and outflows

    In:
    -- in_on: are inflows turned on? (True/False)
    -- out_on: are outflows turned on? (True/False)
    -- in_sfr: inflow rate at time t parameterised by N x SFR (See Rowlands et al 2014 MNRAS 441 1040)
    -- sfr: SFR at time t
    -- m: stellar mass at time t
    '''
    if not in_on:
        gas_inf = 0
    else:
        gas_inf = inflows(sfr,in_sfr)
    if not out_on:
        gas_out = 0.
    else:
        gas_out = outflows_feldmann(sfr, m)
    return gas_inf,gas_out

def metals_inandout(in_on,in_sfr,in_met,out_on,out_met,sfr,Z,m):
    '''
    Derive the metals lost and gained from inflows and outflows

    In:
    -- in_on: are inflows turned on? (True/False)
    -- in_sfr: inflow rate at time t parameterised by N x SFR (See Rowlands et al 2014 MNRAS 441 1040)
    -- in_met: the metallicity of the inflow gas
    -- out_on: are outflows turned on? (True/False)
    -- out_met: is the outflow gas enriched? (True/False)
    -- sfr: SFR at time t
    -- Z: value of metallicity of system at time t
    -- m: stellar mass at time t
    '''
    if in_on == False:
        metal_inf = 0.
    else:
        metal_inf = in_met*inflows(sfr, in_sfr)

    if out_on == False or out_met == False:
        metal_out = 0.
    else:
        metal_out = Z*outflows_feldmann(sfr, m)
    return metal_inf,metal_out

def dust_inandout(in_on,in_sfr,in_md,out_on,out_md,sfr,D,m):
    '''
    Derive the dust mass lost and gained from inflows and outflows

    In:
    -- in_on: are inflows turned on? (True/False)
    -- in_sfr: inflow rate at time t parameterised by N x SFR (See Rowlands et al 2014 MNRAS 441 1040)
    -- in_md: the dust-to-gas ratio of the inflow gas
    -- out_on: are outflows turned on? (True/False)
    -- out_md: does the outflow include dust (True or False)
    -- sfr: SFR at time t
    -- D : dust-to-gas ratio of system at time t
    -- m: stellar mass at time t
    '''
    # If inflows set to False, dust gained in inflows = 0
    if in_on == False:
        dust_inf = 0.
    else:
        dust_inf = in_md*inflows(sfr,in_sfr)
    # if outflows: False or dust outflows False, dust lost in outflows = 0
    if (out_on == False or out_md == False):
        dust_out = 0.
    else:
        dust_out = D*outflows_feldmann(sfr,m)
    return dust_inf,dust_out

def outflows_feldmann(sfr,m):
    '''
    Define outflow rate, parameterised by epsilon_out = 2*f_comb, outflows = epsilon_out x SFR
    See Feldmann et al 2015 MNRAS 449 327 Eq 27 derived from Hopkins et al 2012 MNRAS 421 3522 (Fig 7)
    Here we use same terminology as their paper

    In:
    -- sfr: SFR at time t
    -- m: stellar mass at time t
    '''
    x = 1
    # the function is based on simulations in Hopkins et al 2012 and only go to logM* = 8.1 so we truncate mstar here
    m_low = 1e8
    if (m < m_low):
        outflow_feld = 0.
    else:
        # equation from Feldmann et al 2015 based on Hopkins et al 2012 simulations
        y = (m/1e10)**-0.59
        f_comb = (x+y) - (x**-1+y**-1)**-1
        epsilon_out = 2 * f_comb
        # Feldmann doesnt set a max value of outflow, but Hopkins paper has mean < Feldmann for low M* so
        # Feldmann equation likely overestimates the outflow
        # Here we use Hopkins mean + standard deviation to limit the max outflow rate to be < 30
        # (there are sims above this, but only few outliers)
        if (epsilon_out > 30):
            epsilon_out = 30
        # Feldmann sets outflows to never be less than 2 x SFR, but Hopkins paper has floor at ~1
        elif (epsilon_out < 1):
            epsilon_out = 1
        outflow_feld = sfr * epsilon_out
    return outflow_feld

def mass_integral(choice, delta_lims, reduce_sn, t, metallicity, sfr_lookup, z_lookup, oxy_lookup, imf):
     '''
     This function does the mass integral for:
     - e(t): ejected gas mass em
     - ez(t): ejected metal mass ezm
     - ed(t): ejected dust mass edm
     - Zdiff and SFR diff are also calculated ie at t-taum (when stars that are dying now were born)

     In:
     -- choice: array of dust source choices
     -- delta_lims: efficiency of fresh metals condensing into dust, set by user
     -- t: time in Gyrs
     -- metallicity: metal mass fraction Mz/Mg
     -- sfr_lookup: SFR array (time, SFR) based on previous time steps
     -- z_lookup: metallicity array (time, Z) based on previous time steps
     -- imf: choice of IMF
     '''
     mu = t_lifetime[-1]['mass']

     t_0 = 1e-3
     ezm = 0.
     eom= 0.
     edm = 0.
     em = 0.

     lifetime_cols = {'low_metals':1, 'high_metals':2}
     if metallicity <= 0.008:  # in between Z=0.001 and Z=0.02 files from Schaller et al 1992
        # get the correct lower mass of integral based on age of system
         m_min = lookup_fn(t_lifetime,'lifetime_low_metals',t)['mass']
         # define column choice to speed up taum fn below
         col_choice = lifetime_cols['low_metals']
     else:
         m_min = lookup_fn(t_lifetime,'lifetime_high_metals',t)['mass']
         # define column choice to speed up taum fn below
         col_choice = lifetime_cols['high_metals']

     # checks that m_min does not exceed mu
     if(m_min >= mu):
         m_min = mu

     # increasing the number of steps increases the
     # resolution in the mass integral
     steps = 500
     m = m_min
     dlogm = 0
     logmnew = np.log10(m) + dlogm
     dm = 10**(logmnew)- m
     dlogm = (np.log10(mu)-np.log10(m_min))/steps

     count = 0

     z_near = lambda td : z_lookup[(abs(z_lookup[:,0]-td)).argmin()]
     oxy_near = lambda td : oxy_lookup[(abs(oxy_lookup[:,0]-td)).argmin()]
     sfr_near = lambda td : sfr_lookup[(abs(sfr_lookup[:,0]-td)).argmin()]

     # loop over the full mass range
     while count < steps:
         count += 1
         logmnew = np.log10(m) + dlogm
         dm = 10.0**(logmnew) - m
         mmid = 10.0**((logmnew+np.log10(m))/2.0)

         # pull out lifetime of star of mass m so we can
         # calculate SFR when star was born which is t-lifetime
         taum = lookup_taum(mmid,col_choice)
         tdiff = t - taum
         # only release material after stars die
         if tdiff <= 0:
             ezm = 0
             eom = 0
             edm = 0
             em = 0
         else:
             # get nearest Z and SFR which corresponds to Z and SFR at time=t-taum
             zdiff = z_near(tdiff)[1] # find_nearest(z_lookup,tdiff)[1]
             oxydiff = oxy_near(tdiff)[1] # find_nearest(oxy_lookup,tdiff)[1]
             sfrdiff = sfr_near(tdiff)[1] # find_nearest(sfr_lookup,tdiff)[1]
             ezm += ejected_metal_mass(mmid, sfrdiff, zdiff, metallicity, imf) * dm
             eom += ejected_oxygen_mass(mmid, sfrdiff, oxydiff, metallicity, imf) * dm
             em += ejected_gas_mass(mmid, sfrdiff, imf) * dm
             edm += ejected_dust_mass(choice, delta_lims, reduce_sn, mmid, sfrdiff, zdiff, metallicity, imf) * dm

         #Calculate the next mass value
         mnew = 10**(logmnew)
         m = mnew
     return em, ezm, eom, edm
