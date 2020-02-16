'''
Chemevol - Python package to read in a star formation history file,
input galaxy parameters and run chemical evolution to determine the evolution
of gas, metals and dust in galaxies.

The code is based on Morgan & Edmunds 2003 (MNRAS, 343, 427)
and described in detail in De Vis et al 2017, 2020 (MNRAS, 471, 1743; ).

If you use this code, please do cite the above papers.

Copyright (C) 2015 Haley Gomez, Edward Gomez, Pieter De Vis and Simon Schofield, Cardiff University and LCOGT
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
import math
from numpy import abs, array
import logging
from .lookups import find_nearest, find_yield, find_tot_yield, find_mass_loading, find_reaccretion_time, \
                        dust_mass_sn, dust_mass_sn_eff, t_lifetime, lookup_fn, lookup_taum

logger = logging.Logger('chem')

def extra_sfh_and_inflows(sfh, gamma, tstart):
    '''
    This extrapolates the SFH/SFE and inflows provided to start at 0.001Gyr
    with 100 extra steps between 0.001Gyr and the first non-zero
    entry in the list.
    In:
    -- sfh: combined array time, SFH/SFE and inflows
    -- gamma: parameter that lowers the SHF/SFE and inflows according to power law 

    Returns a new SFH/SFE and inflows list made from joining the
    extrapolated SFH/SFE in this routine to the original input file
    '''
    #to start integral at t_0 regardless of when SFH file starts
    t_0 = tstart # we want it to start at 1e-3
    sfh = sfh[np.where(sfh[:,0]>tstart)]
    tend_sfh=sfh[1][0] # 1st time array after 0
     
    # work out difference between t_0 and [1] entry in SFH
    dlogt = (np.log10(tend_sfh) - np.log10(t_0))/150.  # use larger number instead of 150 if the SFH at the start is not smooth
    norm = sfh[1][1]*(1./np.exp(-1.*gamma*(tend_sfh-tstart)))
    norm_inflow = sfh[1][2]*(1./np.exp(-1.*gamma*(tend_sfh-tstart)))
    t_new = t_0
    newSFH = []
    newinflows = []
    #create new array between 0.001 Gyr and start of SFH data
    while t_new < tend_sfh:
        sfr_new = norm * np.exp(-1.*gamma*(t_new-tstart))
        inflow_new = norm_inflow * np.exp(-1.*gamma*(t_new-tstart))
        newSFH.append([t_new,sfr_new])
        newinflows.append([t_new,inflow_new])
        t_new = 10.**(np.log10(t_new)+dlogt)

    for i in range(1,len(sfh)-1):
        newSFH.append([sfh[i+1][0],sfh[i+1][1]])
        newinflows.append([sfh[i+1][0],sfh[i+1][2]])   
    return newSFH[1::],newinflows[1::]

def astration(parameter,gasmass,sfr):
    '''
    Calculates how much of a given parameter is consumed in the process of star formation.
    E.g. the dust that is mixed with the gas going into a new star will be destroyed.
    '''
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
    '''
    Chabrier (2003, PASP 115 763) IMF
    '''
    if m <= 1.0:
        imf = np.exp(-1.*(np.log10(m)+1.1023729) * (np.log10(m)+1.1023729))
        imf = (0.85*imf)/0.952199/m
    else:
        imf = 0.24*(m**-1.3)/m
    return imf

def imf_topchab(m):
    '''
    Chabrier (2003) IMF that was modified to be more top-heavy.
    It you want to change the slope to be -0.5, you need to change norm factor by 4.72424.
    '''
    if m <= 1.0:
        imf = np.exp(-1.*(np.log10(m)-np.log10(0.079))**2.)
        imf = imf*(0.85/2.21896)/((2.*0.69**2.))/m
    else:
        imf = 0.1081587*(m**-0.8)/m
    return imf

def imf_kroup(m):
    '''
    Kroupa & Weidner (2003, ApJ 598 1076) IMF 
    '''
    if m <= 0.5:
        imf = 0.58*(m**-0.30)/m
    elif (m > 0.5) & (m <= 1.0):
        imf = 0.31*(m**-1.20)/m
    else:
        imf = 0.31*(m**-1.70)/m
    return imf

def imf_salp(m):
    '''
    Salpeter (1955, ApJ 121 161) IMF
    '''
    imf = (0.17/0.990465)*(m**-1.35)/m
    return imf

def initial_mass_function(m, choice):
    '''
    Returns the IMF for a given choice of function and stellar mass.

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

def ejected_metal_masses(t,m, sfrdiff, zdiffs, metallicity, imf, SNyield, AGByield, nisotopes, isotopes, totyields):
    '''
    Calculate the ejected metal masses in (total,Oxygen,Nitrogen) from stars by mass loss/stellar death
    at time t, needs to be integrated from mass corresponding to
    age of system (tau(m)) -- 120 Msolar.

    In:
     -- t: time in Gyrs
     -- m: (initial) mass of star
     -- sfrdiff: SFR at the time the star was formed
     -- zdiffs: metallicity at the time the star was formed
     -- metallicity: metal mass fraction Mz/Mg
     -- imf: choice of IMF
     -- SNyield: string identifier for which metal yield tables is used for SN, set by user
     -- AGByield: string identifier for which metal yield tables is used for AGB, set by user
     -- nisotopes: number of different isotopes/metal budgets that are being tracked (typically 3 for Z, O, N)
     -- isotopes: identifyers for the different isotopes/metal budgets (typically Z, O, N for total metal budget, oxygen and nitrogen)
     -- totyields: boolean whether total metal yield tables are used (True) or fresh metal yields are used (False)

    
    When using total metal yield tables:
    dej(m,t) = mp(m,Z) x SFR(t-tau(m)) x phi(m)
    When using fresh metal yield tables:
    dej(m,t) = (m-m_R(m)*Z(t-taum) + mp(m,Z)) x SFR(t-tau(m)) x phi(m)
    '''


    if m >= 120.0:
        dej = [0.0]*nisotopes

    else:
        if totyields:
            qZs=find_tot_yield(m,metallicity, SNyield, AGByield, nisotopes, isotopes)   # use total yields
            dej = [(m*qZs[iso]) * sfrdiff * imf(m) for iso in range(nisotopes)]

        else:    
            qZs=find_yield(m,metallicity, SNyield, AGByield, nisotopes, isotopes) # use fresh yields
            dej = [((m - remnant_mass(m))*zdiffs[iso] + m*qZs[iso]) * sfrdiff * imf(m) for iso in range(nisotopes)] # add fresh yields to pre-existing metals
        for i in range(len(dej)):
            if dej[i]<0:
                dej[i]=0
 
    return dej

def ejected_dust_mass(choice, delta_lims, reduce_sn, m, mz_ej):
    '''
    Calculate the ejected dust mass from stars by mass loss/stellar death
    at time t, needs to be integrated from mass corresponding to
    age of system (tau(m)) -- 120 Msolar

    Includes contributions from both pre-existing and newly formed metals (together included in mz_ej).
    Dust from massive stars (in SN only) are from Todini & Ferrara 2001 (TF01), converted to efficiencies.
    For both SN and LIMS, a fixed fraction of the expelled metals are converted into dust.

    In:
    -- choice: array of dust source choices, set by user in inits
               0th element = sn value 1 or 0
               1st element = lims value 1 or 0
               2nd element = grain growth value 1 or 0
    -- delta_lims: efficiency of fresh metals condensing into dust in lims
    -- reduce_sn: factor to reduce SN dust contribution, set by user in inits
    -- m: (initial) mass of star
    -- mz_ej: the mass of metals expelled by the star

    '''
    # If LIMS is turned on or off
    if m <= 8.:
        if choice['lims']:
            dej= delta_lims*  mz_ej     

        else:
            dej= 0.

    elif m > 40.: # no dust from stars with m>40Msun.
        dej = 0.0

    else: 
        if choice['sn']:
            dej = reduce_sn**-1*find_nearest(np.array(dust_mass_sn_eff),m)[1]*  mz_ej  
            #dej = reduce_sn**-1*find_nearest(np.array(dust_mass_sn),m)[1] # alternative to work with total dust yields rather than relative to the metal yields.
        else:
            dej =0.0

    return dej

def grow_timescale(on,e,g,sfr,z,d,f_available):
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
    - f_available: fraction of metals that are available for dust grain growth

    '''
    if (on == False or z <= 0 or e <= 0): # set to zero if destroy not turned on
        t_grow = 0
    else:
        t_grow = g/(e*z*sfr)
        t_grow = t_grow/(1.-((d/g)/z)/f_available) #to account for metals already locked up in grains and metals that can never be accreted onto dust grains
 
    return t_grow #units of Gyrs


def graingrowth(on,e,g,sfr,z,md,f_c,f_available):
    '''
    Calculates the grain growth contribution to dust mass, also
    returns grain growth timescale
    Based on Mattsson & Andersen 2012 (MNRAS 423, 38)

    In:
    -- on: is grain growth turned on or off (True or False)
    -- e: the grain growth epsilon factor given by user
    -- g: gas mass at time t in Msolar
    -- sfr: SFR at time t in Msolar per Gyr
    -- z: metallicity of system (Mz/Mg)
    -- md: dust mass at time t in Msolar
    -- f_c: fraction of gas in cold dense clouds
    - f_available: fraction of metals that are available for dust grain growth

    e between 500-1000 appropriate for timescales < 1 Gyr.
    '''

    if (on == False or md == 0 or z == 0 or e == 0): #accounts for 1/z in equation
        mdust_gg = 0.
        time_gg = 0.
    else:
        time_gg = grow_timescale(on,e,g,sfr,z,md,f_available) # units of Gyrs as SFR = Msun/Gyr; the factor (1.-((md/g)/z)) present in Mattsson & Andersen 2012 is a typo and should not be included here (it should only be in t_grow calculation above)
        if time_gg>0.:
            mdust_gg = md * f_c * time_gg**-1 # units of mdust per Gyr; 
        else:
            mdust_gg=0.

    return mdust_gg, time_gg # units of Msolar, Gyrs

def graingrowth_THEMIS_cloud(on,kgg,g,sfr,z,md,f_c,f_available):
    '''
    Calculates the grain growth contribution to dust mass in clouds within the THEMIS (Jones et al., 2013, 2017, 2020/in prep.) dust modelling framework, 
    also returns grain growth timescale.

    In:
    -- on: is grain growth turned on or off (True or False)
    -- kgg: the grain growth factor given by user
    -- g: gas mass at time t in Msolar
    -- sfr: SFR at time t in Msolar per Gyr
    -- z: metallicity of system (Mz/Mg)
    -- md: dust mass at time t in Msolar
    -- f_c: fraction of gas in cold dense clouds
    - f_available: fraction of metals that are available for dust grain growth
    '''

    if (on == False or md == 0 or z == 0 or kgg == 0): #accounts for 1/z in equation
        mdust_gg = 0.
        time_gg = 0.
    else:
        mdmz=md/g/z
        mdust_gg =f_c * kgg * md /g * z/0.0134  * sfr *(1-mdmz/f_available)# units of mdust per Gyr; 0.0134 is the MW metallicity; 0.5 is the fraction of the dense cloud-formed a-C:H mantle material that is processed into refractory a-C dust (e.g., Jones & Ysard (2019, submitted)
        time_gg = md * f_c * mdust_gg**-1

    return mdust_gg, time_gg # units of Msolar, Gyrs

def graingrowth_THEMIS_diffuse(on,kgg,g,sfr,z,md,f_c,f_available):
    '''
    Calculates the grain growth contribution to dust mass in the diffuse ISM within the THEMIS (Jones et al., 2013, 2017, 2020/in prep.) dust modelling framework, 
    also returns grain growth timescale.

    In:
    -- on: is grain growth turned on or off (True or False)
    -- kgg: the grain growth factor given by user
    -- g: gas mass at time t in Msolar
    -- sfr: SFR at time t in Msolar per Gyr
    -- z: metallicity of system (Mz/Mg)
    -- md: dust mass at time t in Msolar
    -- f_c: fraction of gas in cold dense clouds
    - f_available: fraction of metals that are available for dust grain growth
    '''

    if (on == False or md == 0 or z == 0 or kgg == 0): #accounts for 1/z in equation
        mdust_gg = 0.
        time_gg = 0.
    else:
        mdmz=md/g/z
        mdust_gg =(1-f_c) * kgg * 5. * md /g * z/0.0134 * (1-mdmz/f_available) *g # units of mdust per Gyr; 0.0134 is the MW metallicity; 0.5 is the fraction of the dense cloud-formed a-C:H mantle material that is processed into refractory a-C dust (e.g., Jones & Ysard (2019, submitted)
        time_gg = md * f_c * mdust_gg**-1

    return mdust_gg, time_gg # units of Msolar, Gyrs    


def destruction_timescale(on,destruct,g,supernova_rate):
    '''
    Calculates the dust destruction timescale for destruction by SN shocks.
    Based on Dwek, Galliano & Jones 2004 (ApJ, 662, 927).

    In
    -- on: is dust descturction by shocks turned on or off (True or False)
    -- destruct: the amount of gas mass cleared by each supernova event
    -- g: gas mass in Msolar
    -- supernova_rate: supernova rate in N/yr

    '''
    if (supernova_rate <= 0 or on == False or destruct == 0): # set to zero if destroy not turned on
        t_destroy = 0.
    else:
        # sn_rate is in units of N per Gyr
        t_destroy = g/(destruct*supernova_rate)  # units are in Gyrs
 
    return t_destroy # units are in Gyrs

def destroy_dust(on,destruct,gasmass,supernova_rate,md,f_c,sn_eff):
    '''
    Determine how much dust mass is removed by destruction in SN shocks
    Calls destruction_timescale function

    In:
    -- on: is destruction turned on or off (True or False)
    -- destruct: value of destruction parameter MISM
    -- gasmass: gas mass of system in Msolar at time t
    -- supernova_rate: rate of core-collapse SN at time t (in units of Gyr^-1)
    -- md: dust mass at time t
    -- f_c: fraction of gas in cold dense clouds
    -- sn_eff: correction factor to obtain the effective SN rate for dust destruction (to account for previous SN clearing out dust in the vicinity)
    '''

    if (on == False or md == 0 or supernova_rate == 0 or destruct == 0):
        mdust_des = 0
        t_des = 0
    else:
        t_des = destruction_timescale(on,destruct,gasmass,sn_eff*supernova_rate) 
        mdust_des = md*(1-f_c)*t_des**-1
  
    return mdust_des, t_des # in Msolar, Gyrs

def destroy_dust_SN_THEMIS(on,tau_d,gasmass,supernova_rate,md,f_c,sn_eff):
    '''
    Determine how much dust mass is removed by destruction in SN shocks by THEMIS. The prescription is essentially the same as for destroy_dust(), 
    except that it is expressed in terms of the amount of dust that would be destroyed per SN in MW conditions. 
    Essentially this means that the destruct parameter is just reduced by a factor of 135 to find tau_d.
    However this allows for calculations consistent with Bocchio et al (2014), which result in destruct tau_d=10 for silicate dust and tau_d=30 for carbonaceous dust.
    The above numbers have alred folded in the correction factor to obtain the effective SN rate for dust destruction (to account for previous SN clearing out dust in the vicinity)
    In:
    -- on: is this destruction term turned on or off (True or False)
    -- tau_d: dust mass cleared per SN for MW conditions
    -- gasmass: gas mass of system in Msolar at time t
    -- supernova_rate: rate of core-collapse SN at time t (in units of Gyr^-1)
    -- md: dust mass at time t
    -- f_c: fraction of gas in cold dense clouds
    '''

    if (on == False or md == 0 or supernova_rate == 0 or tau_d == 0):
        mdust_des = 0
        t_des = 0
    else:
        t_des = gasmass/(135*tau_d*supernova_rate*sn_eff/0.36)    #135 is a normalisation factor to scale to MW conditions (the THEMIS MW gas-to-dust ratio is 135)
        mdust_des = md*(1-f_c)*t_des**-1

    return mdust_des, t_des # in Msolar, Gyrs   

def destroy_dust_frag_THEMIS(on,tau,gasmass,ssfr,md,f_c):
    '''
    Determine how much dust mass is removed in the diffuse ISM by THEMIS (photo-)fragmentation of 
    large a-C:H/a-C grains and the a-C mantles on the amorphous silicate grains (a-Sil/a-C). 
    See De Vis et al (2020) or Jones et al (2020/in prep.).

    In:
    -- on: is this destruction term turned on or off (True or False)
    -- tau_d: dust mass cleared per SN for MW conditions
    -- gasmass: gas mass of system in Msolar at time t
    -- supernova_rate: rate of core-collapse SN at time t (in units of Gyr^-1)
    -- md: dust mass at time t
    -- f_c: fraction of gas in cold dense clouds
    '''

    if (on == False or md == 0 or tau == 0):
        mdust_des = 0
        t_des = 0
    else:
        mdust_des = (1-f_c) * (1-0.1) * tau * md * ssfr/0.027  # (1-0.1) factor is to account for fragmentation not working on silicate grains
        t_des = md * (1-f_c) * (1-0.1) * mdust_des**-1

    return mdust_des, t_des # in Msolar, Gyrs        

def inflows_SFR(sfr,parameter):
    '''
    Outdated function that is not used in current model but could be of use to others wanting to experiment other inflow prescriptions
    Define inflow rate, parameterised by N x SFR
    See Rowlands et al 2014 (MNRAS 441 1040)

    In:
    -- sfr: SFR at time t
    -- parameter: inflow parameter defined in input dictionary
    '''
    inflow_rate = sfr*parameter
    return inflow_rate

def gas_inandout(redshift,in_on,out_on,in_sfr,sfr,m,reduceout):
    '''
    Derive the gas lost and gained from inflows and outflows

    In:
    -- in_on: are inflows turned on? (True/False)
    -- out_on: are outflows turned on? (True/False)
    -- in_sfr: inflow rate at time t
    -- sfr: SFR at time t
    -- m: stellar mass at time t
    -- reduceout: reduce outflows by a fixed amount
    '''
    if not in_on:
        gas_inf = 0
    else:
        gas_inf = in_sfr
        #gas_inf = inflows_SFR(sfr,in_sfr)  # alternative inflow prescription (inputs need to be edited accordingly)
    if not out_on:
        gas_outflows=0,0,0
        gas_out = 0.
    else:
        gas_outflows= np.array(outflows(redshift,sfr, m))/reduceout
        gas_out = gas_outflows[0]
 
    return gas_inf,gas_out,np.asarray(gas_outflows[1::])

def metals_inandout(in_met,out_met,metallicities,inflow_metalfractions,gas_in,gas_out,nisotopes):
    '''
    Derive the metals (from different isotopes) lost and gained from inflows and outflows

    In:
    -- in_met: the metallicity of the inflow gas
    -- out_met: is the outflow gas enriched? (True/False)
    -- metallicities: the current metallicities of each isotope
    -- inflow_metalfractions: the mass fractions of each isotype, relative to the total metal mass
    -- gas_in: gas inflow rate
    -- gas_out: gas outflow rate
    -- nisotopes: number of isotopes under study
    '''
    
    metals_inf = [inflow_metalfractions[iso]*in_met*gas_in for iso in range(nisotopes)]
    if out_met == False:
        metals_out = np.zeros(nisotopes)
    else:
        metals_out= [metallicities[iso]*gas_out for iso in range(nisotopes)]
 
    return metals_inf,metals_out

def dust_inandout(in_md,out_md,D,gas_in,gas_out):
    '''
    Derive the dust mass lost and gained from inflows and outflows

   In:
    -- in_md: the dust-to-gas ratio of the inflow gas
    -- out_md: does the outflow include dust (True/False)
    -- D : dust-to-gas ratio of system at time t
    -- gas_in: gas inflow rate
    -- gas_out: gas outflow rate
    '''

    dust_inf = in_md*gas_in
    if out_md == False:
        dust_out = 0.
    else:
        dust_out = D*gas_out

    return dust_inf,dust_out

def outflows(redshift,sfr,m):
    '''
    Define outflow rate, parameterised by the mass loading factors epsilon_out, outflows = epsilon_out x SFR.
    The mass loading factors are taken from Nelson et al 2019 (https://arxiv.org/pdf/1902.05554) for 5 different redshift.
    For each mass and redshift, the mass loading factors for three different outflow velocity bins are given (v>0 km/s,v>150 km/s,v>300 km/s). 

    In:
    -- redshift: redshift at time t
    -- sfr: SFR at time t
    -- m: stellar mass at time t
    '''
    
    epsilon_0,epsilon_150,epsilon_300 = find_mass_loading(m,redshift)
    outflows_0,outflows_150,outflows_300= sfr*10**epsilon_0,sfr*10**epsilon_150,sfr*10**epsilon_300
 
    return outflows_0+outflows_150+outflows_300,outflows_0,outflows_150,outflows_300

def outflows_feldmann(sfr,m):
    '''
    Outdated function that is not used in current model but could be of use to others wanting to experiment other outflow prescriptions
    Define outflow rate, parameterised by epsilon_out = 2*f_comb, outflows = epsilon_out x SFR.
    See Feldmann et al 2015 MNRAS 449 327 Eq 27 derived from Hopkins et al 2012 MNRAS 421 3522 (Fig 7).
    Here we use same terminology as their paper.

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



def mass_integral(choice, delta_lims, reduce_sn, t, metallicity, sfr_lookup, z_lookup, imf, SNyield, AGByield, totyields, nisotopes, isotopes):
     '''
     This function does the mass integral for:
     - e(t): ejected gas mass em
     - ez(t): ejected metal mass ezm
     - ed(t): ejected dust mass edm
     - Zdiff and SFR diff are also calculated ie at t-taum (when stars that are dying now were born)

     In:
     -- choice: array of dust source choices
     -- delta_lims: efficiency of AGB metals condensing into dust, set by user
     -- reduce_sn: reduction factor for the efficiency of SN metals condensing into dust, set by user
     -- t: time in Gyrs
     -- metallicity: metal mass fraction Mz/Mg
     -- sfr_lookup: SFR array (time, SFR) based on previous time steps
     -- z_lookup: metallicity array (time, Z) based on previous time steps
     -- imf: choice of IMF
     -- SNyield: string identifier for which metal yield tables is used for SN, set by user
     -- AGByield: string identifier for which metal yield tables is used for AGB, set by user
     -- totyields: boolean whether total metal yield tables are used (True) or fresh metal yields are used (False)
     -- nisotopes: number of different isotopes/metal budgets that are being tracked (typically 3 for Z, O, N)
     -- isotopes: identifyers for the different isotopes/metal budgets (typically Z, O, N for total metal budget, oxygen and nitrogen)
     '''
     mu = t_lifetime[-1]['mass'] # stellar lifetime

     t_0 = 1e-3
     ezm = np.zeros(nisotopes)  #initialise ejected metal masses which include the different isotopes
     edm = 0.                   #initialise ejected dust mass
     em = 0.                    #initialise ejected gas mass

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

     zs_near = lambda td : z_lookup[(abs(z_lookup[:,0]-td)).argmin()] 
     sfr_near = lambda td : sfr_lookup[(abs(sfr_lookup[:,0]-td)).argmin()]
     # loop over the full mass range
     while count < steps:
         count += 1
         logmnew = np.log10(m) + dlogm
         dm = 10.0**(logmnew) - m
         if dm<0: dm=0
         mmid = 10.0**((logmnew+np.log10(m))/2.0)

         # pull out lifetime of star of mass m so we can
         # calculate SFR when star was born which is t-lifetime
         taum = lookup_taum(mmid,col_choice)
         tdiff = t - taum
         # only release material after stars die
         if tdiff <= 0:
             ezm = np.zeros(nisotopes)
             edm = 0
             em = 0
             dezm = np.zeros(nisotopes)
         else:
             # get nearest Z (and O, N, ...) and SFR which corresponds to Z and SFR at time=t-taum
             zdiffs = [zs_near(tdiff)[1+iso] for iso in range(nisotopes)] # find_nearest(z_lookup,tdiff)[1]
             sfrdiff = sfr_near(tdiff)[1] # find_nearest(sfr_lookup,tdiff)[1]
             dezm = ejected_metal_masses(t,mmid, sfrdiff, zdiffs, metallicity, imf,  SNyield, AGByield, nisotopes, isotopes, totyields)     
             ezm= [ ezm[iso]+dezm[iso]* dm for iso in range(nisotopes)]
             em += ejected_gas_mass(mmid, sfrdiff, imf) * dm
             edm += ejected_dust_mass(choice, delta_lims, reduce_sn, mmid, dezm[0]) * dm
         #Calculate the next mass value
         mnew = 10**(logmnew)
         m = mnew
  
     return em, ezm, edm


def recycle(t, dt, redshift, mstars, time, recycle_gas, recycle_dust, recycle_Z, outflows, gas_out, mdust_out, metals_out, escape_probability_perGyr, recycle_time_factor, nisotopes):
    '''
    This function determines when the current outflows (gas, metals and dust) will be recycled and adds the appropriate inflows to the recycled inflows array at that time.
    The arrays recycle_gas, recycle_dust, recycle_Z contain these recycled inflows and are thus modified by this function.
    These arrays contain the total amount of gas to be recycled rather than the rates (i.e. integrated by dt), so that there is conservation of mass even when timesteps are not constant.
    Each of the outflow components (with different outflow velocities) will be recycled at a different time based on their median outflow velocity. 

    In:
    -- t: time in Gyrs
    -- dt: length of timestep in Gyrs
    -- redshift: redshift
    -- mstars: current stellar mass
    -- time: array containing all the timesteps in the model
    -- recycle_gas: input gas recycling flows array before it is modified
    -- recycle_dust: input dust recycling flows array before it is modified
    -- recycle_Z: input metals recycling flows array before it is modified
    -- outflows: current gas outflow rates of each velocity component
    -- gas_out: current total gas outflow rate
    -- mdust_out: current total dust outflow rate
    -- metals_out: current total metal outflow rate
    -- escape_probability_perGyr: probability that a given mass of outflowing material escapes the gravitational well per Gyr spend in the IGM
    -- recycle_time_factor: scaling factor f_{rec} to scale the precomputed recycling times up or down.
    -- nisotopes: number of different isotopes/metal budgets that are being tracked (typically 3 for Z, O, N)
  
    '''

    # lookup precomputed recycling times for each velocity component and modify by factor
    tff_0,tff_150,tff_300=recycle_time_factor*np.asarray(find_reaccretion_time(mstars, redshift)) 
    # how much of the outflows escape the gravitational well
    recycle_fraction_0,recycle_fraction_150,recycle_fraction_300=[np.exp(-tff*escape_probability_perGyr) for tff in [tff_0,tff_150,tff_300]] 
    
    # for the lowest velocity component
    idtcurrent=(np.abs(time-t)).argmin()
    if np.isfinite(tff_0):
        ttest=tff_0
        idt=(np.abs(time-ttest-t)).argmin() #identify timestep that gas will be recycled
        factors=np.zeros(len(recycle_gas))
        # spread out the recycling of this gas so they do not come in all at the same timestep but are rather smeared out over a range of steps 
        # (a discrete burst in outflow will not result in a discrete burst in recycling).
        for ii in range(int(math.floor(idt-(idt-idtcurrent)/2)),int(min(idt+(idt-idtcurrent)/2,len(recycle_gas)))):
            if time[ii]>t: # outflows need to be recycled after they are expelled 
                factors[ii]= np.exp(-(time[ii]-time[idt])**2/ttest**2/4./2.)

        factors=factors/np.sum(factors) # normalise these factors so that mass is conserved

        # add the appropriate amount of material to the appropriate timesteps     
        for ii in range(int(math.floor(idt/2)),min(idt*2,len(recycle_gas))):
            recycle_gas[ii]+=factors[ii]*outflows[0]*dt*recycle_fraction_0
            recycle_dust[ii]+=factors[ii]*outflows[0]*dt*mdust_out/gas_out*recycle_fraction_0
            recycle_Z[ii]+=factors[ii]*outflows[0]*dt*metals_out/gas_out*recycle_fraction_0

    # repeat for other two components        
    if np.isfinite(tff_150):
        ttest=tff_150
        idt=(np.abs(time-ttest-t)).argmin()
        factors=np.zeros(len(recycle_gas))
        for ii in range(int(math.floor(idt-(idt-idtcurrent)/2)),int(min(idt+(idt-idtcurrent)/2,len(recycle_gas)))):
            factors[ii]= np.exp(-(time[ii]-time[idt])**2/ttest**2/4./2.)
        factors=factors/np.sum(factors)    

        for ii in range(int(math.floor(idt/2)),min(idt*2,len(recycle_gas))):
            recycle_gas[ii]+=factors[ii]*outflows[1]*dt*recycle_fraction_150
            recycle_dust[ii]+=0. # dust in fast outflow components is destroyed by shocks 
            recycle_Z[ii]+=factors[ii]*outflows[1]*dt*metals_out/gas_out*recycle_fraction_150
    if np.isfinite(tff_300):
        ttest=tff_300
        idt=(np.abs(time-ttest-t)).argmin()
        factors=np.zeros(len(recycle_gas))
        for ii in range(int(math.floor(idt-(idt-idtcurrent)/2)),int(min(idt+(idt-idtcurrent)/2,len(recycle_gas)))):
            factors[ii]= np.exp(-(time[ii]-time[idt])**2/ttest**2/4./2.)
        factors=factors/np.sum(factors)    

        for ii in range(int(math.floor(idt/2)),min(idt*2,len(recycle_gas))):
            recycle_gas[ii]+=factors[ii]*outflows[2]*dt*recycle_fraction_300
            recycle_dust[ii]+=0. # dust in fast outflow components is destroyed by shocks 
            recycle_Z[ii]+=factors[ii]*outflows[2]*dt*metals_out/gas_out*recycle_fraction_300
            
    return recycle_gas, recycle_dust, recycle_Z
