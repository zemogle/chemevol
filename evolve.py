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

import functions as f
import numpy as np
from lookups import find_nearest, lookup_fn, t_lifetime, lookup_taum
import logging
from datetime import datetime

FORMAT = '%(asctime)-15s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger('chem')


class ChemModel:
    def __init__(self, **inputs):
        '''
        set initial parameters from input dictionary and sort out choices
        for IMF, dust source
        '''
        #f.validate_initial_dict(inputs)
        try:
            self.gasmass_init = inputs['gasmass_init']
            self.gamma = inputs['gamma']
            self.tend = inputs['t_end']
            self.imf_type = inputs['IMF_fn']
            self.dust_source = inputs['dust_source']
            self.reduce_sn = inputs['reduce_sn_dust']
            self.destroy = inputs['destroy']
            self.inflows = inputs['inflows']
            self.outflows = inputs['outflows']
            self.SFH_file = inputs['SFH']
            self.coldfraction = inputs['cold_gas_fraction']
            self.epsilon = inputs['epsilon_grain']
            self.destroy_ism = inputs['destruct']
            # check for SFH file or use Milkway.sfh provided
            if not self.SFH_file:
                self.SFH_file = 'Milkyway.sfh'
            self.sfh_file = self.SFH_file
            self.load_sfh()
        except KeyError:
            logger.error('You must provide initial parameters')
        # Set up IMF Function determined by user, allow for variety of spellings
        if (self.imf_type == "Chab" or self.imf_type == "chab" or self.imf_type == "c"):
            self.imf = f.imf_chab
        elif (self.imf_type == "TopChab" or self.imf_type == 'topchab' or self.imf_type == "tc"):
            self.imf = f.imf_topchab
        elif (self.imf_type == "Kroup" or self.imf_type == "kroup" or self.imf_type == "k"):
            self.imf = f.imf_kroup
        elif (self.imf_type == "Salp" or self.imf_type == "salp" or self.imf_type == "s"):
            self.imf = f.imf_salp
        # Declare if destruction on or off
        if self.reduce_sn == False:
            self.reduce_sn = 1
        if (self.destroy == False):
            self.choice_des = 0
        else:
            self.choice_des = 1
        # set up dust source choice from user: 0 = SN dust on, 1 = LIMS dust on, 2 = GG on
        if (self.dust_source == "ALL" or self.dust_source == "all" or self.dust_source == "All"):
            self.choice_dust = (1, 1, 1)
        elif (self.dust_source == "SN" or self.dust_source == "Sn" or self.dust_source == "sn"):
            self.choice_dust = (1, 0, 0)
        elif (self.dust_source == "LIMS" or self.dust_source == "Lims" or self.dust_source == "lims"):
            self.choice_dust = (0, 1, 0)
        elif (self.dust_source == "SN+LIMS" or self.dust_source == "sn+lims" or \
              self.dust_source == "LIMS+SN" or self.dust_source == "lims+sn"):
            self.choice_dust = (1, 1, 0)
        else:
            print ('oops please check the dust sources are in the right format and try again')
            exit()

    def load_sfh(self):
        '''
        takes in input SFH file and extend backwards to start from 1e-3 Gyr
        '''
        try:
            vals = np.loadtxt(self.SFH_file)
            scale = [1e-9,1e9] # Gyr conversions for time, SFR
            sfh = vals*scale # converts time in Gyr and SFR in Msun/Gyr
            # extrapolates SFH back to 0.001Gyr using SFH file and power law (gamma)
            final_sfh = f.extra_sfh(sfh, self.gamma)
            self.sfh = np.array(final_sfh)
        except:
            logger.error("File '%s' will not parse" % self.SFH_file)
            self.sfh = None

    def sfr(self, t):
        '''
        define sfr as function to look up nearest sfr value at any specified time
        '''
        try:
            vals = find_nearest(self.sfh,t)
            return vals[1]
        except:
            logger.error("No SFH yet")

    def gas_metal_dust_mass(self, sn_rate):
        '''
        Calculates the gas, metal and dust mass from stars
        note mass is only ejected after stars die ie when
        t - taum (where taum is lifetime of star) > 0
        '''
        # initialize
        mg = self.gasmass_init
        mstars = 0
        md = 0
        md_all = 0
        md_stars = 0
        md_gg = 0
        metals = 0
        prev_t = 1e-3
        metals_pre = 0
        dust_sources = []
        timescales = []
        mstars_list = []
        z = []
        z_lookup = []
        sfr_list = []
        sfr_lookup = []
        all_results = []
        # Limit time to less than tend
        time = self.sfh[:,0]
        time = time[time < self.tend]
        now = datetime.now()
        # TIME integral
        for item, t in enumerate(time):
            r_sn = sn_rate [item]
            metallicity = metals/mg

            # start appending arrays for needing later
            z.append([t,metallicity])
            z_lookup = np.array(z)
            sfr_list.append([t,self.sfr(t)])
            sfr_lookup = np.array(sfr_list)

            '''
            STARS: dM_stars = sfr(t) * dt
            '''
            dmstars = self.sfr(t)

            '''
            GAS: dMg = (-sfr(t) + e(t) + inflows(t) - outflows(t)) * dt
            set up astration, inflow, outflow components
            '''
            gas_ast = self.sfr(t)
            gas_inf = f.inflows(self.sfr(t), self.inflows['xSFR'])
            gas_out = f.outflows(self.sfr(t), self.outflows['xSFR'])

            '''
            METALS: dMz = (-Z*sfr(t) + ez(t) + Z*inflows(t) - Z*outflows(t)) * dt
            set up astration, inflow and outflow components
            '''
            metals_ast = f.astration(metals,mg,self.sfr(t))
            if self.outflows['metals']:
                metals_out = metallicity*f.outflows(self.sfr(t), self.outflows['xSFR'])
            else:
                metals_out = 0.
            metals_inf = self.inflows['metals']*f.inflows(self.sfr(t), self.inflows['xSFR'])

            '''
            DUST: dMd = (-Md/Mg*sfr(t) + ed(t) + Md/Mg*inflows(t) - Md/Mg*outflows(t)
                         - (1-f)*Md/t_destroy + f(1-Md/Mg)*Md/t_graingrowth) * dt
            set up astration, inflows, outflows, destruction, grain growth components
            '''
            if self.outflows['dust']:
                mdust_out = (md/mg)*f.outflows(self.sfr(t), self.outflows['xSFR'])
            else:
                mdust_out = 0.
            mdust_inf = self.inflows['dust']*f.inflows(self.sfr(t), self.inflows['xSFR'])
            mdust_ast = f.astration(md,mg,self.sfr(t))

            mdust_gg, t_gg = f.graingrowth(self.choice_dust[2], self.epsilon,mg, self.sfr(t), metallicity, md, self.coldfraction)
            mdust_des, t_des = f.destroy_dust(self.choice_des, self.destroy_ism, mg, r_sn, md, self.coldfraction)

            '''
            Get ejected masses from stars when they die
            gas_ej = e(t): ejected gas mass from stars of mass m at t = taum
            metals_stars = ez(t): ejected metal mass from stars of mass m at t = taum (fresh + recycled)
            mdust_stars = ed(t): ejected dust mass from stars of mass m at t = taum (fresh + recycled)
            '''
            gas_ej, metals_stars, mdust_stars = \
                    f.mass_integral(self.choice_dust, self.reduce_sn, t, metallicity, sfr_lookup, z_lookup, self.imf)

            '''
            integrate over time for gas, metals and stars (mg, metals, md)
            '''
            dmg = -gas_ast + gas_ej + gas_inf - gas_out
            dmetals = -metals_ast + metals_stars + metals_pre + metals_inf - metals_out
            ddust = -mdust_ast + mdust_stars + mdust_inf - mdust_out + mdust_gg - mdust_des
            # dust_source_all separates out the dust sources (Md vs t) wihtout including sinks (Astration etc)
            # and grain growth separately (this is the Md vs time contributed by dust sources)
            dust_source_all = mdust_stars + mdust_gg
            dt = t - prev_t             # calculate  next time step
            prev_t = t
            mstars += dmstars*dt
            mg += dmg*dt # gas mass integral
            if mg <= 0:
                # exit program if all ISM removed
                print ('Oops you have no interstellar medium left')
                break
            metals += dmetals*dt # metal mass integral
            md += ddust*dt # dust mass integral
            md_all += dust_source_all*dt # dust mass sources integral
            md_gg += mdust_gg*dt # dust source from grain growth only
            md_stars += mdust_stars*dt # dust source from stars only
            Z = zip(*z_lookup) # write metallicity to an array
            s_f_r = zip(*sfr_lookup) # write SFR lookup array
            dust_sources.append((md_all, md_stars, md_gg)) # write array of dust sources
            timescales.append((t_des,t_gg)) # write array for grain growth & destruction timescales
            if mg <= 0. or metals <=0:  # write dust/metals ratio
                dust_to_metals = 0.
            else:
                dust_to_metals = md/metals
            all_results.append((t,mg,mstars,metals,metallicity,md,dust_to_metals,self.sfr(t)*1e-9))
        print("Gas, metal and dust mass exterior loop %s" % str(datetime.now()-now))
        return np.array(dust_sources), np.array(timescales), np.array(all_results)

    def supernova_rate(self):
        '''
        Calculates the SN rate at time t by integrating over mass m
        '''
        # initialize
        sn_rate_list = []
        dm = 0.01
        prev_t = 1e-3
        # define time array
        time = self.sfh[:,0]
        time = time[time < self.tend]
        for t in time:
            # need to clear the sn_rates as we don't want them adding up
            sn_rate = 0.
            dsn_rate = 0.
            if t < 0.049:
                m = lookup_fn(t_lifetime,'lifetime_high_metals',t)[0]
            else:
                m = 9.
            while m < 40.:
                if m > 10.:
                    dm = 0.5
                sn_rate += f.initial_mass_function(m, self.imf_type)*dm
                m += dm
            r_sn = self.sfr(t)*sn_rate # units in N per Gyr
            dt = t - prev_t
            prev_t = t
            sn_rate_list.append(r_sn)
        return np.array(sn_rate_list)
