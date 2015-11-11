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
        #f.validate_initial_dict(inputs)
        try:
            self.gasmass_init = inputs['gasmass_init']
            self.gamma = inputs['gamma']
            self.imf_type = inputs['IMF_fn']
            self.dust_source = inputs['dust_source']
            self.destroy = inputs['destroy']
            self.inflows = inputs['inflows']
            self.outflows = inputs['outflows']
            self.SFH_file = inputs['SFH']
            self.coldfraction = inputs['cold_gas_fraction']
            self.epsilon = inputs['epsilon_grain']
            self.destroy_ism = inputs['destruct']
            if not self.SFH_file:
                self.SFH_file = 'Milkyway.sfh'
            self.sfh_file = self.SFH_file
            self.load_sfh()
        except KeyError:
            logger.error('You must provide initial parameters')

    def load_sfh(self):
        try:
            vals = np.loadtxt(self.SFH_file)
            scale = [1e-9,1e9] # this puts time in Gyr and SFR in Msun/Gyr
    #        self.sfh = vals*scale
            sfh = vals*scale
            # extrapolates SFH back to 0.001Gyr using SFH file
            final_sfh = f.extra_sfh(sfh, self.gamma)
            self.sfh = np.array(final_sfh)
        except:
            logger.error("File '%s' will not parse" % self.SFH_file)
            self.sfh = None

    def sfr(self, t):
        try:
            vals = find_nearest(self.sfh,t)

            return vals[1]
        except:
            logger.error("No SFH yet")

    def final_sfr(self, t):
        try:
            vals = find_nearest(self.extra_sfr,t)
            return vals[1]
        except:
            logger.error("No SFH yet")

    def ejected_mass(self, t):
        '''
        Calculates the ejected gas mass from stars em (t)
        for gas mass integral where:

        dmg/dt = -SFR(t) + int em(t)*dm + inflows - outflows

        mass is only ejected after stars die ie when
        t - taum (lifetime of star) > 0
        '''
        #initialize
        mu = t_lifetime[-1]['mass']
        dm = 0.5
        em = 0.
        t_0 = 1e-3
        # we pull out mass corresponding to age of system
        # to get lower limit of integral
        m = lookup_fn(t_lifetime,'lifetime_low_metals',t)['mass']
        # to make taum lookup faster
        lifetime_cols = {'low_metals':1, 'high_metals':2}
        while m <= mu:
            # pull out lifetime of star of mass m so we can
            # calculate SFR when star was born which is t-lifetime
            taum =  taum = lookup_taum(m,lifetime_cols['low_metals'])
            tdiff = t - taum
            # only release metals (ejected_gas_mass) after stars die
            if tdiff <= 0:
                sfr_diff = 0.
            else:
                sfr_diff = self.sfr(tdiff)
            # integral calculation
            em += f.ejected_gas_mass(m, sfr_diff, self.imf_type) * dm
            m += dm
        return em

    def ejected_z_mass(self, t, metallicity):
        '''
        Calculates the ejected metal mass from stars ezm (t)
        for metal mass integral where:
        dmz/dt = -Z * SFR(t) + int ezm*dm + Z,inflows - Z,outflows + Z_i
        '''
        # initialize
        dm = 0.5
        ezm = 0.
        now = datetime.now()
        mu = t_lifetime[-1][0]
        # we pull out mass corresponding to age of system
        # to get lower limit of integral
        m = lookup_fn(t_lifetime,'lifetime_low_metals',t)['mass']
        # to make taum lookup faster
        lifetime_cols = {'low_metals':1, 'high_metals':2}
        while m <= mu:
            # pull out lifetime of star of mass m so we can
            # calculate SFR when star was born which is t-lifetime
            taum =  taum = lookup_taum(m,lifetime_cols['low_metals'])
            tdiff = t - taum
            if tdiff <= 0.:
                zdiff = 0.
                sfr_diff = 0.
            else:
                zdiff = metallicity #needs to be z(t-taum)
                sfr_diff = self.sfr(tdiff)
            ezm += f.ejected_metal_mass(m, sfr_diff, zdiff, self.imf_type) * dm
            m += dm
        #    print m, t, tdiff, metallicity, zdiff, ezm
        return ezm

    def ejected_d_mass(self, t, metallicity):
        '''
        Calculates the ejected dust mass from stars edm (t)
        for dust mass integral where:
        dmd/dt = - md/mg * SFR(t) + int edm*dm + md/mg,inflows - md/mg,outflows
                 + md_graingrowth - md_destroy
        '''
        # initialize
        dm = 0.5
        edm = 0.
        now = datetime.now()
        mu = t_lifetime[-1][0]
        # we pull out mass corresponding to age of system
        # to get lower limit of integral
        m = lookup_fn(t_lifetime,'lifetime_low_metals',t)[0]
        #to make taum lookup faster
        lifetime_cols = {'low_metals':1, 'high_metals':2}
        while m <= mu:
            # pull out lifetime of star of mass m so we can
            # calculate SFR when star was born which is t-lifetime
            taum =  lookup_taum(m,lifetime_cols['low_metals'])
            tdiff = t - taum
            if tdiff <= 0.:
                zdiff = 0.
                sfr_diff = 0.
            else:
                zdiff = metallicity #needs to be z(t-taum)
                sfr_diff = self.sfr(tdiff)
            edm += f.ejected_dust_mass(m, sfr_diff, zdiff, self.imf_type) * dm
            m += dm
#            print m, t, tdiff, em
        return edm

    def supernova_rate(self):
        '''
        Calculates the SN rate at time t by integrating over mass m
        '''
        # initialize
        sn_rate_list = []
        dm = 0.5
        prev_t = 1e-3
        # define time array
        time = self.sfh[:,0]
        time = time[time<20.]
        now = datetime.now()
        for t in time:
            # need to clear the sn_rates as we don't want them adding up
            sn_rate = 0.
            dsn_rate = 0.
            if t < 0.049:
                m = lookup_fn(t_lifetime,'lifetime_low_metals',t)[0]
            else:
                m = 9.
            while m < 40.:
                sn_rate += f.initial_mass_function(m, self.imf_type)*dm
                m += dm
            r_sn = self.sfr(t)*sn_rate # units in N per Gyr
            dt = t - prev_t
            prev_t = t
        #    print t, r_sn
            sn_rate_list.append(r_sn)
        print("SN rate exterior loop %s" % str(datetime.now()-now))
        return np.array(sn_rate_list)

    def gas_mass(self):
        '''
        Calculates the gas mass at time t: mg

        for gas mass integral where:
        dmg/dt = -SFR(t) + int em*dm + inflows - outflows
        '''
        mg = self.gasmass_init
        prev_t = 1e-3
        mg_list = []
    # Limit time to less than 20. Gyrs
        time = self.sfh[:,0]
        time = time[time<20.]
    # set up t-taum array
        now = datetime.now()
        for t in time:
        # set up gas lost due to star formation
            gas_ast = self.sfr(t)

        # set up gas inflows using dictionary inputs
            gas_inf = f.inflows(self.sfr(t), self.inflows['xSFR']).value

        # set up gas outflows using dictionary inputs
            gas_out = f.outflows(self.sfr(t), self.outflows['xSFR']).value

        # call ejected gas mass and t-taum
            gas_ej = self.ejected_mass(t)

        # do the integral
            dmg = - gas_ast \
                  + gas_ej \
                  + gas_inf \
                  - gas_out

            dt = t - prev_t
            prev_t = t
            mg += dmg*dt
            mg_list.append(mg)
        #    print t, gas_ast, gas_ej
        print("Gas mass exterior loop %s" % str(datetime.now()-now))
        return time, np.array(mg_list)

    def stellar_mass(self):
        '''
        Calculates the stellar mass at time t where integral

        dmstars/dt = - SFR(t)
        '''
        mstars = 0.
        prev_t = 1e-3
        mstars_list = []
        # Limit time to less than 20. Gyrs
        time = self.sfh[:,0]
        time = time[time<20.]
        for t in time:
            dmstars = self.sfr(t)
            dt = t - prev_t
            prev_t = t
            mstars += dmstars*dt
            mstars_list.append(mstars)
        # Output time and gas mass as Numpy Arrays
        return time, np.array(mstars_list)

    def metal_mass(self,gasmass):
        '''
        Calculates the metal mass at time t: Mz
        for metal mass integral where:
        dmz/dt = -Z* SFR(t) + int ezm*dm + Z,inflows - Z,outflows + Z_i
        '''
        metals = 0.
        prev_t = 1e-3
        metals_pre = 0. # in case of pre-enrichment needed
        metals_list = []
        # Limit time to less than 20. Gyrs
        time = self.sfh[:,0]
        time = time[time<20.]
        z_lookup = []
        now = datetime.now()
        for item, t in enumerate(time):
            mg = gasmass[item]
            metallicity = metals/mg
            z_lookup.append([t,metallicity])
        # set up time, z(t-taum array)
        # set up metals lost due to astration
            metals_ast = metallicity*self.sfr(t)

        #set up metals ejected by LIMS + SNe
            metals_stars = self.ejected_z_mass(t, metallicity)

        # set up metals lost in outflows with outflow metallicity read from input dictionary
            if self.outflows['metals']:
                outflow_metals = metallicity
                metals_out = outflow_metals*f.outflows(self.sfr(t), self.outflows['xSFR']).value
            else:
                metals_out = 0.

        # set up inflows of metals with inflow metallicity from dictionary
            metals_inf = self.inflows['metals']*f.inflows(self.sfr(t), self.inflows['xSFR']).value

            dmetals = - metals_ast \
                      + metals_stars \
                      + metals_pre \
                      + metals_inf \
                      + metals_out

            dt = t - prev_t
            prev_t = t
            metals += dmetals*dt
            metals_list.append(metals)
        # zip array to make metallicity array containing time + Z
        Z = zip(*z_lookup)
        # np.array(Z[1]) is the metallicity
        print("Metal mass exterior loop %s" % str(datetime.now()-now))
        return time, np.array(metals_list), np.array(Z[1])

    def dust_mass(self,gasmass,metallicity,snrate):
        '''
        Calculates the dust mass at time t: md
        for dust mass integral where:

        dmd/dt = - md/mg * SFR(t) + int edm*dm + md/mg,inflows - md/mg,outflows
                 + md_graingrowth - md_destroy

        Returns 3 arrays:
         -- time (Gyrs)
         -- dust_list (Msolar): 3 columns - total Md, Md_stars only, Md_gg only
         -- dz_ratio_list: dust to metal ratio
        '''
        # initialize
        md_all = 0.
        md_stars = 0.
        md_ism = 0.
        prev_t = 1e-3
        dust_list = []
        dz_ratio_list = []
        # Limit time to less than 20. Gyrs
        time = self.sfh[:,0]
        time = time[time<20.]
        # sort out zdiff
        now = datetime.now()
        for item, t in enumerate(time):
            mg = gasmass[item]
            z = metallicity[item]
            r_sn = snrate[item]

        #set up dust mass from stars (recycled(LIMS) + new (SN+LIMS))
            mdust_stars = self.ejected_d_mass(t, z)

        # set up inflow contribution to dust mass (read from dictionary)
            mdust_inf = self.inflows['dust']*f.inflows(self.sfr(t), self.inflows['xSFR']).value

            if self.outflows['dust']:
                mdust_out = (1./mg)*f.outflows(self.sfr(t), self.outflows['xSFR']).value
            else:
                mdust_out = 0.

        #set up removal of dust via destruction parameter
    #        des = 1e-9*f.destruction_timescale(self.destroy_ism,mg,r_sn).value #in Gyrs
    #        if des <= 0:
    #            mdust_des = 0
    #        else:
    #            mdust_des = md*(1-self.coldfraction)/des

        # Integrate dust mass equation with time - want to do this for separate
        # dust sources to plot later
        # Total dust mass
            ddust_all = - md_all*f.astration(mg,self.sfr(t)) \
                        + mdust_stars \
                        + mdust_inf \
                        - mdust_out*md_all \
                        + f.graingrowth(self.epsilon,mg,self.sfr(t),z,md_all,self.coldfraction)  # - mdust_des
        # Dust mass with stars only
            ddust_stars = - md_stars*f.astration(mg,self.sfr(t)) \
                          + mdust_stars \
                          + mdust_inf \
                          - mdust_out*md_stars  #- mdust_des
        # Dust mass with grain growth only
            ddust_ism = - md_ism*f.astration(mg,self.sfr(t)) \
                        + mdust_inf \
                        - mdust_out*md_ism \
                        + f.graingrowth(self.epsilon,mg,self.sfr(t),z,md_ism,self.coldfraction) # - mdust_des
            dt = t - prev_t
            prev_t = t
            md_all += ddust_all*dt
            md_stars += ddust_stars*dt
            md_ism += ddust_ism*dt
            dust_list.append((md_all,md_stars,md_ism))
            if z <= 0.:
                dust_to_metals = 0.
            else:
                dust_to_metals = (md_all/mg)/z
        #    print t, mdust_stars, mdust_inf, mdust_out*md_all,f.graingrowth(self.epsilon,mg,self.sfr(t),z,md_all,self.coldfraction)  # mg/4.8e10, z, r_sn/1e9#, des*1e9/1e6 #, mdust_des #mdust_ast, mdust_stars, des
            dz_ratio_list.append(dust_to_metals)
        print("Dust mass exterior loop %s" % str(datetime.now()-now))
        # Output time and gas mass as Numpy Arrays
        return time, np.array(dust_list), np.array(dz_ratio_list)
