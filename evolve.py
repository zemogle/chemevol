from scipy.integrate import quad
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
            final_sfh = self.extra_sfh(sfh)
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
        dmg/dt = -SFR(t) + int em*dm + inflows - outflows
        '''
        mu = t_lifetime[-1]['mass']
        # we pull out mass corresponding to age of system
        # to get lower limit of integral
        m = lookup_fn(t_lifetime,'lifetime_low_metals',t)['mass']
        dm = 0.5
        em = 0.
        lifetime_cols = { 'low_metals':0, 'high_metals':1}
        while m <= mu:
            # pull out lifetime of star of mass m so we can
            # calculate SFR when star was born which is t-lifetime
            taum = lookup_taum(m,lifetime_cols['low_metals'])
            tdiff = t - taum
            if tdiff > 0:
                em += f.ejected_gas_mass(m, self.sfr(tdiff), self.imf_type) * dm
            m += dm
        return em, tdiff

    def ejected_z_mass(self, t, metallicity):
        '''
        Calculates the ejected metal mass from stars ezm (t)
        for metal mass integral where:
        dmz/dt = -Z* SFR(t) + int ezm*dm + Z,inflows - Z,outflows + Z_i
        '''
        now = datetime.now()
        mu = t_lifetime[-1][0]
        # we pull out mass corresponding to age of system
        # to get lower limit of integral
        m = lookup_fn(t_lifetime,'lifetime_low_metals',t)[0]
        dm = 0.5
        ezm = 0.
        while m <= mu:
            # pull out lifetime of star of mass m so we can
            # calculate SFR when star was born which is t-lifetime
            tdiff = find_nearest(self.tdiff,t)[1]
            if tdiff <= 0:
                zdiff =0
            else:
                zdiff = metallicity #needs to be done properly t-taum d
            if tdiff > 0:
                ezm += f.ejected_metal_mass(m, self.sfr(tdiff), zdiff, self.imf_type) * dm
            m += dm
#            print m, t, tdiff, em
        # print("Metal mass interior loop %s" % str(datetime.now()-now))
        return ezm

    def ejected_d_mass(self, t, metallicity):
        '''
        Calculates the ejected dust mass from stars edm (t)
        for dust mass integral where:
        dmd/dt = - md/mg * SFR(t) + int edm*dm + md/mg,inflows - md/mg,outflows
                 + md_graingrowth - md_destroy
        '''
        now = datetime.now()
        mu = t_lifetime[-1][0]
        # we pull out mass corresponding to age of system
        # to get lower limit of integral
        m = lookup_fn(t_lifetime,'lifetime_low_metals',t)[0]
        dm = 0.5
        edm = 0.
        while m <= mu:
            # need to change this for sfr(tdiff) and zdiff when fixed
            edm += f.ejected_dust_mass(m, self.sfr(t), metallicity, self.imf_type) * dm
            m += dm
#            print m, t, tdiff, em
        # print("Metal mass interior loop %s" % str(datetime.now()-now))
        return edm


    def extra_sfh(self, sfh):
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
        norm = sfh[1][1]*(1./np.exp(-1.*self.gamma*tend_sfh))
        sfr_extra = norm * np.exp(-1.*self.gamma*t_0)
        sfr_new = sfr_extra
        t_new = t_0
        newlist = []
        #create new array between 0.001 Gyr and start of SFH data
        while t_new < tend_sfh:
            t_new = 10.**(np.log10(t_new)+dlogt)
            sfr_new = norm * np.exp(-1.*self.gamma*t_new)
            newlist.append([t_new,sfr_new])
        # start from [2:] to account for t[0],t[1] repeated entries
        # when new and in old SFHs combined
        final_sfh = newlist + (sfh.tolist()[2:])
        return final_sfh

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
        t_diff = []
        now = datetime.now()
        for t in time:
            ej, tdiff = self.ejected_mass(t)
            dmg = - self.sfr(t) + ej + f.inflows(self.sfr(t), \
                self.inflows['xSFR']).value + f.outflows(self.sfr(t), self.outflows['xSFR']).value
            dt = t - prev_t
            prev_t = t
            mg += dmg*dt
            mg_list.append(mg)
            t_diff.append([t,tdiff])
        self.tdiff = np.array(t_diff)
        print("Gas mass exterior loop %s" % str(datetime.now()-now))
        # Output time and gas mass as Numpy Arrays
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
        # set up time, z(t-taum array)
            z_lookup.append([t,metallicity])
            tdiff_now = find_nearest(self.tdiff, t)
            zdiff = find_nearest(np.array(z_lookup), tdiff_now[1])[1]

        # set up metals lost due to astration
            metals_ast = metallicity*self.sfr(t)

        #set up metals ejected by LIMS + SNe
            metals_stars = self.ejected_z_mass(t, zdiff)
            
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
            print t, metallicity #mg/4.8e10, metals/mg
        # zip array to make metallicity array containing time + Z
        Z = zip(*z_lookup)
        print("Metal mass exterior loop %s" % str(datetime.now()-now))
        # np.array(Z[1]) is the metallicity
        return time, np.array(metals_list), np.array(Z[1])

    def dust_mass(self,gasmass,metallicity):
        '''
        Calculates the dust mass at time t: md
        for dust mass integral where:
        dmd/dt = - md/mg * SFR(t) + int edm*dm + md/mg,inflows - md/mg,outflows
                 + md_graingrowth - md_destroy
        '''
        md = 0.
        prev_t = 1e-3
        dust_list = []
        dz_ratio_list = []
        # Limit time to less than 20. Gyrs
        time = self.sfh[:,0]
        time = time[time<20.]
        z_time = np.array(zip(time, metallicity))
        for item, t in enumerate(time):
            mg = gasmass[item]
            z = metallicity[item]
        # set up time, z(t-taum) array for use in ejected dust mass integral
            tdiff_now = find_nearest(self.tdiff, t)
            zdiff = find_nearest(z_time, tdiff_now[1])[1]

        #set up dust mass from stars (recycled(LIMS) + new (SN+LIMS))
            mdust_stars = self.ejected_d_mass(t, zdiff)

        # set up inflow contribution to dust mass (read from dictionary)
            mdust_inf = self.inflows['dust']*f.inflows(self.sfr(t), self.inflows['xSFR']).value

        # set up outflow contribution to dust mass (read from input dictionary)
            if self.outflows['dust']:
                outflow_dust = md
                mdust_out = (outflow_dust/mg)*f.outflows(self.sfr(t), self.outflows['xSFR']).value
            else:
                mdust_out = 0.

        #set up removal of dust mass via astration
            mdust_ast = (md/mg)*self.sfr(t)

        #set up contribution of grain growth to dust mass
            gg = 1e-9*f.grow_timescale(self.epsilon,mg,self.sfr(t),z,md).value #in Gyrs
            if gg <= 0:
                mdust_gg = 0.
            else:
                mdust_gg = md*self.coldfraction/gg

        #set up removal of dust via destruction parameter
            #des = 1e-9*f.destruction_timescale(md,mg,self.sfr(t),0.).value
            #mdust_des = md*(1-self.coldfraction)/des

        # integral of dust mass equation with time
        # set dust mass to zero if metallicity = zero
            ddust = - mdust_ast \
                    + mdust_stars \
                    + mdust_inf \
                    - mdust_out \
                    + mdust_gg \
                    - 0.#mdust_des

            dt = t - prev_t
            prev_t = t
            md += ddust*dt
            dust_list.append(md)
            if z <= 0.:
                dust_to_metals = 0.
            else:
                dust_to_metals = (md/mg)/z
            dz_ratio_list.append(dust_to_metals)
        #    print t, mg/4.8e10, z, (md/mg)/z,  #mg/4.8e10, metals/mg

        # Output time and gas mass as Numpy Arrays
    #    return time, np.array(dust_list), np.array(dz_ratio_list)
