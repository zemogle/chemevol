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
            self.tend = inputs['t_end']
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

    def gas_metal_dust_mass(self, sn_rate):
        '''
        Calculates the gas, metal and dust mass from stars
        mass is only ejected after stars die ie when
        t - taum (lifetime of star) > 0
        '''
        # initialize
        mg = self.gasmass_init
        md = 0
        md_all = 0
        md_stars = 0
        md_gg = 0
        metals = 0
        prev_t = 1e-3
        metals_pre = 0
        metals_list = []
        dust_list = []
        dust_list_sources = []
        dz_ratio_list = []
        timescales = []
        mg_list = []
        z = []
        z_lookup = []

        # Limit time to less than tend
        time = self.sfh[:,0]
        time = time[time < self.tend]
        now = datetime.now()
        # TIME integral
        for item, t in enumerate(time):
            r_sn = sn_rate [item]
            metallicity = metals/mg
            z.append([t,metallicity])
            z_lookup = np.array(z)
            # GAS
            # astration, inflows, outflows
            gas_ast = self.sfr(t)
            gas_inf = f.inflows(self.sfr(t), self.inflows['xSFR']).value
            gas_out = f.outflows(self.sfr(t), self.outflows['xSFR']).value

            # METALS
            # astration, inflows, outflows
            metals_ast = metallicity*self.sfr(t)
            metals_inf = self.inflows['metals']*f.inflows(self.sfr(t), self.inflows['xSFR']).value
            if self.outflows['metals']:
                outflow_metals = metallicity
                metals_out = outflow_metals*f.outflows(self.sfr(t), self.outflows['xSFR']).value
            else:
                metals_out = 0.

            # DUST
            # astration, inflows, outflows, grain growth, destruction
            mdust_inf = self.inflows['dust']*f.inflows(self.sfr(t), self.inflows['xSFR']).value

            if self.outflows['dust']:
                mdust_out = (1./mg)*f.outflows(self.sfr(t), self.outflows['xSFR']).value
            else:
                mdust_out = 0.

                # destruction timescales + dust mass from grain growth and destruction
                #  t_des = 1e-6*f.destruction_timescale(self.destroy_ism,mg,r_sn).value
            mdust_gg, t_gg = f.graingrowth(self.epsilon,mg,self.sfr(t),metallicity,md,self.coldfraction)
            mdust_des, t_des = f.destroy_dust(self.destroy_ism,mg,r_sn,md,self.coldfraction)

        # MASS integral for gas, metals, dust
        # initialize
            mu = t_lifetime[-1]['mass']
            dm = 0.01
            t_0 = 1e-3
            ezm = 0.
            edm = 0.
            em = 0.
            # we pull out mass corresponding to age of system
            # to get lower limit of integral
            # to make taum lookup faster
            m = lookup_fn(t_lifetime,'lifetime_low_metals',t)['mass']
            lifetime_cols = {'low_metals':1, 'high_metals':2}
            if metallicity < 0.019:
                col_choice = lifetime_cols['low_metals']
            else:
                col_choice = lifetime_cols['high_metals']
            while m <= mu:

                if m > 10.:
                    dm = 0.5
                # pull out lifetime of star of mass m so we can
                # calculate SFR when star was born which is t-lifetime
                taum = lookup_taum(m,col_choice)
                tdiff = t - taum
                # only release metals (ejected_gas_mass) after stars die
                if tdiff <= 0:
                    sfr_diff = 0.
                    zdiff = 0.
                else:
                    sfr_diff = self.sfr(tdiff)
                    # get nearest Z which corresponds to Z at time=t-taum
                    zdiff = find_nearest(z_lookup,tdiff)[1]
                    ezm += f.ejected_metal_mass(m, sfr_diff, zdiff, metallicity, self.imf_type) * dm
                    em += f.ejected_gas_mass(m, sfr_diff, self.imf_type) * dm
                    edm += f.ejected_dust_mass(m, sfr_diff, zdiff, metallicity, self.imf_type) * dm
                m += dm

            gas_ej = em
            metals_stars = ezm
            mdust_stars = edm

            # do the integral for gas mass with time
            dmg = - gas_ast \
                    + gas_ej \
                    + gas_inf \
                    - gas_out

            # do the integral for metal mass with time
            dmetals = - metals_ast \
            + metals_stars \
            + metals_pre \
            + metals_inf \
            + metals_out

            # do the integral for dust mass with time
            ddust = - md*f.astration(mg,self.sfr(t)) \
              + mdust_stars \
              + mdust_inf \
              - md*mdust_out \
              + mdust_gg \
              - mdust_des

            # to plot dust sources
            dust_source_all = mdust_stars + mdust_gg
            # next time step
            dt = t - prev_t
            prev_t = t
            mg += dmg*dt
            mg_list.append(mg)
            metals += dmetals*dt
            metals_list.append(metals)
            Z = zip(*z_lookup) #metallicity
            md += ddust*dt
            md_all += dust_source_all*dt
            md_gg += mdust_gg*dt
            md_stars += mdust_stars*dt
            dust_list.append(md)
            dust_list_sources.append((md_all, md_stars, md_gg))
            timescales.append((t_des,t_gg))
            if metallicity <= 0.:
                dust_to_metals = 0.
            else:
                dust_to_metals = (md/mg)/metallicity
            dz_ratio_list.append(dust_to_metals)
        print("Gas, metal and dust mass exterior loop %s" % str(datetime.now()-now))
        return time, np.array(mg_list), np.array(metals_list), np.array(Z[1]), \
        np.array(dust_list), np.array(dust_list_sources), \
         np.array(dz_ratio_list), np.array(timescales)

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
                if m > 10.:
                    dm = 0.5
                sn_rate += f.initial_mass_function(m, self.imf_type)*dm
                m += dm
            r_sn = self.sfr(t)*sn_rate # units in N per Gyr
            dt = t - prev_t
            prev_t = t
            sn_rate_list.append(r_sn)
        print("SN rate exterior loop %s" % str(datetime.now()-now))
        return np.array(sn_rate_list)

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
        time = time[time < self.tend]
        for t in time:
            dmstars = self.sfr(t)
            dt = t - prev_t
            prev_t = t
            mstars += dmstars*dt
            mstars_list.append(mstars)
        # Output time and gas mass as Numpy Arrays
        return time, np.array(mstars_list)
