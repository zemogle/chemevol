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
        # Set up IMF Function
        if (self.imf_type == "Chab" or self.imf_type == "chab" or self.imf_type == "c"):
            self.imf = f.imf_chab
        elif (self.imf_type == "TopChab" or self.imf_type == 'topchab' or self.imf_type == "tc"):
            self.imf = f.imf_topchab
        elif (self.imf_type == "Kroup" or self.imf_type == "kroup" or self.imf_type == "k"):
            self.imf = f.imf_kroup
        elif (self.imf_type == "Salp" or self.imf_type == "salp" or self.imf_type == "s"):
            self.imf = f.imf_salp
        # Declare if destruction on or off
        if (self.destroy == False):
            self.choice_des = 0
        else:
            self.choice_des = 1
        # set up dust source choice 0 = SN, 1 = LIMS, 2 = GG
        if (self.dust_source == "ALL" or self.dust_source == "all" or self.dust_source == "All"):
            self.choice_dust = (1, 1, 1)
        elif (self.dust_source == "SN" or self.dust_source == "Sn" or self.dust_source == "sn"):
            self.choice_dust = (1, 0, 0)
        elif (self.dust_source == "LIMS" or self.dust_source == "Lims" or self.dust_source == "lims"):
            self.choice_dust = (0, 1, 0)
        elif (self.dust_source == "SN+LIMS" or self.dust_source == "sn+lims" or \
              self.dust_source == "LIMS+SN" or self.dust_source == "lims+sn"):
            self.choice_dust = (1, 1, 0)
        elif (self.dust_source == "GG" or self.dust_source == "gg" or self.dust_source == "Gg"):
            self.choice_dust = (0, 0, 1)
        else:
            print ('oops please check the dust sources are in the right format')

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
        dust_sources = []
        timescales = []
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
            z.append([t,metallicity])
            z_lookup = np.array(z)
            sfr_list.append([t,self.sfr(t)])
            sfr_lookup = np.array(sfr_list)

            # GAS
            # astration, inflows, outflows
            gas_ast = self.sfr(t)
            gas_inf = f.inflows(self.sfr(t), self.inflows['xSFR'])
            gas_out = f.outflows(self.sfr(t), self.outflows['xSFR'])

            # METALS
            # astration, inflows, outflows
            metals_ast = metallicity*self.sfr(t)
            metals_inf = self.inflows['metals']*f.inflows(self.sfr(t), self.inflows['xSFR'])
            if self.outflows['metals']:
                metals_out = metallicity*f.outflows(self.sfr(t), self.outflows['xSFR'])
            else:
                metals_out = 0.

            # DUST
            # astration, inflows, outflows, grain growth, destruction
            mdust_ast = md*f.astration(mg,self.sfr(t))
            mdust_inf = self.inflows['dust']*f.inflows(self.sfr(t), self.inflows['xSFR'])

            if self.outflows['dust']:
                mdust_out = (1./mg)*f.outflows(self.sfr(t), self.outflows['xSFR'])
            else:
                mdust_out = 0.

            # destruction timescales + dust mass from grain growth and destruction
            mdust_gg, t_gg = f.graingrowth(self.choice_dust[2],self.epsilon,mg,self.sfr(t),metallicity,md,self.coldfraction)
            mdust_des, t_des = f.destroy_dust(self.choice_des,self.destroy_ism,mg,r_sn,md,self.coldfraction)

            # do the mass integral to get ejected masses for gas, metals, dust
            gas_ej, metals_stars, mdust_stars = \
                                f.mass_integral(self.choice_dust,t, metallicity, sfr_lookup, z_lookup, self.imf)

            # gas mass integral dmg/dt =
            dmg = - gas_ast + gas_ej + gas_inf - gas_out

            # metal mass integral dMz/dt =
            dmetals = - metals_ast + metals_stars + metals_pre + metals_inf - metals_out

            # dust mass integral dMd/dt =
            ddust = - mdust_ast + mdust_stars + mdust_inf - md*mdust_out + mdust_gg - mdust_des

            dust_source_all = mdust_stars + mdust_gg #dust sources stars + grain growth
            dt = t - prev_t             # calculate  next time step
            prev_t = t
            mg += dmg*dt # gas mass integral
            metals += dmetals*dt # metal mass integral
            md += ddust*dt # dust mass integral
            md_all += dust_source_all*dt # dust mass sources integral
            md_gg += mdust_gg*dt # dust source from grain growth only
            md_stars += mdust_stars*dt # dust source from stars only
            Z = zip(*z_lookup) # write metallicity to an array
            s_f_r = zip(*sfr_lookup) # write SFR lookup array
            dust_sources.append((md_all, md_stars, md_gg)) # write array of dust sources
            timescales.append((t_des,t_gg)) # write array for grain growth & destruction timescales
            if metallicity <= 0.:  # write dust/metals ratio but == 0 when metals  = 0
                dust_to_metals = 0.
            else:
                dust_to_metals = (md/mg)/metallicity
            all_results.append((t,mg,metals,metallicity,md,dust_to_metals,self.sfr(t)*1e-9))
        print("Gas, metal and dust mass exterior loop %s" % str(datetime.now()-now))
        return  np.array(dust_sources), np.array(timescales), np.array(all_results)

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
