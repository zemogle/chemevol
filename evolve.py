from scipy.integrate import quad
import functions as f
import numpy as np
from lookups import find_nearest, lifetime_lookup, t_lifetime
import logging

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
            if not self.SFH_file:
                self.SFH_file = 'Milkyway.sfh'
            self.sfh_file = self.SFH_file
            self.load_sfh()
        except KeyError:
            logger.error('You must provide initial parameters')

    def load_sfh(self):
        try:
            vals = np.loadtxt(self.SFH_file)
            self.sfh = vals
        except:
            logger.error("File '%s' will not parse" % self.SFH_file)
            self.sfh = None

    def sfr(self, t):
        try:
            vals = find_nearest(self.sfh,t)
            return vals[1]
        except:
            logger.error("No SFH yet")

    def ejected_mass(self, t):
        mu = t_lifetime[-1][0]
        m = lifetime_lookup(t_lifetime,'lifetime_low_metals',t)[0]
        dm = 0.5
        em = 0.
        while m <= mu:
            m += dm
            em += f.ejected_gas_mass(m, self.sfr(t), self.imf_type) * dm
        print(em)
        return em

    def gas_mass(self):
        t = self.sfh[0][0] + 0.00001
        t_end = self.sfh[-1][0]
        dlogt=(np.log10(t_end)-np.log10(t))/100.
        mg = self.gasmass_init
        while t <= t_end:
            dmg = - self.sfr(t) + self.ejected_mass(t) + f.inflows(self.sfr(t), self.inflows).value + f.outflows(self.sfr(t), self.outflows).value
            mg += dmg * 10.**dlogt
            t += 10.**(np.log10(t)+dlogt)
            print('dMg', dmg, self.sfr(t), t)
        return mg