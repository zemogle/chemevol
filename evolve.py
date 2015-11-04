from scipy.integrate import quad
import functions as f
import numpy as np
from lookups import find_nearest
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
            return vals
        except:
            logger.error("No SFH yet")



def ejected_mass_integral(t):
    # integrate the ejected mass function

    return ej
