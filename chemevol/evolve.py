'''
Chemevol - Python package to read in a star formation history file,
input galaxy parameters and run a chemical evolution model to determine the evolution
of gas, metals and dust in galaxies.

Running this script will produce a results data file (filename.dat) with file name given by user

The code is based on Morgan & Edmunds 2003 (MNRAS, 343, 427)
and described in detail in De Vis et al 2017, 2021 (MNRAS, 471, 1743; ).

If you use this code, please do cite the above papers.

Copyright (C) 2019 Haley Gomez, Edward Gomez, Pieter De Vis and Simon Schofield, 
Cardiff University and LCOGT. The code has been contributed by Kate Rowlands.

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

********************************************************************************
'''

from chemevol.functions import extra_sfh_and_inflows, astration, imf_chab, imf_topchab, \
    imf_salp, imf_kroup, initial_mass_function, initial_mass_function_integral, \
    lookup_fn, lookup_taum, mass_integral, t_lifetime, graingrowth, destroy_dust, \
    gas_inandout, metals_inandout, dust_inandout, outflows_feldmann, recycle,\
    graingrowth_THEMIS_cloud, graingrowth_THEMIS_diffuse, destroy_dust_SN_THEMIS, destroy_dust_frag_THEMIS

from astropy.table import Table
import numpy as np
import random as rn
from chemevol.lookups import find_nearest, lookup_fn, t_lifetime, lookup_taum
import logging
from datetime import datetime
import os.path
import json
from astropy.cosmology import Planck13, z_at_value
import astropy.units as u

# Set up a logging system so that we can store warnings without them completely filling up the terminal screen and drowning out other usefull output.
logging.captureWarnings(True)
logging.basicConfig(filename='Warnings.log',level=logging.WARNING)
logger = logging.getLogger('chem')

class ChemModel:
    def __init__(self, **inputs):
        '''
        set initial parameters from input dictionary and sort out choices
        for IMF, dust source
        '''
        try:
            self.name = inputs['name']
            self.gasmass_init = inputs['gasmass_init']
            self.starmass_init = inputs['starmass_init']
            self.dustmass_init = inputs['dustmass_init']
            self.Z_init = inputs['Z_init']
            self.gamma = inputs['gamma']
            self.tend = inputs['t_end']
            self.tstart = inputs['t_start']
            self.imf_type = inputs['IMF_fn']
            self.isotopes = inputs['isotopes']
            self.SNyield = inputs['SNyield']
            self.AGByield = inputs['AGByield']
            self.totyields = inputs['totyields']
            self.SFH_file = inputs['SFH']
            self.add_bursts = inputs['add_bursts']
            self.inflow_metalfractions = np.array(inputs['Pristine_isotope_fractions'])
            self.inflows = inputs['inflows']
            self.outflows = inputs['outflows']
            self.outflows['model']='nelson'
            self.recycle = inputs['recycle']
            self.dust_source = inputs['dust_source']
            self.THEMIS = inputs['use_THEMIS']
            self.delta_lims = inputs['delta_lims_fresh']
            self.reduce_sn = inputs['reduce_sn_dust']
            self.destroy = inputs['destroy']
            self.frag = inputs['fragmentgrains']
            self.eff_snrate = inputs['effective_snrate_factor']
            self.graingrowth = inputs['graingrowth']
            self.graingrowth2 = inputs['graingrowth2']
            self.availablefraction = inputs['available_metal_fraction']
            self.coldfraction = inputs['cold_gas_fraction']
            self.sfh_file = self.SFH_file
            self.load_sfh_and_inflow()
        except KeyError:
            print('You must provide initial parameters in the correct format for model %s'%(self.name))
            exit()
        # Set up IMF Function determined by user, allow for variety of spellings
        if (self.imf_type in ["Chab", "chab", "c"]):
            self.imf = imf_chab
        elif (self.imf_type in ["TopChab", "topchab","tc"]):
            self.imf = imf_topchab
        elif (self.imf_type in ["Kroup", "kroup", "k"]):
            self.imf = imf_kroup
        elif (self.imf_type in ["Salp", "salp", "s"]):
            self.imf = imf_salp
        # set it up so that the reduction factor is 1 (i.e. no change) when the reduction is turned off.     
        if not self.reduce_sn['on'] or self.reduce_sn['factor'] == 0:
            self.reduce_sn = 1
        else:
            self.reduce_sn = self.reduce_sn['factor']
        # set up dust source choice from user: sn = True; SN dust on, lims = True; LIMS dust on, gg = Grain Growth
        if self.dust_source in ["ALL", "all", "All"]:
            self.choice_dust = {
                                    'sn' : True,
                                    'lims' : True,
                                    'gg' :True
                                }
        elif self.dust_source in ["SN", "Sn", "sn"]:
            self.choice_dust =  {
                                    'sn' : True,
                                    'lims' : False,
                                    'gg' :False
                                }
        elif self.dust_source in ["LIMS", "Lims", "lims"]:
            self.choice_dust =  {
                                    'sn' : False,
                                    'lims' : True,
                                    'gg' :False
                                }
        elif self.dust_source in ["SN+LIMS", "sn+lims", "LIMS+SN", "lims+sn"]:
            self.choice_dust =  {
                                    'sn' : True,
                                    'lims' : True,
                                    'gg' :False
                                }
        else:
            print('oops please check the dust sources are in the right format and try again for model %s'%(self.name))
            exit()

    def load_sfh_and_inflow(self):
        '''
        takes in input SFH and inflow file and extend backwards to start from 1e-3 Gyr
        '''
        try:
            vals = np.loadtxt(self.SFH_file)
            scale = [1e-9,1e9,1e9] # Gyr conversions for time, SFR (because we want to do dt integral over Gyrs)
            sfh = vals*scale # converts time in Gyr and SFR in Msun/Gyr
            
            # extrapolates SFH back to 0.001Gyr using SFH file and power law (gamma)
            final_sfh, final_inflows = extra_sfh_and_inflows(sfh, self.gamma, self.tstart)
            self.sfh = np.array(final_sfh)
            self.inflowvals= np.array(final_inflows)
            dts=np.diff(self.sfh[:,0])
            dts=np.append(dts[0],dts)
            if self.inflows['mass']>0:
                self.inflowvals[:,1] = self.inflowvals[:,1]/np.sum((self.inflowvals[:,1]*dts)[np.where(self.inflowvals[:,0]<self.tend)])*self.inflows['mass'] # if available, rescale inflow mass to provided value
            else:
                print("Since no inflows are specified, half of the initial gas mass is used as inflows.")    
                self.inflowvals[:,1] = self.inflowvals[:,1]/np.sum((self.inflowvals[:,1]*dts)[np.where(self.inflowvals[:,0]<self.tend)])*self.gasmass_init/2. #if not provided, set inflow mass the same as initial gas mass
                self.gasmass_init = self.gasmass_init/2.
            return vals[1]
        except Exception as e:
            print("File '%s' will not parse %s for model %s" % (self.SFH_file, e, self.name))
            self.sfh = None

    def sfr(self, t):
        '''
        We define sfr as a function to lookup the apporpriate SFR/SFE value at a given time.
        '''
        try:
            # look up nearest SFR/SFE value at any specified time
            vals = find_nearest(self.sfh,t)
            return vals[1]   
        except:
            print('No SFH yet for model %s'%(self.name))

    def inflow(self, t):
        '''
        define inflow as function to look up nearest inflow rate value at any specified time
        '''
        try:
            vals = find_nearest(self.inflowvals,t)
            return vals[1] 

        except:
            print('No inflows yet for model %s'%(self.name))      
                
    def z_at_time(t,times,redshifts):
         nt = len(times)
         it = np.maximum(1,np.minimum(nt-1,bisect.bisect(-times,-t)))
         fnext = np.maximum(0.,np.minimum(1.,(t - times[it-1]) / (times[it] - times[it-1])))
         fprev = 1.-fnext
         redshift =  redshifts[it]*fnext + redshifts[it-1]*fprev
         return redshift

    def gas_metal_dust_mass(self, sn_rate):
            '''
            Calculates the gas, metal and dust mass from the different sources and sinks in the model.
            '''

            # initialize
            mg = self.gasmass_init
            mstars = self.starmass_init
            md_tot = self.dustmass_init
            f_c=self.coldfraction
            md_diff = self.dustmass_init*(1-f_c)
            md_cloud = self.dustmass_init*f_c
            md_all = 0
            md_stars = 0
            md_gg = 0
            md_gg_diff = 0
            md_gg_cloud = 0
            md_cloud = 0
            md_IGM = 0

            # initialise metals as arrays with the same size as self.isotopes 
            nisotopes=len(self.isotopes)
            metals = self.Z_init*mg*self.inflow_metalfractions
            metals_IGM = np.zeros(nisotopes)

            prev_t = self.tstart 
            z_lookup = np.zeros(1+nisotopes)
            sfr_list = []
            sfr_lookup = []
            all_results = []
            
            time = self.sfh[:,0] # sfr is in units of Msun Gyr^-1
            dts = np.diff(time) 
            dts = np.append(dts[0],dts)
            now = datetime.now()
            remainingburststeps=0
            
            gas_inf=0
            gas_rec=0 
            gas_out=0
            mg_IGM=0

            # arrays to store the outflows of gas, dust and metals, which will be used to determine the recycling rates
            recycle_gas = np.zeros(len(time))
            recycle_dust = np.zeros(len(time))
            recycle_Z = np.zeros((len(time),nisotopes))

            redshift_lookup = np.logspace(-3,2.5,100)
            time_lookup = Planck12.age(redshift_lookup).value
            # TIME integral
            for item, t in enumerate(time[(time < self.tend) & (time>self.tstart)]):
                dt = t - prev_t             # calculate  next time step
                try:
                    redshift=z_at_time(t,time_lookup,redshift_lookup)
                except:
                    redshift=0.
                        
                fg=mg/(mg+mstars)
                metallicity = metals[0]/mg
                metallicities = metals/mg
                
                # SET UP SFR
            
                # If the SFH_file ends in .sfe it is given as a star formation efficiency, and thus needs to be multiplied by gas mass.
                # In that case, the SFE prescription described in De Vis et al (2020) Section 3.2 will be used.
                if ".sfe" in self.SFH_file:
                    SFR = self.sfr(t)*(mg)*(max(1e5,mstars)/1e9)**0.25*(1+np.exp(mstars/(10*mg)))**-3*(1+redshift)**-1

                    # If the add_burst parameter is True, add bursts so that there is a 50% chance that there is a burst in any 2 Gyr.
                    if self.add_bursts:
                        if remainingburststeps==0:    # check is burst if currently ongoing, if not check if a new burst starts. 
                            Pburst=1.-0.5**(dt/2.)  # 50% chance in 2 Gyr
                            if rn.random()<Pburst:
                                remainingburststeps=rn.randint(1, 10)
                                facburst=10**(rn.uniform(np.log10(0.004),np.log10(0.1))) # during this burst, between 0.3 and 10% of the total stellar mass is formed during this burst.
                                SFRburst=mstars*facburst/dt/remainingburststeps
                        if remainingburststeps>0:
                            SFR=SFR+SFRburst
                            remainingburststeps+=-1

                # predetermined SFH can also be tested using SFH files with any other extension.             
                else:
                    SFR = self.sfr(t)      

                # SN rate
                r_sn = SFR*sn_rate [item]
                

                # SET UP INFLOWS AND OUTFLOWS

                # How much gas is lost or gained dure to outflows/inflows
                gas_inf,gas_out,outflows = gas_inandout(\
                    redshift,\
                    self.inflows['on'],\
                    self.outflows['on'],\
                    self.outflows['model'],\
                    self.inflow(t),\
                    SFR,\
                    mstars,\
                    self.outflows['reduce'])

                # Reduce outflows so that no more gas can be removed than is currently present
                # We allow a maximum of 50% of the gas to blown out in a single 30 Myr timestep. 
                if gas_out*dt>0.5*mg:
                    logger.warning('SFR has been reduced by factor %s at time= %s Gyr because more than 50%% of the gas was being blown out in a single 30 Myr timestep for model %s'%((gas_out*dt)/(0.5*mg),t,self.name))
                    SFR=SFR*0.5*mg/(gas_out*dt)
                    gas_inf,gas_out,outflows = gas_inandout(\
                        redshift,\
                        self.inflows['on'],\
                        self.outflows['on'],\
                        self.outflows['model'],\
                        self.inflow(t),\
                        SFR,\
                        mstars,\
                        self.outflows['reduce'])
                    r_sn = SFR*sn_rate [item]

                # How much metals are lost or gained due to outflows/inflows
                metals_inf,metals_out = metals_inandout(
                    self.inflows['metals'],\
                    self.outflows['metals'],\
                    metallicities,\
                    self.inflow_metalfractions,\
                    gas_inf,\
                    gas_out,\
                    nisotopes)

                # How much dust is lost or gained due to outflows/inflows
                mdust_inf,mdust_out = dust_inandout(
                    self.inflows['dust'],\
                    self.outflows['dust'],\
                    (md_tot/mg),\
                    gas_inf,\
                    gas_out)

                # Append arrays for later needs
                z_lookup=np.vstack((z_lookup,np.append(t,metallicities)))
                sfr_list.append([t,SFR])
                sfr_lookup = np.array(sfr_list)

               
                # SET UP COMPONENTS OF THE INTEGRALS OF STARS, GAS, METALS AND DUST 
                
                # First, determine the gas, metals and dust comsumed in astration

                gas_ast = SFR # gas lost due to astration
                metals_ast = astration(metals,mg,SFR)
                mdust_ast = astration(md_tot,mg,SFR)
            
                # Second, determine the grain growth and destruction rates according to the THEMIS model or prescriptions from De Vis et al (2017b).
                           
                if self.THEMIS:
                    mdust_des_SN, t_des = destroy_dust_SN_THEMIS(self.destroy['on'], self.destroy['mass'], mg, r_sn, \
                        md_tot, self.coldfraction, self.eff_snrate)
                    mdust_des_frag, t_frag = destroy_dust_frag_THEMIS(self.frag['on'], self.frag['tau'], mg, SFR/mstars, \
                        md_tot, self.coldfraction)

                    mdust_gg_diff, t_gg_diff =  graingrowth_THEMIS_diffuse(self.choice_dust['gg'],self.graingrowth2,mg, SFR, metallicity, md_tot, self.coldfraction, self.availablefraction)
                    if mdust_gg_diff<0:
                        mdust_gg_diff=0
                    mdust_gg_cloud, t_gg_cloud = graingrowth_THEMIS_cloud(self.choice_dust['gg'],self.graingrowth,mg, SFR, \
                        metallicity, md_tot, self.coldfraction, min(self.availablefraction*2.45,1.))

                    mdust_gg=mdust_gg_diff + mdust_gg_cloud
                    
                else:
                    mdust_gg, t_gg = graingrowth(self.choice_dust['gg'],self.graingrowth,mg, SFR, \
                        metallicity, md_tot, self.coldfraction, self.availablefraction)
                    mdust_des, t_des = destroy_dust(self.destroy['on'], 135*self.destroy['mass'], mg, r_sn, \
                        md_tot, self.coldfraction , self.eff_snrate)

                    t_gg_diff=t_gg
                    t_gg_cloud=t_gg
                    t_frag=0

                '''
                Next, get ejected masses from stars when they die
                gas_ej = e(t): ejected gas mass from stars of mass m at t = taum
                metals_stars = ez(t): ejected metal mass from stars of mass m at t = taum (fresh + recycled)
                mdust_stars = ed(t): ejected dust mass from stars of mass m at t = taum (fresh + recycled)
                '''
                gas_ej, metals_stars, mdust_stars = \
                        mass_integral(self.choice_dust, self.delta_lims, self.reduce_sn, t, metallicity, sfr_lookup, z_lookup, self.imf, self.SNyield, self.AGByield, self.totyields, nisotopes, self.isotopes)


                # Calculate when the current outflows will be recycled and store this info in recycle_gas,recycle_dust,recycle_Z for later use
                if self.recycle['on']:
                    recycle_gas,recycle_dust,recycle_Z = \
                        recycle(t, dt, redshift, mstars, time, dts, recycle_gas, recycle_dust, recycle_Z, outflows, gas_out, mdust_out, np.asarray(metals_out), self.recycle['esc_prob_perGyr'], self.recycle['reaccr_time_factor'], nisotopes)
                
                # Read the current recycling rate of gas, dust and metals (using the current timestep given by [item])
                gas_rec, mdust_rec, metals_rec = recycle_gas[item],recycle_dust[item],recycle_Z[item]

                '''
                DO THE INTEGRALS
                integrate over time for gas, metals, dust and stars
                all time units should be in Gyr or per Gyr
                '''

                #STARS
                dmstars = SFR - gas_ej
                mstars += dmstars*dt

                #GAS                
                dmg = -gas_ast + gas_ej + gas_inf - gas_out + gas_rec/dt
                mg += dmg*dt
                
                #METALS
                dmetals = -metals_ast + metals_stars + metals_inf - metals_out + metals_rec/dt
                metals = metals+dmetals*dt

                #DUST
                if self.THEMIS:
                    if mdust_des_frag*dt>md_diff:
                        logger.warning('Rescaling fragmentation rate at time= %s Gyr to avoid destroying more dust than is present in the diffuse ISM for model %s'%(t,self.name))
                        mdust_des_frag=md_diff/dt
                    ddust_diff = mdust_inf - mdust_out*(1-f_c) + mdust_rec/dt + mdust_gg_diff - mdust_des_SN - mdust_des_frag
                    ddust_cloud = -mdust_ast + mdust_stars - mdust_out*f_c + mdust_gg_cloud 

                    md_diff += ddust_diff*dt 
                    md_cloud += ddust_cloud*dt 

                    if(md_diff<0.0) : md_diff = 0.0
                    if(md_cloud<0.0) : md_cloud = 0.0
                    #reduce diffuse dust content if the dust-to-metal ratio is unrealistically high
                    if md_diff/((1-f_c)*metals[0])>self.availablefraction:
                        md_diff=self.availablefraction*(1-f_c)*metals[0]

                    #reduce cloud dust content if the dust-to-metal ratio is unrealistically high
                    if md_cloud/(f_c*metals[0])>self.availablefraction*2.45:
                            correction_dust_depletion=self.availablefraction*metals[0]*f_c*2.45/md_cloud
                            md_cloud=self.availablefraction*metals[0]*f_c*2.45
                    else:
                            correction_dust_depletion=1.        
                   
                    md_tot= md_diff+md_cloud  

                else:
                    ddust = -mdust_ast + mdust_stars + mdust_inf - mdust_out + mdust_rec/dt + mdust_gg - mdust_des
                    md_tot += ddust*dt

                # Dust_source_all separates out the dust sources (Md vs t) wihtout including sinks (Astration etc)
                # and grain growth separately (this is the Md vs time contributed by dust sources)
                dust_source_all = mdust_stars + mdust_gg
                md_all += dust_source_all*dt # dust mass sources integral
                md_gg += mdust_gg*dt # dust source from grain growth only
                md_gg_diff += mdust_gg_diff*dt # dust source from grain growth only
                md_gg_cloud += mdust_gg_cloud*dt # dust source from grain growth only
                md_stars += mdust_stars*dt # dust source from stars only

                # Sanity checks
                if mg/(mstars+mg)<0.0001:
                    # exit program if gas fraction goes below 1%
                    print('Oops you have no interstellar medium left for model %s, model is terminated'%(self.name))
                    break
                
                if metals[0] < 0:
                    # exit program if all ISM removed
                    print('Oops you have no metals left for model %s, model is terminated'%(self.name))
                    break

                if md_tot-0.9*mdust_gg_cloud*dt*correction_dust_depletion<0:   #corrected for cloud dissociation
                    print('Oops you have no dust left for model %s, model is terminated'%(self.name))
                    break

                #INFLOWS AND OUTFLOWS in given timestep to store in output
                mg_outflows = gas_out*dt
                mg_recycled = gas_rec
                mg_inf = gas_inf*dt

                #track mas of gas, metals and dust in the InterGalactic Medium
                mg_IGM += mg_outflows - mg_recycled
                md_IGM += mdust_out*dt - mdust_rec
                metals_IGM += metals_out*dt - metals_rec                

                # determine dust/metals ratio 
                if mg <= 0. or metals[0] <=0:  #This can only be true if something went wrong
                    dust_to_metals = 0.
                else:
                    dust_to_metals = md_tot/metals[0]
                
                #STORE RESULTS
                all_results.append(np.hstack(([t, redshift, mg, mstars, metallicity, \
                                    md_tot, dust_to_metals, SFR*1e-9, \
                                    md_all, md_stars, md_gg, md_gg_diff, md_gg_cloud, t_des, t_frag, t_gg_diff, t_gg_cloud, \
                                    mg_outflows, mg_recycled, mg_inf, mg_IGM, md_IGM, md_diff, md_cloud], 
                                    metals, metals_IGM)))

                #Clouds are dissociated
                md_tot+=-0.9*mdust_gg_cloud*dt*correction_dust_depletion  #90% of the newly accreted dust in clouds evaporates as the clouds are dissociated (ices etc evaporate in the harsher environment)
  
                #Reset cloud and diffuse dust mass as new clouds are formed with the same dust-to-gas ratio as the surrounding diffuse ISM
                md_cloud=md_tot*f_c
                md_diff=md_tot*(1-f_c)
                prev_t=t

            #print("Gas, metal and dust mass exterior loop done after %s for model %s" % (str(datetime.now()-now),self.name))
            return np.array(all_results)

    def supernova_rate(self):
        '''
        Calculates the SN rate at time t by integrating over mass m
        '''
        # initialize
        sn_rate_list = []
        dm = 0.01
        prev_t = 1e-3
        # define time array
        time = self.sfh[:,0] # this is in units of Gyrs
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
                sn_rate += initial_mass_function(m, self.imf_type)*dm
                m += dm
            r_sn = sn_rate # this is in units of Msun Gyr^-1 x Msun --> Gyr^-1
            dt = t - prev_t
            prev_t = t
            sn_rate_list.append(r_sn) # roughly is ~10 per century at early times and <1 per century at late time
        return np.array(sn_rate_list)


'''
Class for running multiple models. 
Parallel processing option is provided in separate file.
'''

class BulkEvolve: 
    def __init__(self, filename):
        if os.path.isfile(filename):
            self.filename = filename
        else:
            print('File {} does not exist'.format(filename))
        return

    def upload_json(self):
        try:
            with open(self.filename) as data_file:
                data = json.load(data_file)
            self.inits = data
        except ValueError:
            print('Cannot read: Are you sure this is a JSON file?')
        return


    def upload_csv(self):
        names = ['name', 'gasmass_init', 'starmass_init', 'dustmass_init', 'Z_init', 'SFH', 'add_bursts', 't_end', 't_start',\
         'gamma', 'IMF_fn', 'dust_source', 'cold_gas_fraction', 'use_THEMIS', 'delta_lims_fresh', 'reduce_sn_dust_on', 'reduce_sn_dust_factor',\
         'destroy_on', 'mass_destroy', 'fragment_on', 'fragment_tau', 'effective_snrate_factor', 'graingrowth', 'graingrowth2',\
         'inflows_on', 'inflows_mass', 'inflows_metals', 'inflows_xSFR', 'inflows_dust', 'outflows_on','outflows_metals', 'outflows_dust', 'outflows_reduce',\
         'recycle_on', 'esc_prob_perGyr', 'reaccr_time_factor', 'available_metal_fraction', 'SNyield', 'AGByield','totyields','isotopes','Pristine_isotope_fractions']

        alttype = np.dtype([('f0','S30'), ('f1', '<f8'), ('f2', '<f8'), ('f3','<f8'), ('f4','<f8'), ('f5','S30'), ('f6','bool'),('f7','<f8'), ('f8','<f8'), \
                    ('f9','<f8'), ('f10','S10'), ('f11','S20'), ('f12','<f8'), ('f13','bool'),('f14','<f8'),('f15','bool'),('f16','<f8'),\
                    ('f17','bool'),('f18','<f8'), ('f19','bool'), ('f20','<f8'), ('f21','<f8'), ('f22','<f8'), ('f23','<f8'),\
                     ('f24','bool'), ('f25','<f8'), ('f26','<f8'), ('f27','<f8'),('f28','<f8'),('f29','bool'), ('f30','bool'), ('f31','bool'), ('f32','<f8'),\
                     ('f33','bool'),('f34','<f8'),('f35','<f8'),('f36','<f8'),('f37','S30'),('f38','S30'),('f39','bool'),('f40','S30'),('f41','S30')])
        try:
            data = np.genfromtxt("modelGrid_dust.csv", dtype=alttype,delimiter=',', autostrip=True, names=names)
        except ValueError:
            print('Cannot read: Are you sure this is a CSV file?')
            stop
        init_list = []

        for i in range(0,len(data)):
            gal_tup = zip(names, data[i])
            gal_data = dict(gal_tup)
            gal_data['reduce_sn_dust'] = {'on': gal_data['reduce_sn_dust_on'],
                                          'factor': gal_data['reduce_sn_dust_factor']}
            gal_data['inflows'] = { 'on': gal_data['inflows_on'],
                                    'mass': gal_data['inflows_mass'],
                                    'metals': gal_data['inflows_metals'],
                                    'xSFR': gal_data['inflows_xSFR'],
                                    'dust': gal_data['inflows_dust']}
            gal_data['outflows'] ={ 'on': gal_data['outflows_on'],
                                    'metals': gal_data['outflows_metals'],
                                    'dust': gal_data['outflows_dust'],
                                    'reduce': gal_data['outflows_reduce']}
            gal_data['recycle'] = { 'on': gal_data['recycle_on'],
                                    'esc_prob_perGyr': gal_data['esc_prob_perGyr'],
                                    'reaccr_time_factor': gal_data['reaccr_time_factor']}                        
            gal_data['destroy'] = { 'on': gal_data['destroy_on'],
                                    'mass': gal_data['mass_destroy']}
            gal_data['fragmentgrains'] = {'on': gal_data['fragment_on'],
                                          'tau': gal_data['fragment_tau']}                                                            
            gal_data['isotopes'] = gal_data['isotopes'].split(";")
            gal_data['Pristine_isotope_fractions' ] = [float(stri) for stri in gal_data['Pristine_isotope_fractions'].split(";")]                           
            
            init_list.append(gal_data)
        self.inits = init_list
        return


    def evolve_all(self):
        '''
        call modules to run the model:
        snrate:         SN rate at each time step - this also sets time array
                        so ch.supernova_rate() must be called first to set
                        time array for the entire code

        all results:    contains all relevant chemical evolution parameters to be outputted
        '''
        snrate = []
        all_results = []
        galaxies = []

        for item in self.inits:
            print('Starting run on {}'.format(item['name']))
            ch = ChemModel(**item)

            snrate = ch.supernova_rate()
            all_results = ch.gas_metal_dust_mass(snrate)
            
            # write all the results to a dictionary
            params = {'time' : all_results[:,0],
            'z' : all_results[:,1],
            'mgas' : all_results[:,2],
            'mstars' : all_results[:,3],
            'metallicity' : all_results[:,4],
            'mdust' : all_results[:,5],
            'dust_metals_ratio' : all_results[:,6],
            'sfr' : all_results[:,7],
            'dust_all' : all_results[:,8],
            'dust_stars' : all_results[:,9],
            'dust_ism' : all_results[:,10],
            'dust_diff' : all_results[:,11],
            'dust_cloud' : all_results[:,12],
            'time_destroy' : all_results[:,13],
            'time_fragment' : all_results[:,14],
            'time_gg_diffuse' : all_results[:,15],
            'time_gg_cloud' : all_results[:,16],
            'mgas_outflow' : all_results[:,17],
            'mgas_recycled' : all_results[:,18],
            'mgas_inflow' : all_results[:,19],
            'mgas_IGM' : all_results[:,20],
            'mdust_IGM' : all_results[:,21],
            'mdust_diffuse' : all_results[:,22],
            'mdust_cloud' : all_results[:,23]}   

            #compute additional parameters           'dust_ism' : all_results[:,9],
            params['fg'] = params['mgas']/(params['mgas']+params['mstars'])
            params['ssfr'] = params['sfr']/params['mstars']
            
            paramsorder=['time','z','fg','mgas','mstars','mdust','mdust_diffuse',\
                        'mdust_cloud','metallicity','dust_metals_ratio','sfr',\
                        'ssfr','dust_all','dust_stars','dust_ism','dust_diff','dust_cloud','time_destroy',\
                        'time_fragment','time_gg_diffuse','time_gg_cloud','mgas_outflow',\
                        'mgas_recycled','mgas_inflow','mgas_IGM','mdust_IGM']

            #properties for the isotypes specified in input
            for iso,nameiso in enumerate(item['isotopes']):
                Miso=all_results[:,24+iso]
                params['M'+nameiso]=Miso
                Miso_IGM=all_results[:,24+iso+len(item['isotopes'])]
                params['M'+nameiso+'_IGM']=Miso_IGM
                paramsorder+=['M'+nameiso,'M'+nameiso+'_IGM']

            # write out to file based on 'name' identifier
            name = item['name']
            
            ascii.write(params,str(name+'.csv'),format='csv',names=paramsorder)
            
            # if you want an array including every inits entry:
            galaxies.append(params)
        self.results = galaxies
        return

