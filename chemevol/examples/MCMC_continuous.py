from chemevol import ChemModel
import numpy as np
import emcee
from multiprocessing import Pool
import re
from copy import deepcopy
from scipy.interpolate import interp1d
from contextlib import contextmanager

@contextmanager
def terminating(thing):
    try:
        yield thing
    finally:
        thing.terminate()

# Grid parameters (these are only used for setting the limits on the priors)
Minigrid=[1e7,1e8,1e9,1e10,1e11,3e11]
SFHgrid=["average.sfe","fast.sfe","slow.sfe","average_bursts.sfe","fast_bursts.sfe","slow_bursts.sfe"]
SNgrid=[1,5,20,80]
destroygrid=[0,15,30]
fragmentgrid=[0.005,0.05,0.5,1,5]
ggcloudgrid=[1000,2000,4000,8000,16000]
ggdiffusegrid=[0,5,10]
favailablegrid=[0.1,0.2,0.3,0.4]

# Marginalisation values for initial mass and time
Minigrid_int=10**(np.arange(7,11.5,0.25))
Tendgrid=np.arange(0,13.8,0.3)

# Read in the data files
dat_all=np.load("data_dust_LTG.npy")
err_all=np.load("data_err_dust_LTG.npy")

# Read in model grid and then take the first line to use as a template for the model input dictionary
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

data = np.genfromtxt('./modelGrid_dust2.csv', dtype=alttype,delimiter=',', autostrip=True, names=names)

#take the first line of the grid and make the dictionary
gal_tup = zip(names, data[0])
gal_data = dict(gal_tup)
gal_data['reduce_sn_dust'] = {'on': gal_data['reduce_sn_dust_on'],
                                'factor': gal_data['reduce_sn_dust_factor']}
gal_data['inflows'] = { 'on': gal_data['inflows_on'],
                                'mass': gal_data['inflows_mass'],
                                'metals': gal_data['inflows_metals'],
                                'xSFR': gal_data['inflows_xSFR'],
                                'dust': gal_data['inflows_dust']}
gal_data['outflows'] = {'on': gal_data['outflows_on'],
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

def runModel(item):
    ch = ChemModel(**item)
    '''
    call modules to run the model:
    snrate:         SN rate at each time step - this also sets time array
                    so ch.supernova_rate() must be called first to set
                    time array for the entire code

    all results:    t, mg, m*, mz, Z, md, md/mz, sfr,
                    dust_source(all), dust_source(stars),
                    dust_source(ism), destruction_time, graingrowth_time,
                    oxygen_mass
    '''
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

    #compute additional parameters
    params['fg'] = params['mgas']/(params['mgas']+params['mstars'])
    params['ssfr'] = params['sfr']/params['mstars']

    paramsorder=['time','z','fg','mgas','mstars','mdust','mdust_diffuse',\
                'mdust_cloud','metallicity','dust_metals_ratio','sfr',\
                'ssfr','dust_all','dust_stars','dust_ism','time_destroy',\
                'time_fragment','time_gg_diffuse','time_gg_cloud','mgas_outflow',\
                'mgas_recycled','mgas_inflow','mgas_IGM','mdust_IGM']

    #properties for the isotypes specified in input
    for iso,nameiso in enumerate(item['isotopes']):
        Miso=all_results[:,24+iso]
        params['M'+nameiso]=Miso
        Miso_IGM=all_results[:,24+iso+len(item['isotopes'])]
        params['M'+nameiso+'_IGM']=Miso_IGM
        paramsorder+=['M'+nameiso,'M'+nameiso+'_IGM']


    return params

def currentModel(Mini,SFH,theta):
    # Change the relevant parameters in the input dictionary and then run the model for these parameters
    curModel=np.zeros((6,len(Tendgrid)))
    gal_data2=deepcopy(gal_data)
    gal_data2['gasmass_init']=Mini
    if "bursts" in SFH:
        gal_data2['SFH']=re.sub('_bursts', '', SFH)
        gal_data2['bursts']= True
    else:
        gal_data2['SFH']=SFH
        gal_data2['bursts']= False
    gal_data2['reduce_sn_dust_factor']=10**theta[0]
    gal_data2['mass_destroy']=theta[1]
    gal_data2['fragment_tau']=10**theta[2]
    gal_data2['epsilon_grain']=10**theta[3]
    gal_data2['epsilon_grain2']=theta[4]
    gal_data2['available_metal_fraction']=theta[5]
    gal_data2['reduce_sn_dust'] = {'on': True, 'factor': gal_data2['reduce_sn_dust_factor']}
    gal_data2['destroy'] = { 'on': True,'mass': gal_data2['mass_destroy']}
    gal_data2['fragmentgrains'] = {'on': True,'tau': gal_data2['fragment_tau']}
    Model=runModel(gal_data2)
    for it in range(len(Tendgrid)):
        itime=(np.abs(Model['time']-Tendgrid[it])).argmin()
        curModel[0,it] = np.log10(Model['mgas'][itime])
        curModel[1,it] = np.log10(Model['mstars'][itime])
        curModel[2,it] = np.log10(Model['sfr'][itime])
        curModel[3,it] = 12+np.log10(((Model['MO']-0.238*Model['mdust'])/16)/(Model['mgas']/1.36))[itime]
        curModel[4,it] = -np.log10(((Model['MO']-0.238*Model['mdust'])/16)/(Model['MN']/14.))[itime]
        curModel[5,it] = np.log10(Model['mdust'][itime])
    del gal_data2,Model
    return curModel


def findchimarg_all(theta):
    Ptots=np.zeros(len(dat_all[0]))
    for SFH_i in range(6):
        Models_i=np.empty((len(Minigrid),6,len(Tendgrid)))
        for ig in range(len(Minigrid)): # run models for the 6 initial masses, and interpolate between them to get the other initial masses
            Models_i[ig]=currentModel(Minigrid[ig],SFHgrid[SFH_i],theta)
        intf=interp1d(Minigrid,Models_i,axis=0)
        for Mini_i in range(len(Minigrid_int)):
            mod=intf(Minigrid_int[Mini_i])
            for t in range(len(mod[0])):
                for it in range(len(dat_all[0])):
                    chisum=0    
                    for i in range(len(dat_all)):
                        ob2=mod[i][t]
                        if dat_all[i][it]>-98:
                            chi=(ob2-dat_all[i][it])**2/err_all[i][it]**2
                        if np.isfinite(chi):
                            chisum+=chi
                    Ptots[it]+=np.exp(-0.5*chisum)
    return -2*np.log(Ptots)        

def lnlike(theta):
    return -0.5*(np.sum(findchimarg_all(theta)))

def lnprior(theta):
    if np.log10(min(SNgrid)) <= theta[0] <= np.log10(max(SNgrid)) and min(destroygrid) <=theta[1] <= max(destroygrid) and np.log10(min(fragmentgrid)) <=theta[2] <= np.log10(max(fragmentgrid)) and np.log10(min(ggcloudgrid)) <=theta[3] <= np.log10(max(ggcloudgrid)) and min(ggdiffusegrid) <=theta[4] <= 10. and min(favailablegrid) <=theta[5] <= max(favailablegrid):
        return 0.0
    return -np.inf

def lnprob(theta):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    print(theta,lp+lnlike(theta))
    return lp + lnlike(theta)

ndimw, nwalkers = 6, 50
pos = [np.array([np.random.uniform(np.log10(min(SNgrid)),np.log10(max(SNgrid)),1)[0],np.random.uniform(min(destroygrid),max(destroygrid),1)[0],(np.random.uniform(np.log10(min(fragmentgrid)),np.log10(max(fragmentgrid)),1)[0]),(np.random.uniform(np.log10(min(ggcloudgrid)),np.log10(max(ggcloudgrid)),1)[0]),np.random.uniform(min(ggdiffusegrid),10,1)[0],np.random.uniform(min(favailablegrid),max(favailablegrid),1)[0]]) for i in range(nwalkers)]

with terminating(Pool(processes=15)) as pool:
        sampler = emcee.EnsembleSampler(nwalkers, ndimw, lnprob, pool=pool)
        sampler.run_mcmc(pos, 50, progress=True)

samples = sampler.chain[:, :, :].reshape((-1, ndimw))

np.save("samples_dust_cont_2.npy",samples)
np.savetxt("samples_dust_cont_2.csv",samples, delimiter=",")

g1_mcmc, g2_mcmc, g3_mcmc, g4_mcmc , g5_mcmc, g6_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))

