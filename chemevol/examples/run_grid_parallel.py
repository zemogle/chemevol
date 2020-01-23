'''
This code has somewhat the same functionallity as the BulkEvolve class in evolve.py, but adds parallel processing.
It is based on the upload_csv and evolve_all functions in evolve.py.
It essentially reads in a precomputed csv file with all the inputs, and runs a model for each line in the csv.
The models are processed in parallel, with each model being processed independently, which leads to significant time gains when used with many cores.
'''


from chemevol import ChemModel
from astropy.table import Table
import time
from multiprocessing import Pool
import numpy as np
import os
import astropy.units as u
from astropy.io import ascii

'''
Function to run a single model for a given input set.
'''
def runModel(item):
    print("starting %s on ID: %s"%(item['name'],os.getpid()))
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
           'mgas' : all_results[:,1],
           'mstars' : all_results[:,2],
           'metallicity' : all_results[:,3],
           'mdust' : all_results[:,4],
           'dust_metals_ratio' : all_results[:,5],
           'sfr' : all_results[:,6],
           'dust_all' : all_results[:,7],
           'dust_stars' : all_results[:,8],
           'dust_ism' : all_results[:,9],
           'time_destroy' : all_results[:,10],
           'time_fragment' : all_results[:,11],
           'time_gg_diffuse' : all_results[:,12],
           'time_gg_cloud' : all_results[:,13],
           'mgas_outflow' : all_results[:,14],
           'mgas_recycled' : all_results[:,15],
           'mgas_inflow' : all_results[:,16],
           'mgas_IGM' : all_results[:,17],
           'mdust_IGM' : all_results[:,18],
           'mdust_diffuse' : all_results[:,19],
           'mdust_cloud' : all_results[:,20]}

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
        Miso=all_results[:,21+iso]
        params['M'+nameiso]=Miso
        Miso_IGM=all_results[:,21+iso+len(item['isotopes'])]
        params['M'+nameiso+'_IGM']=Miso_IGM
        paramsorder+=['M'+nameiso,'M'+nameiso+'_IGM']

    # write out to file based on 'name' identifier
    name = item['name']
    
    ascii.write(params,str(name+'.csv'),format='csv',names=paramsorder)

    # clear some memory
    del params
    del ch      
    print(name+" done")

def main():
    '''
    Read csv file with all the input parameters, format them into a list of input dictionaries.
    Process the models in parallel using these list of inputs.
    '''
    t = time.time()

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
        data = np.genfromtxt("modelGrid.csv", dtype=alttype,delimiter=',', autostrip=True, names=names)
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
        
        init_list.append(gal_data)

    procs = []

    pool = Pool(processes=4)  # set up parallel processing pool and number of processors to use.
    pool.map(runModel, init_list) # do the actual running of the models in parallel.

    print("total time for grid: %s s"%(time.time()-t))

if __name__ == "__main__":
    main()    
