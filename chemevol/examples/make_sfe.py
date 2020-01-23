'''
Script to make some of the .sfh and .sfe files.
these files have as columns:
- the time that has ellapsed since formation of the galaxy
- the star formation rate or star formation efficiency for .sfh and .sfe files respectively
- the inflow rate 
'''

import numpy as np

# We generate an array with the timesteps for the chemical evolution models.
# The chemev code will add steps at the start to better resolve the abrubt changes early in the galaxy's evolution.
# This is where the size of the timesteps is set for all the models.
# We extend this array to an age of 15 Gyr, but this age will not be reached as the models will stop after a time specified in the input file.
time=np.arange(0.03,15,0.03)*10**9

# We first give an example of how to generate a .sfh file.

# parameters to define the shape of delayed sfh taken from De Vis et al (2017b)
sfrmax=4.4
tau=6.9*10**9

# parameters to define the inflow rate
tinf=2*10**9
MG=4*10**10  # Mg can be overwritten within the inputs of the models

sfr=np.zeros(len(time))
inflows=np.zeros(len(time))

for i in range(len(time)):
	sfr[i]=time[i]/tau*2*np.exp(-time[i]/tau)*sfrmax/np.max(time/tau*2*np.exp(-time/tau))   #delayed sfh shape, normalised so that sfrmax is the maximum SFR value reached.
	inflows[i]=np.exp(-time[i]/tinf)/(tinf*(1-np.exp(-13.8*10**9/tinf)))

datfile=open("delayed.sfh","w")
for i in range(len(time)):
	datfile.write("%s %s %s \n"%(time[i],sfr[i],inflows[i]))
datfile.close()


# The following is to generate .sfe files used in De Vis et al (2020).

sfe=np.zeros(len(time))
inflows=np.zeros(len(time))

sfelog=[-8.5,-9,-9.5]
name=["fast.sfe","inflows.sfe","slow.sfe"]
for kk in range(3):
	for i in range(len(time)):
		sfe[i]=10**sfelog[kk]
		inflows[i]=np.exp(-time[i]/tinf)/(tinf*(1-np.exp(-13.8*10**9/tinf)))

	datfile=open(name[kk],"w")
	for i in range(len(time)):
		datfile.write("%s %s %s \n"%(time[i],sfe[i],inflows[i]))
	datfile.close()
