'''
This makes a quick look up figure and writes to data
'''
import matplotlib.pyplot as plt

def writedata(time, mgas, mstars, sfr, ssfr, mdust, metalmass, metallicity, gasfraction):
    # write data to a file for plotting
    f = open('results.dat', 'w')
    text = "# Time(Gyr) Mg(Msol) Ms(Msol) SFR(Msol/yr) SSFR(/yr) Md(Msol) M_Z(Msol) Z, f_g \n"
    f.write(text)
    for i in range(0,len(mgas)):
        t = "%s %.2f %.2f %.4f %s %.2f %.2f %s %.4f \n" % \
            (time[i], mgas[i], mstars[i], sfr[i], ssfr[i], mdust[i], \
            metalmass[i], metallicity[i], gasfraction[i])
        f.write(t)
    f.close()

def figure(time,gas,stars,metals,metallicity,dust,dust_metals_ratio,gasfraction,dust_in,timescale):
    # plot some figures to check out results
    plt.figure(figsize = (20,8))
    f1 = plt.subplot(2,3,1)
    f1.semilogy(time,gas,color='black',linestyle='-',linewidth=2,label='Gas')
    f1.semilogy(time,stars,color='black',linestyle='--',linewidth=2,label='Stars')
    f1.set_xlim(0.01,20)
    f1.set_ylim(5e8,5e11)
    f1.set_ylabel("Mass (Msun)", fontsize='16')
    f1.set_xlabel("Time (Gyrs)", fontsize='16')
    f1.legend(frameon=False,loc='lower right',fontsize='11')

    f2 = plt.subplot(2,3,2)
    f2.semilogy(time,metals,color='black',linewidth=2,label='Metals')
    f2.semilogy(time,dust,color='purple',linestyle='-',linewidth=2,label='Dust All')
    f2.legend(frameon=False,loc='lower right',fontsize='11')
    f2.set_xlim(0.01,20)
    f2.set_ylim(1e6,5e9)
    f2.set_ylabel("Mass (Msun)", fontsize='16')
    f2.set_xlabel("Time (Gyrs)", fontsize='16')

    f3 = plt.subplot(2,3,3)
    f3.semilogy(time,metals,color='black',linewidth=2,label='Metals')
    f3.semilogy(time,dust,color='purple',linestyle='-',linewidth=2,label='Dust All')
    f3.legend(frameon=False, loc='upper right', fontsize='11')
    f3.set_xlim(0.001,0.1)
    f3.set_ylim(1e2,1e9)
    f3.set_ylabel("Mass (Msun)", fontsize='16')
    f3.set_xlabel("Time (Gyrs)", fontsize='16')

    f4 = plt.subplot(2,3,4)
    f4.semilogy(time,dust_in[:,0],color='purple',linestyle='-',linewidth=2,label='Dust In (Stars+ISM)')
    f4.semilogy(time,dust_in[:,1],color='purple',linestyle='-.',linewidth=2,label='Dust Stars')
    f4.semilogy(time,dust_in[:,2],color='purple',linestyle='--',linewidth=2,label='Dust ISM')
    f4.legend(frameon=False, loc='lower right', fontsize='11')
    f4.set_xlim(0.01,5)
    f4.set_ylim(5e3,1e10)
    f4.set_ylabel("Mass (Msun)", fontsize='16')
    f4.set_xlabel("Time (Gyrs)", fontsize='16')

    f5 = plt.subplot(2,3,5)
    f5.semilogy(gasfraction,metallicity,color='black',linewidth=2,label='Metallicity')
    f5.semilogy(gasfraction,dust_metals_ratio,color='purple',linestyle='-',linewidth=2,label='Dust to Metals')
    f5.legend(frameon=False, loc='lower right', fontsize='11')
    f5.set_xlim(1.0,0.1)
    f5.set_ylim(1e-5,1.5)
    f5.set_ylabel("Fraction", fontsize='16')
    f5.set_xlabel("Gas Fraction", fontsize='16')

    f6 = plt.subplot(2,3,6)
    f6.plot(time,timescale[:,0],color='black',linewidth=2,label='destruction')
    f6.plot(time,timescale[:,1],color='purple',linewidth=2,label='grain growth')
    f6.legend(frameon=False, loc='upper right', fontsize='11')
    f6.set_xlim(0.01,20)
    f6.set_ylim(0.01,2)
    f6.set_ylabel("Timescale (Gyr)", fontsize='16')
    f6.set_xlabel("Time (Gyr)", fontsize='16')
    plt.show()
