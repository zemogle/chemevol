'''
make some quick look up figures
'''
import matplotlib.pyplot as plt
def figure(time,gas,stars,metals,dust):
    plt.figure(figsize = (15,10))
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
    f2.semilogy(time,dust[:,0],color='purple',linestyle='-',linewidth=2,label='Dust')
    f2.semilogy(time,dust[:,1],color='purple',linestyle='--',linewidth=2,label='Stars')
    f2.semilogy(time,dust[:,2],color='purple',linestyle='-.',linewidth=2,label='ISM')
    f2.legend(frameon=False,loc='lower right',fontsize='11')
    f2.set_xlim(0.01,20)
    f2.set_ylim(5e5,5e9)
    f2.set_ylabel("Mass (Msun)", fontsize='16')
    f2.set_xlabel("Time (Gyrs)", fontsize='16')

    f3 = plt.subplot(2,3,3)
    f3.semilogy(time,metals,color='black',linewidth=2,label='Metals')
    f3.semilogy(time,dust[:,0],color='purple',linestyle='-',linewidth=2,label='Dust All')
    f3.semilogy(time,dust[:,1],color='purple',linestyle='--',linewidth=2,label='Dust Stars')
    f3.semilogy(time,dust[:,2],color='purple',linestyle='-.',linewidth=2,label='Dust ISM')
    f3.legend(frameon=False, loc='upper right', fontsize='11')
    f3.set_xlim(0.001,0.1)
    f3.set_ylim(1e2,1e9)
    f3.set_ylabel("Mass (Msun)", fontsize='16')
    f3.set_xlabel("Time (Gyrs)", fontsize='16')
    plt.show()
