fig,axs=plt.subplots(1,2,figsize=(12,4),constrained_layout=True)
midz=np.abs(xyz[:,2]-p.H/2)<dx/2+1e-10
sc=axs[0].scatter(xyz[midz,0]*1e3,xyz[midz,1]*1e3,c=counts[midz],s=5,cmap='viridis')
axs[0].set(xlabel='Beam length x (mm)',ylabel='Width y (mm)',title='Peridynamic family size — middle layer')
fig.colorbar(sc,ax=axs[0],label='Neighbor count')
axs[1].hist(LEN/dx,bins=np.arange(0.5,3.51,0.1),edgecolor='none')
axs[1].set(xlabel='Bond length / particle spacing',ylabel='Bond frequency',title='Full 3-D bond-length distribution')
plt.show()
