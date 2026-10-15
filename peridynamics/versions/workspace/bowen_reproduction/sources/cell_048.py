fig,axs=plt.subplots(1,2,figsize=(12,4),constrained_layout=True)
axs[0].scatter(xyz[::30,0]*1e3,xyz[::30,2]*1e3,s=2,alpha=.12,label='Concrete particles')
for iy,yy in enumerate(bar_y):
    q=steel_xyz[np.abs(steel_xyz[:,1]-yy)<1e-10]
    axs[0].plot(q[:,0]*1e3,q[:,2]*1e3,lw=2,label=f'Steel rod {iy+1}')
axs[0].set(xlabel='Beam x (mm)',ylabel='Beam z (mm)',title='Provisional reinforcement layout')
axs[0].legend()
axs[1].hist(sample_w.ravel(),bins=np.linspace(0,1,26))
axs[1].set(xlabel='Eight-node interpolation weight',ylabel='Frequency',title='Interface interpolation weights')
plt.show()
