fig,axs=plt.subplots(3,1,figsize=(11,10),constrained_layout=True)
# Side view
a=axs[0]; a.add_patch(plt.Rectangle((0,0),L,H,facecolor='#dae9f2',edgecolor='black'))
for sx in support_x: a.plot(sx,-0.007,marker='^',color='#c45a25',ms=12)
a.plot([bar_x[0],bar_x[-1]],[bar_z,bar_z],color='#a53030',lw=3,label='steel rod height (assumed)')
theta=np.linspace(0,2*np.pi,300);a.plot(impact_center[0]+p.impactor_radius*np.cos(theta),impact_center[2]+p.impactor_radius*np.sin(theta),'k--',lw=1.2,label='spherical head projection')
a.set(xlim=(-0.02,0.62),ylim=(-0.025,0.23),xlabel='x (m)',ylabel='z (m)',title='Side view: beam, supports and impactor (schematic)');a.legend(loc='upper right')
# Top view
a=axs[1];a.add_patch(plt.Rectangle((0,0),L,W,facecolor='#dae9f2',edgecolor='black'))
for yy in bar_y: a.plot(bar_x,[yy,yy],color='#a53030',lw=3)
a.plot(*impact_xy,marker='x',color='black',ms=10,mew=2);a.set(xlim=(-0.02,0.62),ylim=(-0.008,0.078),xlabel='x (m)',ylabel='y (m)',title='Plan view: transverse bar positions are ASSUMED');a.set_aspect('equal',adjustable='box')
# Cross section
a=axs[2];a.add_patch(plt.Rectangle((0,0),W,H,facecolor='#dae9f2',edgecolor='black'))
for yy in bar_y:a.add_patch(plt.Circle((yy,bar_z),p.bar_diameter/2,color='#a53030'))
a.set(xlim=(-0.01,0.08),ylim=(-0.012,0.115),xlabel='y (m)',ylabel='z (m)',title='Cross-section: two Ø6 rods, centers to verify');a.set_aspect('equal')
fig.savefig(OUT/'Bowen_reference_geometry.png',dpi=220);plt.show()
