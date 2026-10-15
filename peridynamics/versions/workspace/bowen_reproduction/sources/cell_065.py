center_xy=hammer_xy
radial2=(xyz[:,0]-center_xy[0])**2+(xyz[:,1]-center_xy[1])**2
eligible=radial2<HAMMER_RADIUS**2
# For a sphere above the beam, the first point-contact height is z_i+sqrt(R^2-r^2).
first_center_height=np.max(xyz[eligible,2]+np.sqrt(HAMMER_RADIUS**2-radial2[eligible]))
travel=max(hammer_z_initial-first_center_height,0.0)
first_t=travel/abs(HAMMER_V0)
first_step=int(np.ceil(first_t/p.dt))
print(f'First geometric contact with undeformed points: t={first_t*1e6:.3f} microseconds; step {first_step}')
print(f'First-point equivalent travel: {travel*1e3:.5f} mm')
assert first_step>=1
Fc0,Fh0,Uc0,n0,pen0=sphere_contact(xyz,hammer_z_initial)
assert n0==0 and np.linalg.norm(Fh0)==0.0
print('PASS: initial contact is zero; predicted geometric first-contact step computed')
