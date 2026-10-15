import time
# Reconstructed Chapter 3 impact configuration (source: thesis).
HAMMER_MASS = p.impactor_mass
HAMMER_RADIUS = p.impactor_radius
HAMMER_V0 = -p.impactor_speed
# Assumed clearance to beam TOP envelope; top particle layer is z=H-dx/2.
INITIAL_GAP = 0.001  # m; provisional, not specified as an exact number by thesis
hammer_z_initial = float(np.max(xyz[:,2]) + HAMMER_RADIUS + INITIAL_GAP)
hammer_xy = np.array([p.L/2,p.W/2])
# Stiffness is an explicitly provisional penalty; diagnose dt and resulting penetration.
CONTACT_STIFFNESS = 2.0e8  # N/m PER PARTICLE; not Bowen's exact contact model
CONTACT_DAMPING = 0.0  # N s/m, set zero to expose energy drift
USE_GRAVITY = False    # constant impact velocity is imposed as INITIAL condition, not prescribed motion
# Support patches on the underside, at x = 50 and 550 mm, spanning beam width.
support_at_x = np.min(np.abs(xyz[:,0,None]-support_x[None,:]),axis=1) <= dx/2 + 1e-12
support_at_z = xyz[:,2] <= dx/2+1e-12
support_mask = support_at_x & support_at_z
SUPPORT_DOF_Z_ONLY=True  # roller-like vertical supports at BOTH supports: provisional simplification
# To eliminate horizontal rigid-body drift, fix x motion for ONE support particle and y for another.
anchor_candidates = np.flatnonzero(support_mask)
assert len(anchor_candidates) >= 2
anchor_x = anchor_candidates[0]; anchor_y = anchor_candidates[1]
print('Hammer center (initial):',hammer_z_initial,'m; initial gap:',INITIAL_GAP*1e3,'mm')
print('Support constrained bottom-layer particles:',int(support_mask.sum()))
print('Provisional particle contact stiffness:',CONTACT_STIFFNESS,'N/m')
assert abs(HAMMER_V0)==10.8 and len(anchor_candidates)>0
