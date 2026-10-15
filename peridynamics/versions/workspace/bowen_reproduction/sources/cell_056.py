CONCRETE_NODE_MASS = p.rho_concrete*dx**3
# Lump half of each element's mass to each adjacent steel node.
STEEL_MASS = np.zeros(Nsteel)
rod_mass_per_element=p.rho_steel*steel_A*steel_lengths
np.add.at(STEEL_MASS,rod_pairs[:,0],0.5*rod_mass_per_element)
np.add.at(STEEL_MASS,rod_pairs[:,1],0.5*rod_mass_per_element)
assert np.all(STEEL_MASS>0)
print('Concrete mass (kg):',N*CONCRETE_NODE_MASS)
print('Steel mass (kg):',STEEL_MASS.sum())
print('Hammer mass (kg):',HAMMER_MASS)
print('Initial hammer kinetic energy (J):',0.5*HAMMER_MASS*HAMMER_V0**2)
