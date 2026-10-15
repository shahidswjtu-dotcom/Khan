radii=np.linspace(0,delta,20001)
c_const=micromodulus(radii,improved=False)
c_impr=micromodulus(radii,improved=True)
# Relative fracture energy = integral c(r)*r^3 dr, common s0,sc factors cancel.
I_const=np.trapezoid(c_const*radii**3,radii)
I_impr=np.trapezoid(c_impr*radii**3,radii)
ratio=I_impr/I_const
print(f'Weighted radial integral ratio (improved / constant) = {ratio:.6f}')
print(f'Implied energy scale if Gf=105 N/m calibrated for constant c: {p.Gf_concrete*ratio:.2f} N/m')
print('NOTE: This is an analytic/numerical consistency diagnostic, not a measured beam fracture energy.')
assert 0.0<ratio<5.0
