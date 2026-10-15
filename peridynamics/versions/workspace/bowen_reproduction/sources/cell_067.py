pt=np.array([[p.L/2+0.005, p.W/2+0.002, 0.0]])
z_test=HAMMER_RADIUS-0.002
fcont,fh,ep,nn,dp=sphere_contact(pt,z_test)
assert nn==1
for d in range(3):
    eps=1e-7
    pplus=pt.copy();pminus=pt.copy()
    pplus[0,d]+=eps;pminus[0,d]-=eps
    plus=sphere_contact(pplus,z_test)[2]
    minus=sphere_contact(pminus,z_test)[2]
    grad=(plus-minus)/(2*eps)
    np.testing.assert_allclose(fcont[0,d],-grad,atol=1e-3,rtol=1e-5)
print('PASS: finite-difference contact force equals negative gradient of potential')
