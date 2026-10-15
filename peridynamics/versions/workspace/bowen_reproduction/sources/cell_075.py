class AppendixCPDROD:
    """Source-shaped translation of the Appendix C nonlinear steel-to-concrete force.

    Sample is a virtual material point whose displacement is trilinearly interpolated
    from eight concrete points. Forces on the real concrete nodes are -N_a*F_steel.
    Persistent broken flags prevent crack healing on subsequent evaluations.
    """
    def __init__(self, steel_reference, virtual_reference, indices, weights, steel_ids, params):
        self.xs=np.asarray(steel_reference,dtype=float)
        self.xv=np.asarray(virtual_reference,dtype=float)
        self.indices=np.asarray(indices,dtype=np.int32)
        self.weights=np.asarray(weights,dtype=float)
        self.steel_ids=np.asarray(steel_ids,dtype=np.int32)
        self.params=params
        self.reference=self.xv-self.xs[self.steel_ids]
        self.lengths=np.linalg.norm(self.reference,axis=1)
        if self.indices.ndim!=2 or self.indices.shape[1]!=8: raise ValueError('eight interpolation weights required')
        if np.any(self.lengths<=0): raise ValueError('zero-length steel–concrete interface bond')
        if np.any(self.lengths>=params.horizon): raise ValueError('all coupling samples must lie inside horizon')
        np.testing.assert_allclose(self.weights.sum(axis=1),1.0,atol=1e-12)
        if not 0<params.s01<params.sc1: raise ValueError('invalid interface stretches')
        self.broken=np.zeros(len(self.lengths),dtype=bool)
        self.prev_stretch=np.zeros(len(self.lengths))

    def evaluate(self, concrete_u, steel_u, commit=False):
        concrete_u=np.asarray(concrete_u,dtype=float)
        steel_u=np.asarray(steel_u,dtype=float)
        interp=np.sum(concrete_u[self.indices]*self.weights[:,:,None],axis=1)
        deformed=self.reference+interp-steel_u[self.steel_ids]
        new_lengths=np.linalg.norm(deformed,axis=1)
        if np.any(new_lengths<=1e-16): raise FloatingPointError('collapsed interface bond')
        stretch=(new_lengths-self.lengths)/self.lengths
        prm=self.params
        c=prm.bcc*(1-(self.lengths/prm.horizon)**2)**2
        scalar=np.where(stretch<=prm.s01, c*stretch,
                    c*prm.s01*np.clip((prm.sc1-stretch)/(prm.sc1-prm.s01),0,1))
        newly_broken=stretch>prm.sc1
        effective_broken=self.broken | newly_broken
        scalar=np.where(effective_broken,0.,scalar)
        # Appendix C expresses the steel interaction as displacement toward the virtual point.
        sample_force=(scalar*prm.quadrature_volume)[:,None]*deformed/new_lengths[:,None]
        fsteel=np.zeros_like(steel_u)
        fconcrete=np.zeros_like(concrete_u)
        np.add.at(fsteel,self.steel_ids,sample_force)
        for a in range(8):
            np.add.at(fconcrete,self.indices[:,a],-self.weights[:,a,None]*sample_force)
        if commit:
            self.broken |= newly_broken
            self.prev_stretch[:]=stretch
        return fconcrete,fsteel,stretch,effective_broken

    def reset(self):
        self.broken.fill(False)
        self.prev_stretch.fill(0.)

pdrod_c=AppendixCPDROD(steel_xyz,sample_pos,sample_idx,sample_w,sample_steel,DEMO_INTERFACE)
print('Appendix C interactions:',len(pdrod_c.lengths),'initial broken:',pdrod_c.broken.sum())
