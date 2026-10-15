@njit(cache=False)
def _fast_pdrod_kernel(uc,us,ref,ref_len,st_ids,indices,weights,broken,prev,
                       bcc,s01,sc1,horizon,quadrature_volume,commit):
    Fcon=np.zeros((uc.shape[0],3),dtype=np.float64)
    Fsteel=np.zeros((us.shape[0],3),dtype=np.float64)
    failed_count=0
    for k in range(ref_len.size):
        sid=st_ids[k]
        dxv=ref[k,0]-us[sid,0]
        dyv=ref[k,1]-us[sid,1]
        dzv=ref[k,2]-us[sid,2]
        for a in range(8):
            idx=indices[k,a]
            wt=weights[k,a]
            dxv+=wt*uc[idx,0]
            dyv+=wt*uc[idx,1]
            dzv+=wt*uc[idx,2]
        length=np.sqrt(dxv*dxv+dyv*dyv+dzv*dzv)
        if length<=1.e-15:
            raise ValueError('collapsed PDROD interface bond')
        stretch=(length-ref_len[k])/ref_len[k]
        failed=broken[k] or stretch>sc1
        if failed:
            scalar=0.0
            failed_count+=1
        else:
            c=bcc*(1.0-(ref_len[k]/horizon)**2)**2
            if stretch<=s01:
                scalar=c*stretch
            else:
                scalar=c*s01*(sc1-stretch)/(sc1-s01)
        amp=scalar*quadrature_volume/length
        fx=amp*dxv;fy=amp*dyv;fz=amp*dzv
        Fsteel[sid,0]+=fx;Fsteel[sid,1]+=fy;Fsteel[sid,2]+=fz
        for a in range(8):
            idx=indices[k,a];wt=weights[k,a]
            Fcon[idx,0]-=wt*fx;Fcon[idx,1]-=wt*fy;Fcon[idx,2]-=wt*fz
        if commit:
            broken[k]=failed
            prev[k]=stretch
    return Fcon,Fsteel,failed_count

class FastInterface:
    def __init__(self, original):
        self.original=original
        self.broken=np.zeros_like(original.broken)
        self.prev_stretch=np.zeros_like(original.prev_stretch)
    def reset(self):
        self.broken[:]=False;self.prev_stretch[:]=0.
    def evaluate(self,uc,us,commit=False):
        it=self.original;p0=it.params
        return _fast_pdrod_kernel(uc,us,it.reference,it.lengths,it.steel_ids,it.indices,
            it.weights,self.broken,self.prev_stretch,p0.bcc,p0.s01,p0.sc1,
            p0.horizon,p0.quadrature_volume,commit)

fast_interface=FastInterface(interface8)
rng10=np.random.default_rng(120)
uc=rng10.normal(size=xyz.shape)*1.e-8
us=rng10.normal(size=steel_xyz.shape)*1.e-8
interface8.reset();fast_interface.reset()
a,b,_,_=interface8.evaluate(uc,us,commit=False)
c,d,nfailed=fast_interface.evaluate(uc,us,commit=False)
np.testing.assert_allclose(a,c,rtol=1.e-5,atol=2.e-4)
np.testing.assert_allclose(b,d,rtol=1.e-5,atol=2.e-4)
print('PASS: fast eight-node nonlinear interface agrees with Stage 8')