@njit(cache=False)
def _stage18_pdrod_conservative(uc,us,ref,ref_len,sids,ids,w,broken,prev,
                               bcc,s01,sc1,horizon,vol,steel_cv,commit):
    fc=np.zeros((uc.shape[0],3),dtype=np.float64)
    fs=np.zeros((us.shape[0],3),dtype=np.float64)
    nf=0
    for k in range(ref_len.size):
        sid=sids[k]
        ax=ref[k,0]-us[sid,0]; ay=ref[k,1]-us[sid,1]; az=ref[k,2]-us[sid,2]
        for a in range(8):
            j=ids[k,a]; ww=w[k,a]
            ax+=ww*uc[j,0];ay+=ww*uc[j,1];az+=ww*uc[j,2]
        curr=np.sqrt(ax*ax+ay*ay+az*az)
        if curr<1e-15: raise ValueError('collapsed PDROD bond')
        st=(curr-ref_len[k])/ref_len[k]
        fail=broken[k] or st>sc1
        if fail:
            scalar=0.0;nf+=1
        else:
            c=bcc*(1.0-(ref_len[k]/horizon)**2)**2
            if st<=s01: scalar=c*st
            else: scalar=c*s01*(sc1-st)/(sc1-s01)
        # Appendix-C shaped bond force density q (N/m^3 IF bcc units N/m^6)
        # Convert to physical force using steel volume = lumped mass/rho.
        amp=scalar*vol*steel_cv[sid]/curr
        fx=amp*ax;fy=amp*ay;fz=amp*az
        fs[sid,0]+=fx;fs[sid,1]+=fy;fs[sid,2]+=fz
        for a in range(8):
            j=ids[k,a];ww=w[k,a]
            fc[j,0]-=ww*fx;fc[j,1]-=ww*fy;fc[j,2]-=ww*fz
        if commit:
            broken[k]=fail;prev[k]=st
    return fc,fs,nf

def conservative_interface(uc,us,commit=False):
    obj=fast_interface.original; pp=obj.params
    return _stage18_pdrod_conservative(uc,us,obj.reference,obj.lengths,
        obj.steel_ids,obj.indices,obj.weights,fast_interface.broken,
        fast_interface.prev_stretch,pp.bcc,pp.s01,pp.sc1,pp.horizon,
        pp.quadrature_volume,STEEL_CONTROL_VOLUME,commit)