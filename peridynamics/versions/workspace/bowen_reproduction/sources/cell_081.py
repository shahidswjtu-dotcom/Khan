# Stage 8 uses the existing Numba conservative assembly and Stage 7 constitutive functions.
# Concrete forces are integrated over unique (i,j) bonds; PDROD returns actual nodal forces.
class SourceShapeConcrete(DiscretePDConcrete):
    def evaluate(self, displacement, dt, commit=False, override_rate=None):
        deformed=self.xyz+np.asarray(displacement,float)
        d=deformed[self.J]-deformed[self.I]
        cur=np.linalg.norm(d,axis=1)
        if np.any(cur<1e-15): raise FloatingPointError('Collapsed concrete bond')
        stretch=(cur-self.lengths)/self.lengths
        rate=(stretch-self.prev_stretch)/dt if override_rate is None else np.broadcast_to(override_rate,stretch.shape)
        scalar,failed=bowen_concrete_force_with_modulus_DIF(stretch,rate,self.lengths,self.broken)
        amplitude=scalar*(V**2)*self.weights
        f=accumulate_unique_pairs(self.I,self.J,d/cur[:,None],amplitude,len(self.xyz))
        if commit:
            self.broken|=failed
            self.prev_stretch[:]=stretch
        return f,stretch,failed

concrete8=SourceShapeConcrete(xyz,I,J,LEN,PAIR_WEIGHT)
interface8=AppendixCPDROD(steel_xyz,sample_pos,sample_idx,sample_w,sample_steel,DEMO_INTERFACE)

def reset_stage8():
    concrete8.broken.fill(False)
    concrete8.prev_stretch.fill(0.)
    interface8.reset()

def source_shaped_forces(uc,us,zh,dt,commit=False):
    pd_force, stretch, damaged=concrete8.evaluate(uc,dt,commit=commit)
    interface_c, interface_s, iface_stretch, iface_failed=interface8.evaluate(uc,us,commit=commit)
    rod_f=rod_forces(us)
    contact_c,contact_h,contact_U,ncontact,maxpen=sphere_contact(xyz+uc,zh)
    concrete_force=pd_force+interface_c+contact_c
    steel_force=rod_f+interface_s
    return concrete_force,steel_force,float(contact_h[2]),contact_U,ncontact,maxpen,int(damaged.sum()),int(iface_failed.sum())

print('Stage 8 integrated source-shaped forces READY')
print('WARNING: interface constants and contact law remain PROVISIONAL')