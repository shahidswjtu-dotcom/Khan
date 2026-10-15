"""Explicit probe locations in reference coordinates; metres internally."""
import numpy as np

def probes(m,u):
    def at(q):
        q=np.asarray(q,float)
        a=q/m.dx-.5
        shape=np.array([m.Nx,m.Ny,m.Nz])
        if np.any(a < -1e-12) or np.any(a > shape-1+1e-12):
            raise ValueError("Probe outside particle-centre grid")
        base=np.minimum(np.maximum(np.floor(a).astype(int),0),shape-2)
        frac=a-base
        ids=[];w=[]
        for i in (0,1):
            for j in (0,1):
                for k in (0,1):
                    ids.append(np.ravel_multi_index(base+[i,j,k],shape))
                    w.append(np.prod(np.where([i,j,k],frac,1-frac)))
        ids=np.array(ids);w=np.array(w)
        return float(np.dot(w,u[ids,2]))
    # No extrapolation to physical surface. Report centre-layer offset explicitly.
    bottom=at([m.p.L/2,m.p.W/2,m.dx/2])
    supports=[at([x,m.p.W/2,m.dx/2]) for x in m.support_x]
    mid=(abs(m.xyz[:,0]-m.p.L/2)<=m.dx/2+1e-12)&(m.xyz[:,2]>=m.p.H/2)
    return dict(legacy_half_depth_mm=-1000*np.mean(u[mid,2]),bottom_center_layer_mm=-1000*bottom,
        support_relative_bottom_mm=-1000*(bottom-np.mean(supports)),
        neutral_axis_mm=-1000*at([m.p.L/2,m.p.W/2,m.p.H/2]),
        top_center_layer_mm=-1000*at([m.p.L/2,m.p.W/2,m.p.H-m.dx/2]))
