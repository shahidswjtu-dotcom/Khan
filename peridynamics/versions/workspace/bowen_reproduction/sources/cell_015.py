metadata={
 'source':'Tu Bowen 2025 thesis, Chapter 3 section 3.1.1 (source-based geometry setup)',
 'parameters':asdict(p), 'number_particles':int(N),'grid_shape':[Nx,Ny,Nz],
 'support_x_m':support_x.tolist(),
 'assumptions':{'bar_y_centers_m':bar_y.tolist(),'bar_z_m':bar_z,
                'bar_x_endpoints_m':bar_x.tolist(),'illustrative_impactor_clearance_m':illustrative_gap,
                'grid_origin':"cell-centered; author's grid origin not explicitly specified",
                'supports':'markers only; no constraint/contact formulation implemented',
                'steel':'marker nodes only; no PDROD mechanics implemented',
                'impactor':'geometry only; no contact implementation'},
 'reference_targets_not_simulated':{'peak_force_N':359800.0,'peak_midspan_deflection_m':0.02270}
}
np.savez_compressed(OUT/'Bowen_Chapter3_geometry.npz', concrete_xyz=xyz,
                    support_x=support_x, steel_marker_xyz=steel_xyz, impactor_center=impact_center)
with (OUT/'Bowen_Chapter3_geometry_metadata.json').open('w',encoding='utf-8') as f: json.dump(metadata,f,indent=2)
print('Saved:',[str(v) for v in OUT.iterdir()])
