# Explicitly demonstrate the known key behavioral difference after a bond breaks.
appendix_c_beam.reset();fast_interface.reset()
# Deform only the steel node associated with one sample, temporarily;
# many interactions on this node may fail, which is expected.
uf=np.zeros_like(xyz);sf=np.zeros_like(steel_xyz)
first=int(sample_steel[0]);sf[first,0]=1.0e-3
appendix_c_beam.evaluate(uf,sf,mode='literal',commit=True)
fast_interface.evaluate(uf,sf,commit=True)
assert np.count_nonzero(appendix_c_beam.failfail==0)>0
uf.fill(0.);sf.fill(0.);sf[first,0]=1.0e-7
literal_after=appendix_c_beam.evaluate(uf,sf,mode='literal')['steel_force'][first].copy()
corrected_after=appendix_c_beam.evaluate(uf,sf,mode='irreversible_corrected')['steel_force'][first].copy()
fc_prev,fs_prev,_=fast_interface.evaluate(uf,sf,False)
print('Broken interface samples:',np.count_nonzero(appendix_c_beam.failfail==0))
print('Re-loading nodal steel force, literal (N):',literal_after)
print('Re-loading nodal steel force, history-corrected (N):',corrected_after)
print('Previously accelerated irreversible code (N):',fs_prev[first])
appendix_c_beam.reset();fast_interface.reset()
