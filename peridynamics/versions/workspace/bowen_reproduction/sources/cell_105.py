# Bridge the source-shaped Stage 13 class into Stage 12/10 beam arrays.
# IMPORTANT: sample_pos, sample_idx, sample_w and DEMO_INTERFACE are PROVISIONAL
# from earlier setup. Appendix C does not give their original Chapter 3 values.
SOURCE_C_PARAMS=AppendixCParameters(
    bcc=DEMO_INTERFACE.bcc,
    s01=DEMO_INTERFACE.s01,
    sc1=DEMO_INTERFACE.sc1,
    delta=DEMO_INTERFACE.horizon,
    vol=DEMO_INTERFACE.quadrature_volume,
)
appendix_c_beam=AppendixCPublishedKernel(
    steel_xyz,xyz,sample_steel,sample_idx,sample_w,sample_pos,SOURCE_C_PARAMS)
assert len(appendix_c_beam.st)==len(sample_steel)
print('Stage 14 Appendix C beam samples:',len(sample_steel),' UNVERIFIED constants:',SOURCE_C_PARAMS)
