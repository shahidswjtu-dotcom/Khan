"""Independent numerical checks for the supplied ORIGINAL Fortran benchmark.

Reads the original source and completed nodal output; never changes the solver.
The exact reference rod branch is reproduced in an isolated two-node calculation.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'Console2.original.f90'
DX = 19e-3
VOL = DX**3
VOL1 = DX*530.66e-6
RHO_C = 2500.
RHO_S = 7850.
DAMP = 3e6
DT = 3e-6
CAP = 416.*530.66
E_S = 208e9
A1 = np.pi*(26e-3)**2/4

def trilinear_weights(xi):
    corners = np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],
                        [-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]])
    return np.prod(1+corners*np.asarray(xi), axis=1)/8

def original_yielded_rod_density():
    f = np.zeros((2,3))
    f[0,0] += CAP / VOL1
    f[0,1] += 0.
    f[0,2] += 0.
    f[1,0] -= CAP / VOL1
    # Literal original lines 746 and 747, not a hypothesized formulation.
    f[1,1] = f[1,0] + 0.
    f[1,2] = f[1,0] + 0.
    return f

def verify_source():
    src = SOURCE.read_text()
    for exact in ['pforce22(ccnode,2)=pforce22(ccnode,1)+0.0d0',
                  'pforce22(ccnode,3)=pforce22(ccnode,1)+0.0d0',
                  'vol1=dx*(530.66d-6)', 'cn=3.0d6']:
        assert exact in src
    return src

def run_checks():
    verify_source()
    force_density = original_yielded_rod_density()
    physical = force_density*VOL1
    assert np.allclose(physical.sum(axis=0), [0,-CAP,-CAP])
    corrected = physical.copy(); corrected[1,1:] = 0.
    assert np.allclose(corrected.sum(axis=0), 0.)

    # One weighted PDROD link. Virtual interpolation and force transfer are adjoint
    # in nodal density space. Physical volume factors decide force/work balance.
    weights = trilinear_weights([.23,-.41,.61])
    assert np.all(weights > 0) and np.isclose(weights.sum(),1)
    dforce = np.array([1e9,-2e9,3e9])  # N/m^3, arbitrary nonzero test vector
    concrete_density = -weights[:,None]*dforce
    assert np.allclose(dforce+concrete_density.sum(axis=0),0,atol=1e-6)
    net = dforce*VOL1+concrete_density.sum(axis=0)*VOL
    assert np.allclose(net, dforce*(VOL1-VOL))
    # With the author's density on steel preserved, multiply each concrete density
    # by VOL1/VOL. This restores both physical action-reaction and virtual work.
    concrete_corrected_density = concrete_density*(VOL1/VOL)
    assert np.allclose(dforce*VOL1+concrete_corrected_density.sum(axis=0)*VOL,0,atol=1e-10)
    rng = np.random.default_rng(1749)
    velocity_s = rng.normal(size=3)
    velocity_c = rng.normal(size=(8,3))
    power_expected = VOL1*np.dot(dforce, velocity_s-weights@velocity_c)
    power_fixed = VOL1*np.dot(dforce,velocity_s)+VOL*np.sum(concrete_corrected_density*velocity_c)
    assert np.isclose(power_expected,power_fixed,rtol=1e-14)

    # The active time update is semi-implicit Euler with explicit linear damping:
    # v[n+1]=(1-gamma*dt)*v[n]-lambda*dt*u[n], u[n+1]=u[n]+dt*v[n+1].
    # Its scalar oscillator stability condition is lambda*dt^2+2*gamma*dt < 4.
    gamma_s = DAMP/RHO_S
    gamma_c = DAMP/RHO_C
    lambda_rod = 4*E_S*A1/(RHO_S*VOL1*DX)  # infinite axial chain high-frequency end
    stable_dt_rod = (np.sqrt(gamma_s**2+4*lambda_rod)-gamma_s)/lambda_rod
    v_equilibrium = force_density[1]/DAMP
    assert np.isclose(v_equilibrium[1], -7298.245614035087)
    # An isolated clamped position retains nonzero constrained velocity because
    # the source zeros displacement only. It approaches f_density / cn.
    v_clamped = 0.
    for _ in range(20000):
        v_clamped += DT*(1e6-DAMP*v_clamped)/RHO_C
        u_unconstrained = v_clamped*DT
        u_clamped = 0.
    assert u_clamped == 0 and np.isclose(v_clamped,1e6/DAMP)

    # The reference-linear truss matrix is not objective under finite rigid rotation.
    angle = .1
    reference_linear_rotation_force = E_S*A1*(np.cos(angle)-1)
    assert abs(reference_linear_rotation_force)>5e5

    final_path = ROOT/'runs/literal_100000/final_state.txt'
    corroboration = None
    if final_path.exists():
        # Columns: index, reference xyz, displacement xyz, velocity xyz, damage.
        state = np.loadtxt(final_path)
        target = v_equilibrium[1]
        vy_match = np.isclose(state[:,8], target, rtol=1e-10, atol=1e-7)
        vz_match = np.isclose(state[:,9], target, rtol=1e-10, atol=1e-7)
        corroboration = {
            'node_2_displacement_m': state[1,4:7].tolist(),
            'node_2_velocity_m_s': state[1,7:10].tolist(),
            'predicted_wrong_transverse_terminal_velocity_m_s': float(target),
            'node_2_y_velocity_difference_m_s': float(state[1,8]-target),
            'node_2_z_velocity_difference_m_s': float(state[1,9]-target),
            'nodes_with_same_y_and_z_wrong_terminal_velocity': int(np.count_nonzero(vy_match&vz_match)),
            'maximum_displacement_norm_m': float(np.linalg.norm(state[:,4:7],axis=1).max()),
        }
    return {
        'source': str(SOURCE),
        'checks_passed': ['literal_yielded_rod_transverse_force_bug',
                          'corrected_axial_rod_has_zero_internal_resultant',
                          'PDROD_density_space_adjoint_transfer',
                          'physical_PDROD_volume_normalization_counterexample',
                          'volume_corrected_interface_force_and_virtual_work',
                          'displacement_only_support_retains_constrained_velocity',
                          'reference_linear_rod_finite_rotation_counterexample'],
        'constants': {'dx_m':DX, 'dt_s':DT, 'concrete_volume_m3':VOL,
                      'steel_nodal_volume_m3':VOL1, 'volume_ratio':VOL1/VOL,
                      'original_rod_A1_m2':float(A1), 'volume_area_m2':530.66e-6,
                      'cap_force_N':CAP, 'steel_mass_kg':RHO_S*VOL1,
                      'concrete_mass_kg':RHO_C*VOL},
        'yielded_rod': {'physical_forces_N':physical.tolist(),
                       'net_internal_force_N':physical.sum(axis=0).tolist(),
                       'wrong_right_node_terminal_velocity_m_s':v_equilibrium.tolist()},
        'interface': {'weight_sum':float(weights.sum()), 'test_density_N_m3':dforce.tolist(),
                      'physical_internal_resultant_N':net.tolist(),
                      'required_density_reaction_factor_preserving_steel_force':VOL1/VOL,
                      'physical_normalization_is_conditional_on_declared_vol1_steel_mass':True,
                      'alternative_uniform_vol_interpretation':'Conserves interface force but scales physical rod stiffness/cap by vol/vol1 and gives a steel nodal mass inconsistent with the rod area.'},
        'time_integration': {'gamma_concrete_per_s':gamma_c,'gamma_steel_per_s':gamma_s,
                             'gamma_concrete_times_dt':gamma_c*DT,
                             'gamma_steel_times_dt':gamma_s*DT,
                             'isolated_axial_chain_lambda_per_s2':float(lambda_rod),
                             'isolated_axial_chain_critical_dt_s':float(stable_dt_rod),
                             'isolated_axial_chain_stability_lhs':float(lambda_rod*DT**2+2*gamma_s*DT),
                             'full_mesh_and_contact_stability_not_established':True,
                             'damping_recurrence_is_stable_at_original_dt':True},
        'support_counterexample':{'clamped_displacement_m':u_clamped,
                                 'retained_constrained_velocity_m_s':v_clamped,
                                 'last_unprojected_increment_m':u_unconstrained},
        'rod_rotation': {'rigid_rotation_rad':angle,
                         'exact_finite_length_strain':0.,
                         'reference_linear_axial_force_N':float(reference_linear_rotation_force)},
        'completed_run_corroboration':corroboration,
    }

if __name__=='__main__':
    result = run_checks()
    output = Path(__file__).with_suffix('.json')
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
