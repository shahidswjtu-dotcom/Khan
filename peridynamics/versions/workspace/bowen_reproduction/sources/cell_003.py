@dataclass(frozen=True)
class Parameters:
    L: float=0.600
    W: float=0.070
    H: float=0.100
    dx: float=0.005
    horizon_ratio: float=3.0
    support_span: float=0.500
    E_concrete: float=25e9
    nu_concrete: float=0.25
    rho_concrete: float=2400.0
    ft_concrete: float=4.6e6
    fc_concrete: float=60e6
    Gf_concrete: float=105.0
    E_steel: float=200e9
    rho_steel: float=7750.0
    fy_steel: float=482e6
    bar_diameter: float=0.006
    n_bars: int=2
    impactor_mass: float=38.4
    impactor_speed: float=10.8
    impactor_radius: float=0.060
    dt: float=3e-7

p = Parameters()
L, W, H, dx = p.L, p.W, p.H, p.dx
Nx, Ny, Nz = [round(t/dx) for t in (L,W,H)]
assert all(np.isclose(np.array([Nx,Ny,Nz])*dx, [L,W,H]))
assert (Nx, Ny, Nz)==(120,14,20)
print(f'Grid: {Nx} × {Ny} × {Nz} = {Nx*Ny*Nz:,} concrete material points')
print('Horizon =',p.horizon_ratio*dx,'m; reference time step =',p.dt,'s')
