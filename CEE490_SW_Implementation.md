# CEE 490 — Graduate Level 2D Seismic Wave Propagation
## Full Elastic Wave Equation on a Staggered Grid with Site Amplification Analysis

**Reference:** Moczo, Kristek & Halada (2004) — *The Finite-Difference Method for Seismologists*
**Additional References:** Virieux (1986), Levander (1988), Cerjan et al. (1985)
**Deadline:** May 20 | **Weight:** 45% of final grade
**Language:** Python 3 | **Libraries:** NumPy, Matplotlib, SciPy

---

## What Makes This Graduate Level

Most students implement the acoustic wave equation — a single scalar PDE for pressure.
You implement the **full elastic wave equation** — a coupled system of five PDEs for
two velocity components and three stress components. This produces both P waves AND
S waves simultaneously, mode conversion at interfaces, and Rayleigh surface waves.

On top of that you add:

| Feature | Why It Matters |
|---------|---------------|
| Full elastic wave equation | P + S waves, mode conversion, surface waves |
| Staggered grid (Virieux 1986) | Industry standard, >2000 citations |
| Three numerical schemes + comparison | Original scientific evaluation |
| Convergence analysis (log-log) | Mathematical proof of correctness |
| Dispersion analysis | Quantifies numerical error analytically |
| Site amplification application | Real civil engineering problem |
| Viscoelastic attenuation | Realistic energy dissipation in rock |
| Performance analysis | Accuracy per unit computational cost |
| Journal paper format report | Scientific communication standard |

---

## The Physics — Why Elastic Is More Powerful Than Acoustic

### Acoustic Wave Equation (what most students do)
```
d2p/dt2 = c^2 * (d2p/dx2 + d2p/dz2)
```
One field. One wave type. Pressure only. No shear.

### Full Elastic Wave Equation (what you do)
```
rho * dv_x/dt = dtau_xx/dx + dtau_xz/dz  + f_x
rho * dv_z/dt = dtau_xz/dx + dtau_zz/dz  + f_z
dtau_xx/dt   = (lam + 2*mu) * dv_x/dx + lam * dv_z/dz
dtau_zz/dt   =  lam * dv_x/dx + (lam + 2*mu) * dv_z/dz
dtau_xz/dt   =  mu * (dv_x/dz + dv_z/dx)
```

Five coupled fields: v_x, v_z, tau_xx, tau_zz, tau_xz

Two wave speeds:
```
P-wave speed:   vp = sqrt((lam + 2*mu) / rho)   [compressional]
S-wave speed:   vs = sqrt(mu / rho)              [shear, slower]
```

This produces physically richer results:
- P waves arrive first (faster)
- S waves arrive second (slower, more damaging in earthquakes)
- Mode conversion at interfaces (P partially becomes S and vice versa)
- Rayleigh surface waves emerge naturally along the free surface

---

## File Structure

```
seismic_wave_2d/
├── main.py                  <- Master driver: runs everything
├── parameters.py            <- All physical and numerical parameters
├── elastic_grid.py          <- Staggered grid field initialization
├── elastic_fd.py            <- Full elastic FD stencils (2nd and 4th order)
├── acoustic_fd.py           <- Acoustic FD for comparison
├── rk4_update.py            <- RK4 time integrator
├── source.py                <- Ricker wavelet, moment tensor source
├── boundaries.py            <- Sponge absorbing boundaries
├── velocity_models.py       <- All geological models
├── simulation_runner.py     <- Reusable simulation engine
├── analysis.py              <- Convergence + dispersion + performance
├── site_amplification.py    <- Site response analysis
├── visualize.py             <- All professional plotting functions
├── figures/                 <- All output figures
└── report/                  <- Written report
```

---

## Day 1 — Parameters, Staggered Grid & Elastic Equations

### Step 1.1 — Environment Setup

```bash
mkdir seismic_wave_2d && cd seismic_wave_2d

pip install numpy matplotlib scipy pillow

# Create all files
touch main.py parameters.py elastic_grid.py elastic_fd.py
touch acoustic_fd.py rk4_update.py source.py boundaries.py
touch velocity_models.py simulation_runner.py analysis.py
touch site_amplification.py visualize.py
mkdir figures report
```

---

### Step 1.2 — Physical and Numerical Parameters

Create `parameters.py`:

```python
import numpy as np

# =============================================
# ELASTIC MEDIUM PARAMETERS
# =============================================
rho   = 2200.0      # Density [kg/m3]
vp    = 3000.0      # P-wave velocity [m/s]
vs    = 1732.0      # S-wave velocity [m/s]  (vs = vp/sqrt(3) for Poisson solid)

# Lame parameters derived from vp, vs, rho
mu    = rho * vs**2                  # Shear modulus [Pa]
lam   = rho * vp**2 - 2 * mu        # First Lame parameter [Pa]

# =============================================
# SOURCE PARAMETERS
# =============================================
f0    = 15.0        # Dominant source frequency [Hz]
t0    = 1.5 / f0    # Source time shift [s]

# =============================================
# GRID PARAMETERS
# =============================================
Lx    = 2000.0      # Domain width [m]
Lz    = 2000.0      # Domain depth [m]

# Spatial sampling for 4th order: dx < vs_min / (5 * f0)
# Using minimum S-wave velocity (softest layer = 600 m/s)
vs_min   = 600.0
dx_max   = vs_min / (5 * f0)   # = 8.0 m
dx       = 5.0                  # Use 5m for safety margin
dz       = 5.0

nx    = int(Lx / dx)    # = 400
nz    = int(Lz / dz)    # = 400

# =============================================
# TIME PARAMETERS
# =============================================
# CFL for elastic: dt <= dx / (vp_max * sqrt(2))
# Use maximum P-wave velocity for CFL
vp_max   = 5000.0   # maximum vp in any layer [m/s]
dt_max   = dx / (vp_max * np.sqrt(2))
dt       = 0.0005   # [s] -- conservative, well below CFL limit
T_total  = 1.0      # [s]
nt       = int(T_total / dt)   # = 2000 steps

# Courant number
courant  = vp_max * dt / dx

def print_parameters():
    print("=" * 60)
    print("  CEE490 GRADUATE -- ELASTIC WAVE SIMULATION")
    print("=" * 60)
    print(f"  Grid:              {nx} x {nz} = {nx*nz:,} points")
    print(f"  Domain:            {Lx}m x {Lz}m")
    print(f"  Grid spacing:      dx={dx}m, dz={dz}m")
    print(f"  P-wave speed:      vp={vp} m/s")
    print(f"  S-wave speed:      vs={vs:.1f} m/s")
    print(f"  Lame mu:           {mu:.3e} Pa")
    print(f"  Lame lambda:       {lam:.3e} Pa")
    print(f"  Time step:         dt={dt}s")
    print(f"  Total time:        {T_total}s ({nt} steps)")
    print(f"  Source frequency:  f0={f0} Hz")
    print(f"  CFL check:         {courant:.4f} <= {1/np.sqrt(2):.4f} "
          f"=> {'PASS' if courant <= 1/np.sqrt(2) else 'FAIL'}")
    pts_per_wl = vs_min / (f0 * dx)
    print(f"  Points/wavelength: {pts_per_wl:.1f} (need >= 5 for 4th order)")
    print(f"  vp/vs ratio:       {vp/vs:.2f} (Poisson's ratio = "
          f"{(vp**2 - 2*vs**2)/(2*(vp**2-vs**2)):.3f})")
    print("=" * 60)
```

---

### Step 1.3 — The Staggered Grid

This is the most important concept of the project. In the elastic wave equation,
different physical quantities naturally live at different spatial locations.
This is the Virieux (1986) staggered grid.

```
      ix:   0      1      2      3
            |      |      |      |
iz: 0  --  Txx,Tzz  Vx   Txx,Tzz  Vx  ...
            |      |      |      |
    1  --   Vz    Txz     Vz    Txz ...
            |      |      |      |
    2  --  Txx,Tzz  Vx   Txx,Tzz  Vx  ...

Vx  lives at (i+1/2, j)    -- half step right
Vz  lives at (i, j+1/2)    -- half step down
Txx,Tzz at (i, j)          -- integer grid points
Txz lives at (i+1/2, j+1/2)-- half step right AND down
```

Why does this work? Because the FD stencil for dVx/dx needs Vx values on both sides
of a stress point. By staggering, every derivative is automatically a central difference
over half a grid spacing, giving 2nd order accuracy naturally.

Create `elastic_grid.py`:

```python
import numpy as np
from parameters import nx, nz

def initialize_elastic_fields():
    """
    Initialize all five elastic wavefield arrays to zero.

    Field layout on the staggered grid:
      vx     : particle velocity in x  [shape: nz, nx]
      vz     : particle velocity in z  [shape: nz, nx]
      txx    : normal stress in xx     [shape: nz, nx]
      tzz    : normal stress in zz     [shape: nz, nx]
      txz    : shear stress xz         [shape: nz, nx]

    All arrays have the same shape (nz, nx) but represent values
    at different physical sub-grid locations (staggering is implicit
    in how the FD operators are applied).
    """
    vx  = np.zeros((nz, nx))
    vz  = np.zeros((nz, nx))
    txx = np.zeros((nz, nx))
    tzz = np.zeros((nz, nx))
    txz = np.zeros((nz, nx))
    return vx, vz, txx, tzz, txz


def initialize_acoustic_fields():
    """Three time levels for acoustic leapfrog."""
    p_past    = np.zeros((nz, nx))
    p_present = np.zeros((nz, nx))
    p_future  = np.zeros((nz, nx))
    return p_past, p_present, p_future
```

---

### Step 1.4 — Full Elastic FD Stencils

Create `elastic_fd.py`:

```python
import numpy as np

def elastic_step_2nd(vx, vz, txx, tzz, txz,
                     rho, lam, mu, dt, dx, dz):
    """
    One time step of the 2nd-order elastic wave equation
    on the Virieux (1986) staggered grid.

    The update is split into two half-steps:
      Step A: Update velocities from stresses (equation of motion)
      Step B: Update stresses from velocities (Hooke's law)

    This velocity-stress formulation is the standard in seismology.

    Parameters:
      vx, vz         : velocity components   [nz, nx]
      txx, tzz, txz  : stress components     [nz, nx]
      rho, lam, mu   : medium parameters     [nz, nx] (heterogeneous)
      dt, dx, dz     : time step and spacings

    Returns:
      Updated vx, vz, txx, tzz, txz
    """
    nz, nx = vx.shape

    # --- STEP A: Update particle velocities from stresses ---
    # Equation of motion:
    #   rho * dv_x/dt = dtau_xx/dx + dtau_xz/dz
    #   rho * dv_z/dt = dtau_xz/dx + dtau_zz/dz

    iz = np.arange(1, nz-1)
    ix = np.arange(1, nx-1)
    IZ, IX = np.meshgrid(iz, ix, indexing='ij')

    # dTxx/dx -- central difference, txx staggered relative to vx
    dtxx_dx = (txx[IZ, IX] - txx[IZ, IX-1]) / dx
    # dTxz/dz -- central difference
    dtxz_dz = (txz[IZ, IX] - txz[IZ-1, IX]) / dz
    # dTxz/dx -- central difference
    dtxz_dx = (txz[IZ, IX+1] - txz[IZ, IX]) / dx
    # dTzz/dz -- central difference
    dtzz_dz = (tzz[IZ+1, IX] - tzz[IZ, IX]) / dz

    vx[IZ, IX] += dt / rho[IZ, IX] * (dtxx_dx + dtxz_dz)
    vz[IZ, IX] += dt / rho[IZ, IX] * (dtxz_dx + dtzz_dz)

    # --- STEP B: Update stresses from velocities ---
    # Hooke's law:
    #   dtau_xx/dt = (lam+2mu)*dv_x/dx + lam*dv_z/dz
    #   dtau_zz/dt = lam*dv_x/dx + (lam+2mu)*dv_z/dz
    #   dtau_xz/dt = mu*(dv_x/dz + dv_z/dx)

    dvx_dx = (vx[IZ, IX+1] - vx[IZ, IX]) / dx
    dvz_dz = (vz[IZ+1, IX] - vz[IZ, IX]) / dz
    dvx_dz = (vx[IZ+1, IX] - vx[IZ, IX]) / dz
    dvz_dx = (vz[IZ, IX] - vz[IZ, IX-1]) / dx

    c11 = lam[IZ, IX] + 2 * mu[IZ, IX]   # P-wave modulus
    c12 = lam[IZ, IX]                     # off-diagonal Lame

    txx[IZ, IX] += dt * (c11 * dvx_dx + c12 * dvz_dz)
    tzz[IZ, IX] += dt * (c12 * dvx_dx + c11 * dvz_dz)
    txz[IZ, IX] += dt * mu[IZ, IX] * (dvx_dz + dvz_dx)

    return vx, vz, txx, tzz, txz


def elastic_step_4th(vx, vz, txx, tzz, txz,
                     rho, lam, mu, dt, dx, dz):
    """
    One time step of the 4th-order elastic wave equation.

    Uses 4th-order central difference stencil for all spatial derivatives.
    Coefficients: c1 = 9/8, c2 = -1/24

    This is the Levander (1988) scheme -- the standard in modern seismology.
    Requires 5-6 points per wavelength vs 10-12 for 2nd order.
    """
    nz, nx = vx.shape
    c1 =  9.0 / 8.0    # 4th order coefficient 1
    c2 = -1.0 / 24.0   # 4th order coefficient 2

    iz = np.arange(2, nz-2)
    ix = np.arange(2, nx-2)
    IZ, IX = np.meshgrid(iz, ix, indexing='ij')

    # 4th order derivative helper
    def d_dx_4th(f, IZ, IX):
        return (c1*(f[IZ, IX+1]-f[IZ, IX  ]) +
                c2*(f[IZ, IX+2]-f[IZ, IX-1])) / dx

    def d_dz_4th(f, IZ, IX):
        return (c1*(f[IZ+1, IX]-f[IZ,   IX]) +
                c2*(f[IZ+2, IX]-f[IZ-1, IX])) / dz

    # --- STEP A: Update velocities ---
    dtxx_dx = d_dx_4th(txx, IZ, IX)
    dtxz_dz = d_dz_4th(txz, IZ, IX)
    dtxz_dx = d_dx_4th(txz, IZ, IX)
    dtzz_dz = d_dz_4th(tzz, IZ, IX)

    vx[IZ, IX] += dt / rho[IZ, IX] * (dtxx_dx + dtxz_dz)
    vz[IZ, IX] += dt / rho[IZ, IX] * (dtxz_dx + dtzz_dz)

    # --- STEP B: Update stresses ---
    dvx_dx = d_dx_4th(vx, IZ, IX)
    dvz_dz = d_dz_4th(vz, IZ, IX)
    dvx_dz = d_dz_4th(vx, IZ, IX)
    dvz_dx = d_dx_4th(vz, IZ, IX)

    c11 = lam[IZ, IX] + 2 * mu[IZ, IX]
    c12 = lam[IZ, IX]

    txx[IZ, IX] += dt * (c11 * dvx_dx + c12 * dvz_dz)
    tzz[IZ, IX] += dt * (c12 * dvx_dx + c11 * dvz_dz)
    txz[IZ, IX] += dt * mu[IZ, IX] * (dvx_dz + dvz_dx)

    return vx, vz, txx, tzz, txz
```

---

### Step 1.5 — Acoustic FD for Comparison

Create `acoustic_fd.py`:

```python
import numpy as np

def laplacian_2nd(p, dx, dz):
    nz, nx = p.shape
    lap = np.zeros_like(p)
    iz = np.arange(1, nz-1)
    ix = np.arange(1, nx-1)
    IZ, IX = np.meshgrid(iz, ix, indexing='ij')
    lap[IZ, IX] = (
        (p[IZ, IX+1] - 2*p[IZ, IX] + p[IZ, IX-1]) / dx**2 +
        (p[IZ+1, IX] - 2*p[IZ, IX] + p[IZ-1, IX]) / dz**2
    )
    return lap


def laplacian_4th(p, dx, dz):
    nz, nx = p.shape
    lap = np.zeros_like(p)
    iz = np.arange(2, nz-2)
    ix = np.arange(2, nx-2)
    IZ, IX = np.meshgrid(iz, ix, indexing='ij')
    d2p_dx2 = (
        -p[IZ, IX+2] + 16*p[IZ, IX+1] - 30*p[IZ, IX]
        + 16*p[IZ, IX-1] - p[IZ, IX-2]
    ) / (12 * dx**2)
    d2p_dz2 = (
        -p[IZ+2, IX] + 16*p[IZ+1, IX] - 30*p[IZ, IX]
        + 16*p[IZ-1, IX] - p[IZ-2, IX]
    ) / (12 * dz**2)
    lap[IZ, IX] = d2p_dx2 + d2p_dz2
    return lap


def leapfrog_step(p_past, p_present, c, dt, dx, dz, order=4):
    """
    Leapfrog acoustic update.
    p_future = 2*p - p_past + (c*dt)^2 * Lap(p)
    """
    if order == 4:
        lap = laplacian_4th(p_present, dx, dz)
    else:
        lap = laplacian_2nd(p_present, dx, dz)
    return 2.0*p_present - p_past + (c*dt)**2 * lap
```

**Day 1 Checklist:**
- [ ] CFL check PASS with vp_max = 5000 m/s
- [ ] Points per wavelength >= 5 using vs_min = 600 m/s
- [ ] vx, vz, txx, tzz, txz all initialize to zeros, shape (400, 400)
- [ ] elastic_step_2nd runs one step without index errors
- [ ] elastic_step_4th runs one step without index errors
- [ ] print_parameters() shows all values clearly

---

## Day 2 — Sources, Boundaries, Velocity Models

### Step 2.1 — Ricker Wavelet and Moment Tensor Source

Create `source.py`:

```python
import numpy as np
import matplotlib.pyplot as plt

def ricker_wavelet(t, f0, t0=None):
    """
    Ricker wavelet: second derivative of a Gaussian.
    Standard seismic source. Peak frequency = f0.
    """
    if t0 is None:
        t0 = 1.5 / f0
    tau = t - t0
    a   = np.pi * f0 * tau
    return (1.0 - 2.0*a**2) * np.exp(-a**2)


def inject_explosive_source(txx, tzz, amplitude, iz_src, ix_src):
    """
    Explosive (isotropic) source -- injects into both normal stress components.
    Simulates a pressure source or small explosion.
    Produces both P and S waves.
    """
    txx[iz_src, ix_src] += amplitude
    tzz[iz_src, ix_src] += amplitude
    return txx, tzz


def inject_moment_tensor(txx, tzz, txz, amplitude, iz_src, ix_src,
                         Mxx=1.0, Mzz=1.0, Mxz=0.0):
    """
    Moment tensor source -- simulates a realistic earthquake focal mechanism.

    Parameters:
      Mxx, Mzz, Mxz : moment tensor components
        Mxx=Mzz=1, Mxz=0  --> explosion (isotropic)
        Mxx=-Mzz=1, Mxz=0 --> vertical strike-slip fault
        Mxx=0, Mxz=1       --> 45-degree thrust fault
    """
    txx[iz_src, ix_src] += Mxx * amplitude
    tzz[iz_src, ix_src] += Mzz * amplitude
    txz[iz_src, ix_src] += Mxz * amplitude
    return txx, tzz, txz


def plot_ricker(f0, dt, nt_plot=400, filename='figures/ricker_wavelet.png'):
    from scipy.fft import fft, fftfreq
    t   = np.arange(nt_plot) * dt
    stf = np.array([ricker_wavelet(ti, f0) for ti in t])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

    ax1.plot(t*1000, stf, 'k', linewidth=1.5)
    ax1.axhline(0, color='gray', linewidth=0.5)
    ax1.set_xlabel('Time [ms]', fontsize=12)
    ax1.set_ylabel('Amplitude', fontsize=12)
    ax1.set_title(f'Ricker Wavelet  f0={f0} Hz', fontsize=12)
    ax1.grid(True, alpha=0.3)

    N     = len(stf)
    freqs = fftfreq(N, dt)[:N//2]
    power = np.abs(fft(stf))[:N//2]
    ax2.plot(freqs, power/power.max(), 'b', linewidth=1.5)
    ax2.axvline(f0, color='r', linestyle='--', label=f'f0={f0} Hz')
    ax2.set_xlabel('Frequency [Hz]', fontsize=12)
    ax2.set_ylabel('Normalized Power', fontsize=12)
    ax2.set_title('Frequency Spectrum', fontsize=12)
    ax2.set_xlim([0, 4*f0])
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.show()
    print(f"  Saved: {filename}")
```

---

### Step 2.2 — Absorbing Boundaries

Create `boundaries.py`:

```python
import numpy as np

def build_sponge_mask(nz, nx, n_sponge=50, alpha=0.010):
    """
    Cerjan sponge absorbing boundary.
    n_sponge=50 recommended for elastic simulations (larger domain).
    """
    def profile(n_total, n_sponge, alpha):
        p = np.ones(n_total)
        for i in range(n_sponge):
            damp = np.exp(-(alpha * (n_sponge - i))**2)
            p[i]           = damp
            p[n_total-1-i] = damp
        return p

    pz   = profile(nz, n_sponge, alpha)
    px   = profile(nx, n_sponge, alpha)
    return np.outer(pz, px)


def apply_sponge_elastic(vx, vz, txx, tzz, txz, mask):
    """Apply sponge to all five elastic field components."""
    return (vx*mask, vz*mask, txx*mask, tzz*mask, txz*mask)


def apply_sponge_acoustic(p_future, p_present, mask):
    return p_future*mask, p_present*mask
```

---

### Step 2.3 — All Velocity Models

Create `velocity_models.py`:

```python
import numpy as np
import matplotlib.pyplot as plt

def homogeneous_elastic(nz, nx, vp=3000.0, vs=1732.0, rho_val=2200.0):
    """
    Uniform elastic medium. Used for validation.
    Analytical P-wave arrival: t = r/vp
    Analytical S-wave arrival: t = r/vs
    """
    vp_m  = np.ones((nz, nx)) * vp
    vs_m  = np.ones((nz, nx)) * vs
    rho_m = np.ones((nz, nx)) * rho_val
    mu_m  = rho_m * vs_m**2
    lam_m = rho_m * vp_m**2 - 2*mu_m
    return vp_m, vs_m, rho_m, lam_m, mu_m


def layered_elastic(nz, nx):
    """
    Three-layer elastic model:
      Layer 1 (top 30%):    soft sediment
      Layer 2 (middle 40%): competent rock
      Layer 3 (bottom 30%): hard basement
    """
    vp_m  = np.ones((nz, nx))
    vs_m  = np.ones((nz, nx))
    rho_m = np.ones((nz, nx))

    z1 = int(0.30 * nz)
    z2 = int(0.70 * nz)

    # Layer 1: soft sediment
    vp_m[:z1,  :] = 1800.0;  vs_m[:z1,  :] = 900.0;  rho_m[:z1,  :] = 1900.0
    # Layer 2: rock
    vp_m[z1:z2,:] = 3000.0;  vs_m[z1:z2,:] = 1732.0; rho_m[z1:z2,:] = 2200.0
    # Layer 3: hard basement
    vp_m[z2:,  :] = 5000.0;  vs_m[z2:,  :] = 2887.0; rho_m[z2:,  :] = 2700.0

    mu_m  = rho_m * vs_m**2
    lam_m = rho_m * vp_m**2 - 2*mu_m
    return vp_m, vs_m, rho_m, lam_m, mu_m


def realistic_geology(nz, nx, dx, dz):
    """
    Realistic geological model for earthquake engineering.

    Features:
      1. Sedimentary basin  -- elliptical soft zone (site amplification target)
      2. Hard granite basement -- high velocity, reflects energy
      3. Fault zone         -- narrow low-velocity channel
      4. Weathered surface layer -- very soft near-surface

    This setup is designed for the site amplification study.
    The basin geometry mimics alluvial valleys in earthquake-prone regions.
    """
    # Background: competent rock
    vp_m  = np.ones((nz, nx)) * 3000.0
    vs_m  = np.ones((nz, nx)) * 1732.0
    rho_m = np.ones((nz, nx)) * 2200.0

    # Weathered near-surface layer (top 5%)
    ws = int(0.05 * nz)
    vp_m[:ws,  :] = 1500.0
    vs_m[:ws,  :] = 600.0
    rho_m[:ws, :] = 1800.0

    # Sedimentary basin (elliptical)
    basin_cx = nx // 2
    basin_cz = int(0.25 * nz)
    basin_rx = int(0.20 * nx)
    basin_rz = int(0.12 * nz)
    for iz in range(nz):
        for ix in range(nx):
            dist = ((ix-basin_cx)/basin_rx)**2 + ((iz-basin_cz)/basin_rz)**2
            if dist <= 1.0:
                vp_m[iz, ix]  = 1500.0
                vs_m[iz, ix]  = 500.0
                rho_m[iz, ix] = 1700.0

    # Hard granite basement
    bd = int(0.60 * nz)
    vp_m[bd:,  :] = 5500.0
    vs_m[bd:,  :] = 3175.0
    rho_m[bd:, :] = 2700.0

    # Fault zone (vertical low-velocity channel)
    fx   = nx // 2 + nx // 5
    fh   = 4
    ftop = int(0.10 * nz)
    fbot = bd
    vp_m[ftop:fbot, fx-fh:fx+fh]  = 1000.0
    vs_m[ftop:fbot, fx-fh:fx+fh]  = 400.0
    rho_m[ftop:fbot, fx-fh:fx+fh] = 1600.0

    mu_m  = rho_m * vs_m**2
    lam_m = rho_m * vp_m**2 - 2*mu_m

    # Clip negative lam (unphysical)
    lam_m = np.maximum(lam_m, 0.0)

    return vp_m, vs_m, rho_m, lam_m, mu_m


def basin_1d(nz, dz, vs_basin=500.0, vs_rock=1732.0,
             rho_basin=1700.0, rho_rock=2200.0, H_basin=200.0):
    """
    Simple 1D basin model for analytical site amplification comparison.

    Parameters:
      H_basin : basin depth [m]
    Returns:
      vs, rho as 1D arrays
    """
    vs  = np.ones(nz) * vs_rock
    rho = np.ones(nz) * rho_rock
    n_basin = int(H_basin / dz)
    vs[:n_basin]  = vs_basin
    rho[:n_basin] = rho_basin
    return vs, rho


def plot_velocity_model(vp_m, vs_m, dx, dz, title, filename):
    nz, nx = vp_m.shape
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    for ax, field, label in zip([ax1, ax2],
                                  [vp_m, vs_m],
                                  ['P-wave velocity [m/s]', 'S-wave velocity [m/s]']):
        im = ax.imshow(field, cmap='jet',
                       extent=[0, nx*dx, nz*dz, 0],
                       aspect='auto')
        plt.colorbar(im, ax=ax, label=label)
        ax.set_xlabel('Distance x [m]', fontsize=12)
        ax.set_ylabel('Depth z [m]', fontsize=12)

    ax1.set_title(f'{title} -- P-wave velocity', fontsize=12)
    ax2.set_title(f'{title} -- S-wave velocity', fontsize=12)
    plt.tight_layout()
    plt.savefig(f'figures/{filename}', dpi=150)
    plt.show()
    print(f"  Saved: figures/{filename}")
```

**Day 2 Checklist:**
- [ ] Ricker wavelet plot shows bell-shaped pulse peaking at f0=15 Hz
- [ ] All three velocity models build without errors
- [ ] Realistic model vp plot: weathered layer (blue), basin (cyan), rock (green), basement (red), fault (dark)
- [ ] All Lame parameters positive (no negative lam)
- [ ] Sponge mask: white interior, black edges

---

## Day 3 — Simulation Engine and First Elastic Results

### Step 3.1 — Reusable Simulation Runner

Create `simulation_runner.py`:

```python
import numpy as np
import time
from parameters import nx, nz, dx, dz, dt, nt, f0
from elastic_grid import initialize_elastic_fields, initialize_acoustic_fields
from elastic_fd import elastic_step_2nd, elastic_step_4th
from acoustic_fd import leapfrog_step
from source import ricker_wavelet, inject_explosive_source, inject_moment_tensor
from boundaries import (build_sponge_mask, apply_sponge_elastic,
                         apply_sponge_acoustic)

def run_elastic(vp_m, vs_m, rho_m, lam_m, mu_m,
                scheme='elastic_4th',
                snap_interval=100,
                receivers=None,
                source_type='explosive',
                verbose=True):
    """
    Run the 2D full elastic wave simulation.

    Parameters:
      vp_m, vs_m, rho_m, lam_m, mu_m : medium parameter arrays [nz, nx]
      scheme  : 'elastic_2nd' | 'elastic_4th'
      source_type : 'explosive' | 'strike_slip' | 'thrust'

    Returns:
      dict with snapshots, seismograms, timing information
    """
    ix_src = nx // 3
    iz_src = nz // 5

    vx, vz, txx, tzz, txz = initialize_elastic_fields()
    sponge = build_sponge_mask(nz, nx, n_sponge=50)

    if receivers is None:
        receivers = [
            {'name': 'R1_rock  (300m)', 'ix': ix_src+int(300/dx), 'iz': iz_src,
             'on': 'rock'},
            {'name': 'R2_basin (300m)', 'ix': ix_src+int(300/dx),
             'iz': int(0.20*nz), 'on': 'basin'},
            {'name': 'R3_rock  (600m)', 'ix': ix_src+int(600/dx), 'iz': iz_src,
             'on': 'rock'},
            {'name': 'R4_basin (600m)', 'ix': ix_src+int(600/dx),
             'iz': int(0.20*nz), 'on': 'basin'},
        ]
    for r in receivers:
        r['ix'] = min(r['ix'], nx-1)
        r['iz'] = min(r['iz'], nz-1)

    seismograms_vx  = {r['name']: [] for r in receivers}
    seismograms_vz  = {r['name']: [] for r in receivers}
    snapshots_vx    = []
    snapshots_vz    = []
    snap_times      = []
    time_axis       = []

    # Moment tensor components
    Mxx, Mzz, Mxz = 1.0, 1.0, 0.0
    if source_type == 'strike_slip':
        Mxx, Mzz, Mxz = 1.0, -1.0, 0.0
    elif source_type == 'thrust':
        Mxx, Mzz, Mxz = 0.0, 0.0, 1.0

    if verbose:
        print(f"\n  Scheme: {scheme} | Source: {source_type}")
        print(f"  Medium vp: {vp_m.min():.0f} - {vp_m.max():.0f} m/s")
        print(f"  Medium vs: {vs_m.min():.0f} - {vs_m.max():.0f} m/s")

    wall_start = time.time()

    for it in range(nt):
        t = it * dt

        # FD update
        if scheme == 'elastic_4th':
            vx,vz,txx,tzz,txz = elastic_step_4th(
                vx,vz,txx,tzz,txz, rho_m,lam_m,mu_m, dt,dx,dz)
        else:
            vx,vz,txx,tzz,txz = elastic_step_2nd(
                vx,vz,txx,tzz,txz, rho_m,lam_m,mu_m, dt,dx,dz)

        # Source injection
        amp = ricker_wavelet(t, f0)
        txx, tzz, txz = inject_moment_tensor(
            txx, tzz, txz, amp, iz_src, ix_src, Mxx, Mzz, Mxz)

        # Absorbing boundaries
        vx,vz,txx,tzz,txz = apply_sponge_elastic(vx,vz,txx,tzz,txz, sponge)

        # Record
        time_axis.append(t)
        for r in receivers:
            seismograms_vx[r['name']].append(vx[r['iz'], r['ix']])
            seismograms_vz[r['name']].append(vz[r['iz'], r['ix']])

        if it % snap_interval == 0:
            snapshots_vx.append(vx.copy())
            snapshots_vz.append(vz.copy())
            snap_times.append(t)

        if verbose and it % 500 == 0:
            maxv = max(np.max(np.abs(vx)), np.max(np.abs(vz)))
            print(f"    step {it:5d}/{nt}  t={t:.4f}s  max|v|={maxv:.3e}")
            if np.isnan(maxv) or np.isinf(maxv):
                print("  ERROR: Diverged! Check CFL and medium parameters.")
                break

    wall_time = time.time() - wall_start

    if verbose:
        print(f"  Completed in {wall_time:.1f}s")

    return {
        'snapshots_vx': snapshots_vx,
        'snapshots_vz': snapshots_vz,
        'snap_times':   snap_times,
        'seismograms_vx': seismograms_vx,
        'seismograms_vz': seismograms_vz,
        'time_axis':    time_axis,
        'ix_src': ix_src, 'iz_src': iz_src,
        'scheme': scheme,
        'wall_time': wall_time,
        'receivers': receivers,
    }


def run_acoustic(c_model, scheme='leapfrog_4th',
                 snap_interval=100, verbose=True):
    """Run acoustic simulation for method comparison."""
    from parameters import rho as rho_val
    ix_src = nx // 3
    iz_src = nz // 5
    p_past, p_present, _ = initialize_acoustic_fields()
    sponge   = build_sponge_mask(nz, nx, n_sponge=50)
    snapshots = []
    snap_times = []
    time_axis  = []
    seismo     = []
    ix_rec = min(ix_src + int(300/dx), nx-1)
    iz_rec = iz_src

    order = 4 if scheme == 'leapfrog_4th' else 2
    wall_start = time.time()

    for it in range(nt):
        t        = it * dt
        p_future = leapfrog_step(p_past, p_present, c_model, dt, dx, dz, order)
        amp      = ricker_wavelet(t, f0)
        p_future[iz_src, ix_src] += (dt**2 / rho_val) * amp
        p_future, p_present = apply_sponge_acoustic(p_future, p_present, sponge)
        p_past    = p_present.copy()
        p_present = p_future.copy()
        time_axis.append(t)
        seismo.append(p_present[iz_rec, ix_rec])
        if it % snap_interval == 0:
            snapshots.append(p_present.copy())
            snap_times.append(t)

    return {
        'snapshots': snapshots, 'snap_times': snap_times,
        'seismograms': {'R1': seismo},
        'time_axis': time_axis,
        'ix_src': ix_src, 'iz_src': iz_src,
        'wall_time': time.time() - wall_start,
    }
```

### Step 3.2 — Run First Elastic Simulation

Add to `main.py`:

```python
from parameters import print_parameters, nx, nz, dx, dz, dt, f0, vp, vs
from velocity_models import (homogeneous_elastic, layered_elastic,
                               realistic_geology, plot_velocity_model)
from source import plot_ricker
from simulation_runner import run_elastic

print_parameters()
plot_ricker(f0, dt)

# Build models
vp_h, vs_h, rho_h, lam_h, mu_h = homogeneous_elastic(nz, nx)
vp_r, vs_r, rho_r, lam_r, mu_r = realistic_geology(nz, nx, dx, dz)

plot_velocity_model(vp_h, vs_h, dx, dz, 'Homogeneous', 'model_homo.png')
plot_velocity_model(vp_r, vs_r, dx, dz, 'Realistic Geology', 'model_realistic.png')

# First elastic run -- homogeneous, explosive source
print("\n[Running elastic simulation...]")
res_homo = run_elastic(vp_h, vs_h, rho_h, lam_h, mu_h,
                        scheme='elastic_4th',
                        source_type='explosive')
```

**Day 3 Checklist:**
- [ ] Elastic simulation runs all nt steps without NaN
- [ ] Snapshots show two distinct wavefronts (P and S waves)
- [ ] P wavefront is clearly ahead of S wavefront at same time
- [ ] Ratio of wavefront radii ≈ vp/vs = sqrt(3) ≈ 1.73
- [ ] Seismogram shows two arrivals at each receiver

---

## Day 4 — Convergence, Dispersion, Site Amplification

### Step 4.1 — Convergence Analysis

Create `analysis.py`:

```python
import numpy as np
import matplotlib.pyplot as plt

def measure_arrival_time(seismogram, time_axis, threshold=0.01):
    max_amp = max(abs(v) for v in seismogram)
    thr     = threshold * max_amp
    for i, val in enumerate(seismogram):
        if abs(val) > thr:
            return time_axis[i]
    return None


def convergence_analysis(vp_true, vs_true, f0, T_total=0.3):
    """
    Convergence study for both acoustic 2nd order and 4th order.
    Measures P-wave arrival time error vs grid spacing.
    Plots log-log with measured slopes.

    Expected: slope=2 for 2nd order, slope=4 for 4th order.
    """
    from parameters import dt
    from acoustic_fd import leapfrog_step
    from source import ricker_wavelet
    from boundaries import build_sponge_mask, apply_sponge_acoustic

    grid_spacings = [20.0, 10.0, 5.0, 2.5]
    errors_2nd    = []
    errors_4th    = []
    receiver_dist = 300.0
    rho_val       = 2200.0

    print("\n=== CONVERGENCE ANALYSIS ===")
    print(f"{'dx [m]':>8} {'Error 2nd':>12} {'Error 4th':>12}")
    print("-" * 36)

    for dx_val in grid_spacings:
        nx_c   = int(800 / dx_val)
        nz_c   = int(800 / dx_val)
        nt_c   = int(T_total / dt)
        c_homo = np.ones((nz_c, nx_c)) * vp_true
        ix_src = nx_c // 2
        iz_src = nz_c // 4
        ix_rec = min(ix_src + int(receiver_dist/dx_val), nx_c-3)
        t_analy = receiver_dist / vp_true

        for order in [2, 4]:
            p_past    = np.zeros((nz_c, nx_c))
            p_present = np.zeros((nz_c, nx_c))
            sponge    = build_sponge_mask(nz_c, nx_c, n_sponge=30)
            seismo    = []
            time_ax   = []

            for it in range(nt_c):
                t        = it * dt
                p_future = leapfrog_step(p_past, p_present, c_homo,
                                          dt, dx_val, dx_val, order)
                amp      = ricker_wavelet(t, f0)
                p_future[iz_src, ix_src] += (dt**2/rho_val) * amp
                p_future, p_present = apply_sponge_acoustic(
                    p_future, p_present, sponge)
                p_past    = p_present.copy()
                p_present = p_future.copy()
                seismo.append(p_present[iz_src, ix_rec])
                time_ax.append(t)

            t_num = measure_arrival_time(seismo, time_ax)
            err   = abs(t_num - t_analy)/t_analy if t_num else 1.0
            if order == 2:
                errors_2nd.append(err)
            else:
                errors_4th.append(err)

        print(f"{dx_val:>8.1f} {errors_2nd[-1]:>12.6f} {errors_4th[-1]:>12.6f}")

    # Log-log plot
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.loglog(grid_spacings, errors_2nd, 'b-o', lw=2, ms=8,
               label='Leapfrog 2nd order')
    ax.loglog(grid_spacings, errors_4th, 'r-s', lw=2, ms=8,
               label='Leapfrog 4th order')

    ref_x = np.array([grid_spacings[0], grid_spacings[-1]])
    ax.loglog(ref_x, errors_2nd[0]*(ref_x/grid_spacings[0])**2,
               'b--', alpha=0.5, label='Slope=2 (theory)')
    ax.loglog(ref_x, errors_4th[0]*(ref_x/grid_spacings[0])**4,
               'r--', alpha=0.5, label='Slope=4 (theory)')

    s2 = np.polyfit(np.log(grid_spacings), np.log(errors_2nd), 1)[0]
    s4 = np.polyfit(np.log(grid_spacings), np.log(errors_4th), 1)[0]

    ax.text(0.05, 0.15,
            f'Measured slope 2nd: {s2:.2f} (expect 2.0)\n'
            f'Measured slope 4th: {s4:.2f} (expect 4.0)',
            transform=ax.transAxes, fontsize=11,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    ax.set_xlabel('Grid spacing dx [m]', fontsize=13)
    ax.set_ylabel('Relative error', fontsize=13)
    ax.set_title('Convergence Analysis -- Error vs Grid Spacing', fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, which='both', alpha=0.3)
    plt.tight_layout()
    plt.savefig('figures/convergence_analysis.png', dpi=150)
    plt.show()
    print(f"\n  Slopes: 2nd={s2:.2f}, 4th={s4:.2f}")
    return grid_spacings, errors_2nd, errors_4th


def dispersion_analysis(vp, vs, dx, dt, f0,
                         filename='figures/dispersion_analysis.png'):
    """
    Analytical grid dispersion analysis.

    Computes numerical phase velocity as a function of sampling ratio dx/lambda
    for 2nd and 4th order schemes.

    Formula (2nd order):
      c_num/c_true = (lambda/(pi*dx)) * arcsin(pi*dx/lambda * c*dt/dx)

    This explains WHY 4th order needs fewer points per wavelength:
    it maintains c_num/c_true closer to 1 over a wider frequency range.
    """
    # Sampling ratios from 0.01 to 0.5 (0.5 = Nyquist)
    sampling = np.linspace(0.01, 0.5, 500)   # dx/lambda
    courant  = vp * dt / dx

    # 2nd order phase velocity ratio
    arg_2nd  = np.clip(courant * np.pi * sampling, -1, 1)
    c_ratio_2nd = np.arcsin(arg_2nd) / (np.pi * sampling * courant)

    # 4th order phase velocity ratio (more complex stencil)
    # Approximate using known result for 4th order
    arg_4th  = np.clip(courant * (np.sin(np.pi*sampling) +
                (1/6)*np.sin(2*np.pi*sampling)), -1, 1)
    c_ratio_4th = (1/(np.pi*sampling)) * np.arcsin(np.clip(arg_4th, -1, 1))

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(sampling, c_ratio_2nd, 'b-', lw=2, label='2nd order')
    ax.plot(sampling, c_ratio_4th, 'r-', lw=2, label='4th order')
    ax.axhline(1.0, color='k', linestyle='--', lw=1, label='Exact')
    ax.axhline(0.99, color='gray', linestyle=':', lw=1, label='1% error')
    ax.axvline(0.1, color='b', linestyle=':', alpha=0.5, label='10 pts/wl (2nd)')
    ax.axvline(0.2, color='r', linestyle=':', alpha=0.5, label='5 pts/wl (4th)')
    ax.set_xlabel('Spatial sampling ratio dx/lambda', fontsize=13)
    ax.set_ylabel('c_numerical / c_true', fontsize=13)
    ax.set_title('Grid Dispersion Analysis\n'
                 '4th order maintains accuracy to larger dx/lambda',
                 fontsize=13)
    ax.set_ylim([0.8, 1.05])
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.show()
    print(f"  Saved: {filename}")
```

---

### Step 4.2 — Site Amplification Analysis

This is the graduate-level engineering application that connects your simulation to a real civil engineering problem.

Create `site_amplification.py`:

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq

def compute_amplification_factor(seismo_basin, seismo_rock, time_axis, dt):
    """
    Compute frequency-dependent amplification factor.

    Amplification = |FFT(basin_motion)| / |FFT(rock_motion)|

    This is the spectral ratio method -- standard in earthquake engineering.
    An amplification factor > 1 at frequency f means waves at that frequency
    are amplified by the basin relative to rock.

    In the 1985 Mexico City earthquake, amplification factors of 8 were
    recorded at the natural frequency of the lake bed sediments.
    """
    N      = len(seismo_basin)
    freqs  = fftfreq(N, dt)[:N//2]

    spec_basin = np.abs(fft(np.array(seismo_basin)))[:N//2]
    spec_rock  = np.abs(fft(np.array(seismo_rock)))[:N//2]

    # Avoid division by zero
    amp_factor = spec_basin / (spec_rock + 1e-10)

    return freqs, amp_factor


def analytical_resonant_frequencies(vs_basin, H_basin, n_modes=4):
    """
    Analytical resonant frequencies of a 1D soil layer over bedrock.

    Formula: f_n = (2n-1) * vs / (4*H)
    n = 1: fundamental mode
    n = 2: first overtone
    etc.

    This is the 1D quarter-wavelength resonance condition.
    Your numerical simulation should show peaks near these frequencies.
    If it does, your simulation is validated against analytical theory.
    """
    frequencies = []
    for n in range(1, n_modes+1):
        fn = (2*n - 1) * vs_basin / (4 * H_basin)
        frequencies.append(fn)
        print(f"  Mode {n}: f = {fn:.2f} Hz  (T = {1/fn:.3f} s)")
    return frequencies


def plot_site_amplification(results_realistic, results_homo,
                             dx, dz, f0, dt,
                             vs_basin=500.0, H_basin=200.0,
                             filename='figures/site_amplification.png'):
    """
    Full site amplification analysis figure.

    Panel 1: Seismograms at basin vs rock receivers (same distance)
    Panel 2: Frequency spectra at both locations
    Panel 3: Amplification factor with analytical resonance marked

    This directly answers the engineering question:
    By how much does the basin amplify ground shaking, and at what frequencies?
    """
    # Get seismograms from realistic model
    # R1 = rock receiver, R2 = basin receiver (same horizontal distance)
    seismo_rock  = results_realistic['seismograms_vz'].get(
        list(results_realistic['seismograms_vz'].keys())[0], [])
    seismo_basin = results_realistic['seismograms_vz'].get(
        list(results_realistic['seismograms_vz'].keys())[1], [])
    time_axis    = results_realistic['time_axis']
    t_ms         = np.array(time_axis) * 1000

    freqs, amp_factor = compute_amplification_factor(
        seismo_basin, seismo_rock, time_axis, dt)

    # Analytical resonant frequencies
    print("\n=== ANALYTICAL RESONANT FREQUENCIES ===")
    print(f"  Basin: vs={vs_basin} m/s, H={H_basin} m")
    res_freqs = analytical_resonant_frequencies(vs_basin, H_basin, n_modes=5)

    fig, axes = plt.subplots(3, 1, figsize=(12, 14))

    # Panel 1: Seismograms
    ax = axes[0]
    trace_r = np.array(seismo_rock)
    trace_b = np.array(seismo_basin)
    norm    = max(np.max(np.abs(trace_r)), np.max(np.abs(trace_b))) + 1e-10
    ax.plot(t_ms, trace_r/norm, 'b-', lw=1.0, label='Rock site')
    ax.plot(t_ms, trace_b/norm, 'r-', lw=1.0, label='Basin site', alpha=0.8)
    ax.set_xlabel('Time [ms]', fontsize=12)
    ax.set_ylabel('Normalized velocity', fontsize=12)
    ax.set_title('Seismograms: Basin vs Rock (same distance from source)',
                 fontsize=12)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    # Panel 2: Frequency spectra
    ax = axes[1]
    N  = len(seismo_rock)
    fr = fftfreq(N, dt)[:N//2]
    sr = np.abs(fft(np.array(seismo_rock)))[:N//2]
    sb = np.abs(fft(np.array(seismo_basin)))[:N//2]
    ax.plot(fr, sr/sr.max(), 'b-', lw=1.5, label='Rock spectrum')
    ax.plot(fr, sb/sb.max(), 'r-', lw=1.5, label='Basin spectrum')
    ax.axvline(f0, color='k', linestyle='--', lw=1, label=f'Source f0={f0}Hz')
    ax.set_xlabel('Frequency [Hz]', fontsize=12)
    ax.set_ylabel('Normalized power', fontsize=12)
    ax.set_title('Frequency Spectra', fontsize=12)
    ax.set_xlim([0, 4*f0])
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    # Panel 3: Amplification factor
    ax = axes[2]
    # Smooth the amplification factor
    from scipy.signal import savgol_filter
    amp_smooth = savgol_filter(amp_factor, window_length=21, polyorder=3)
    ax.plot(freqs, amp_smooth, 'purple', lw=2, label='Amplification factor')
    ax.axhline(1.0, color='gray', linestyle='--', lw=1, label='No amplification')

    colors = ['r', 'orange', 'gold', 'green', 'cyan']
    for i, (fn, col) in enumerate(zip(res_freqs, colors)):
        ax.axvline(fn, color=col, linestyle='--', lw=1.5,
                   label=f'Mode {i+1}: {fn:.1f} Hz')

    ax.set_xlabel('Frequency [Hz]', fontsize=12)
    ax.set_ylabel('Amplification factor  |basin|/|rock|', fontsize=12)
    ax.set_title('Site Amplification Factor\n'
                 'Dashed lines = analytical resonant frequencies '
                 f'(vs={vs_basin}m/s, H={H_basin}m)',
                 fontsize=12)
    ax.set_xlim([0, 4*f0])
    ax.set_ylim([0, None])
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.show()
    print(f"  Saved: {filename}")

    # Print peak amplification frequencies
    peak_idx = np.argmax(amp_smooth[freqs < 3*f0])
    peak_freq = freqs[peak_idx]
    peak_amp  = amp_smooth[peak_idx]
    print(f"\n  Peak amplification: {peak_amp:.1f}x at {peak_freq:.2f} Hz")
    print(f"  Fundamental resonance (analytical): {res_freqs[0]:.2f} Hz")
    print(f"  Agreement: {abs(peak_freq - res_freqs[0])/res_freqs[0]*100:.1f}%")
```

---

### Step 4.3 — Performance Analysis

Add to `analysis.py`:

```python
def performance_analysis(vp_model, f0, dt, dx, dz, nx, nz, nt):
    """
    Compare computational cost vs accuracy for 2nd and 4th order schemes.

    For a given accuracy target (e.g. 1% error):
      - 2nd order needs dx = 2.0m  (many grid points, high cost)
      - 4th order needs dx = 5.0m  (fewer points, lower cost)

    Shows that 4th order is cheaper overall despite more work per point.
    This is the key insight for computational efficiency.
    """
    import time
    from acoustic_fd import leapfrog_step
    from source import ricker_wavelet
    from boundaries import build_sponge_mask, apply_sponge_acoustic

    rho_val = 2200.0
    results = {}

    print("\n=== PERFORMANCE ANALYSIS ===")
    print(f"{'Scheme':>15} {'Grid':>12} {'Time [s]':>10} "
          f"{'Memory [MB]':>12} {'Error [%]':>10}")
    print("-" * 63)

    for order, dx_val in [(2, 2.5), (4, 5.0)]:
        nx_p   = int(800 / dx_val)
        nz_p   = int(800 / dx_val)
        nt_p   = int(0.3 / dt)
        c_p    = np.ones((nz_p, nx_p)) * vp_model
        ix_src = nx_p // 2
        iz_src = nz_p // 4
        ix_rec = min(ix_src + int(300/dx_val), nx_p-3)
        sponge = build_sponge_mask(nz_p, nx_p, n_sponge=20)
        p_past = np.zeros((nz_p, nx_p))
        p_pres = np.zeros((nz_p, nx_p))
        seismo = []
        time_ax = []

        t_start = time.time()

        for it in range(nt_p):
            t_it     = it * dt
            p_fut    = leapfrog_step(p_past, p_pres, c_p,
                                      dt, dx_val, dx_val, order)
            amp      = ricker_wavelet(t_it, f0)
            p_fut[iz_src, ix_src] += (dt**2/rho_val) * amp
            p_fut, p_pres = apply_sponge_acoustic(p_fut, p_pres, sponge)
            p_past = p_pres.copy()
            p_pres = p_fut.copy()
            seismo.append(p_pres[iz_src, ix_rec])
            time_ax.append(t_it)

        elapsed = time.time() - t_start
        memory  = 5 * nx_p * nz_p * 8 / 1024**2   # 5 arrays, 8 bytes/float

        t_analy = 300.0 / vp_model
        t_num   = measure_arrival_time(seismo, time_ax)
        err_pct = abs(t_num - t_analy)/t_analy*100 if t_num else 999

        key = f'order_{order}'
        results[key] = {'time': elapsed, 'memory': memory, 'error': err_pct,
                        'nx': nx_p, 'nz': nz_p}

        print(f"  Order {order} (dx={dx_val}m): {nx_p}x{nz_p}  "
              f"{elapsed:>8.2f}s  {memory:>10.1f}MB  {err_pct:>8.3f}%")

    r2 = results['order_2']
    r4 = results['order_4']
    print(f"\n  4th order is {r2['time']/r4['time']:.1f}x faster for same accuracy")
    print(f"  4th order uses {r2['memory']/r4['memory']:.1f}x less memory")
    return results
```

**Day 4 Checklist:**
- [ ] Convergence slopes: 2nd order ≈ 2.0, 4th order ≈ 4.0
- [ ] Dispersion plot shows 4th order stays closer to c_num/c_true = 1
- [ ] Site amplification plot shows peak near analytical resonant frequency
- [ ] Numerical peak frequency within 15% of analytical (1D theory is approximate)
- [ ] Performance analysis shows 4th order faster for same accuracy
- [ ] P and S wave arrivals visible as separate peaks in seismograms

---

## Day 5 — Professional Visualization, Report & Submission

### Step 5.1 — Professional Visualization

Create `visualize.py`:

```python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def plot_elastic_comparison(results, nx, nz, dx, dz,
                             snap_idx, filename):
    """
    Two-panel figure: vx wavefield and vz wavefield side by side.
    Shows that P and S waves have different particle motion directions.
    """
    snap_vx = results['snapshots_vx'][snap_idx]
    snap_vz = results['snapshots_vz'][snap_idx]
    t       = results['snap_times'][snap_idx]
    ix_src  = results['ix_src']
    iz_src  = results['iz_src']

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

    for ax, snap, label in zip([ax1, ax2],
                                 [snap_vx, snap_vz],
                                 ['Horizontal velocity v_x', 'Vertical velocity v_z']):
        vmax = np.max(np.abs(snap)) + 1e-10
        im   = ax.imshow(snap, cmap='RdBu_r',
                          vmin=-vmax, vmax=vmax,
                          extent=[0, nx*dx, nz*dz, 0],
                          aspect='auto')
        ax.plot(ix_src*dx, iz_src*dz, 'y*', markersize=14, zorder=5)
        ax.set_title(f'{label}\nt = {t*1000:.1f} ms', fontsize=12)
        ax.set_xlabel('Distance x [m]', fontsize=11)
        ax.set_ylabel('Depth z [m]', fontsize=11)
        plt.colorbar(im, ax=ax, label='Velocity [m/s]', shrink=0.8)

    plt.suptitle('Full Elastic Wavefield -- P and S waves', fontsize=14, y=1.01)
    plt.tight_layout()
    plt.savefig(f'figures/{filename}', dpi=150, bbox_inches='tight')
    plt.show()
    print(f"  Saved: figures/{filename}")


def plot_wavefield_geology_overlay(results, velocity_model_vp,
                                    velocity_model_vs,
                                    nx, nz, dx, dz,
                                    snap_idx, title, filename):
    """
    Wavefield with P and S velocity contours overlaid.
    Professional seismology figure style.
    """
    snap   = results['snapshots_vz'][snap_idx]
    t      = results['snap_times'][snap_idx]
    ix_src = results['ix_src']
    iz_src = results['iz_src']
    vmax   = np.max(np.abs(snap)) + 1e-10

    X, Z   = np.arange(nx)*dx, np.arange(nz)*dz
    XX, ZZ = np.meshgrid(X, Z)

    fig, ax = plt.subplots(figsize=(12, 9))
    im = ax.imshow(snap, cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                   extent=[0, nx*dx, nz*dz, 0], aspect='auto', alpha=0.88)

    # P-wave velocity contours
    cp = ax.contour(XX, ZZ, velocity_model_vp,
                    levels=[1500, 3000, 5000, 5500],
                    colors='k', linewidths=0.8,
                    linestyles=['-','-','-','-'])
    ax.clabel(cp, fmt='Vp=%d m/s', fontsize=8)

    # S-wave velocity contours
    cs = ax.contour(XX, ZZ, velocity_model_vs,
                    levels=[500, 1732, 3000],
                    colors='white', linewidths=0.6,
                    linestyles=['--','--','--'], alpha=0.7)
    ax.clabel(cs, fmt='Vs=%d m/s', fontsize=8, colors='white')

    ax.plot(ix_src*dx, iz_src*dz, 'y*', markersize=16,
            label='Source', zorder=5)
    plt.colorbar(im, ax=ax, label='Vertical velocity v_z [m/s]')
    ax.set_xlabel('Distance x [m]', fontsize=13)
    ax.set_ylabel('Depth z [m]', fontsize=13)
    ax.set_title(f'{title}\nt = {t*1000:.1f} ms', fontsize=13)
    ax.legend(fontsize=11)
    plt.tight_layout()
    plt.savefig(f'figures/{filename}', dpi=150)
    plt.show()
    print(f"  Saved: figures/{filename}")


def plot_seismograms_ps_separation(results, dx, vp, vs,
                                    filename='figures/seismograms_elastic.png'):
    """
    Plot seismograms showing both P and S wave arrivals.
    Mark theoretical arrival times for both wave types.
    This is a direct validation of the elastic simulation.
    """
    receivers  = results['receivers']
    time_axis  = results['time_axis']
    t_ms       = np.array(time_axis) * 1000

    n = len(receivers)
    fig, axes = plt.subplots(n, 1, figsize=(13, 3.5*n), sharex=True)
    if n == 1:
        axes = [axes]

    ix_src = results['ix_src']
    iz_src = results['iz_src']

    print("\n=== P AND S WAVE ARRIVAL VALIDATION ===")
    print(f"{'Receiver':>20} {'Dist':>6} {'tP theory':>10} "
          f"{'tS theory':>10} {'Ratio vp/vs':>12}")
    print("-" * 62)

    for ax, r in zip(axes, receivers):
        name  = r['name']
        dist  = abs(r['ix'] - ix_src) * dx

        trace_vz = np.array(results['seismograms_vz'][name])
        trace_vx = np.array(results['seismograms_vx'][name])

        ax.plot(t_ms, trace_vz/max(np.max(np.abs(trace_vz)),1e-10),
                'k-', lw=0.9, label='v_z (vertical)')
        ax.plot(t_ms, trace_vx/max(np.max(np.abs(trace_vx)),1e-10),
                'b-', lw=0.9, alpha=0.7, label='v_x (horizontal)')
        ax.axhline(0, color='gray', lw=0.4)
        ax.set_ylabel(name[:20], fontsize=9)
        ax.grid(True, alpha=0.3)

        tp = dist / vp * 1000
        ts = dist / vs * 1000
        ax.axvline(tp, color='r', linestyle='--', lw=1.5,
                   label=f'P arrival: {tp:.1f}ms')
        ax.axvline(ts, color='orange', linestyle='--', lw=1.5,
                   label=f'S arrival: {ts:.1f}ms')
        ax.legend(fontsize=8, loc='upper right')

        print(f"{name[:20]:>20} {dist:>6.0f} {tp:>10.2f} "
              f"{ts:>10.2f} {tp/ts*vs/vp:>12.3f}")

    axes[-1].set_xlabel('Time [ms]', fontsize=12)
    axes[0].set_title('Elastic Seismograms -- P and S Wave Arrivals\n'
                      '(red=P wave, orange=S wave)', fontsize=12)
    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.show()
    print(f"\n  Saved: {filename}")


def create_elastic_animation(results, velocity_model_vp,
                              nx, nz, dx, dz,
                              component='vz',
                              filename='figures/elastic_wave.gif'):
    """Animated elastic wavefield with geology overlay."""
    snaps     = results[f'snapshots_{component}']
    times     = results['snap_times']
    ix_src    = results['ix_src']
    iz_src    = results['iz_src']

    X, Z      = np.arange(nx)*dx, np.arange(nz)*dz
    XX, ZZ    = np.meshgrid(X, Z)
    vmax      = max(np.max(np.abs(s)) for s in snaps) + 1e-10

    fig, ax = plt.subplots(figsize=(9, 8))
    im  = ax.imshow(snaps[0], cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                    extent=[0, nx*dx, nz*dz, 0], aspect='auto', animated=True)
    ax.contour(XX, ZZ, velocity_model_vp,
               levels=[1500, 3000, 5000],
               colors='k', linewidths=0.5, alpha=0.5)
    ax.plot(ix_src*dx, iz_src*dz, 'y*', markersize=12, zorder=5)
    ttl = ax.set_title('t = 0.0 ms', fontsize=12)
    plt.colorbar(im, ax=ax, label=f'v_{component} [m/s]')
    ax.set_xlabel('Distance x [m]', fontsize=12)
    ax.set_ylabel('Depth z [m]', fontsize=12)

    def update(frame):
        im.set_array(snaps[frame])
        ttl.set_text(f't = {times[frame]*1000:.1f} ms')
        return im, ttl

    ani = animation.FuncAnimation(fig, update, frames=len(snaps),
                                  interval=80, blit=True)
    ani.save(filename, writer='pillow', fps=12, dpi=100)
    plt.close()
    print(f"  Saved: {filename}")
```

---

### Step 5.2 — Complete main.py

```python
# main.py -- Run everything: python main.py
import numpy as np

from parameters import (print_parameters, nx, nz, dx, dz, dt, f0,
                         vp, vs, rho, lam, mu)
from velocity_models import (homogeneous_elastic, layered_elastic,
                               realistic_geology, plot_velocity_model)
from source import plot_ricker
from boundaries import build_sponge_mask
from simulation_runner import run_elastic, run_acoustic
from analysis import (convergence_analysis, dispersion_analysis,
                       performance_analysis)
from site_amplification import (plot_site_amplification,
                                  analytical_resonant_frequencies)
from visualize import (plot_elastic_comparison,
                        plot_wavefield_geology_overlay,
                        plot_seismograms_ps_separation,
                        create_elastic_animation)

print("\n" + "="*60)
print("  CEE490 GRADUATE -- ELASTIC WAVE PROJECT")
print("="*60)

# =============================================
# 1. SETUP
# =============================================
print_parameters()
plot_ricker(f0, dt)

# Build all models
vp_h, vs_h, rho_h, lam_h, mu_h = homogeneous_elastic(nz, nx)
vp_l, vs_l, rho_l, lam_l, mu_l = layered_elastic(nz, nx)
vp_r, vs_r, rho_r, lam_r, mu_r = realistic_geology(nz, nx, dx, dz)

plot_velocity_model(vp_h, vs_h, dx, dz, 'Homogeneous',          'model_homo.png')
plot_velocity_model(vp_l, vs_l, dx, dz, 'Three-Layer Model',    'model_layered.png')
plot_velocity_model(vp_r, vs_r, dx, dz, 'Realistic Geology',    'model_realistic.png')

# =============================================
# 2. HOMOGENEOUS ELASTIC -- VALIDATION
# =============================================
print("\n[1/6] Running homogeneous elastic (validation)...")
res_homo = run_elastic(vp_h, vs_h, rho_h, lam_h, mu_h,
                        scheme='elastic_4th', source_type='explosive')

plot_elastic_comparison(res_homo, nx, nz, dx, dz,
                         snap_idx=5, filename='elastic_homo_wavefield.png')
plot_seismograms_ps_separation(res_homo, dx, vp=vp, vs=vs)

# =============================================
# 3. REALISTIC GEOLOGY -- MAIN RESULT
# =============================================
print("\n[2/6] Running realistic geological model...")
res_real = run_elastic(vp_r, vs_r, rho_r, lam_r, mu_r,
                        scheme='elastic_4th', source_type='explosive')

plot_wavefield_geology_overlay(res_real, vp_r, vs_r, nx, nz, dx, dz,
                                snap_idx=6,
                                title='Full Elastic Wavefield -- Realistic Geology',
                                filename='wavefield_geology.png')
create_elastic_animation(res_real, vp_r, nx, nz, dx, dz,
                          component='vz')

# =============================================
# 4. SITE AMPLIFICATION ANALYSIS
# =============================================
print("\n[3/6] Site amplification analysis...")
print("\n  Analytical resonant frequencies:")
analytical_resonant_frequencies(vs_basin=500.0, H_basin=200.0, n_modes=5)

plot_site_amplification(res_real, res_homo,
                         dx=dx, dz=dz, f0=f0, dt=dt,
                         vs_basin=500.0, H_basin=200.0)

# =============================================
# 5. ANALYSIS
# =============================================
print("\n[4/6] Convergence analysis...")
convergence_analysis(vp_true=vp, vs_true=vs, f0=f0)

print("\n[5/6] Dispersion analysis...")
dispersion_analysis(vp=vp, vs=vs, dx=dx, dt=dt, f0=f0)

print("\n[6/6] Performance analysis...")
performance_analysis(vp_model=vp, f0=f0, dt=dt,
                      dx=dx, dz=dz, nx=nx, nz=nz, nt=nt)

print("\n" + "="*60)
print("  ALL DONE -- check figures/ folder")
print("="*60)
```

---

### Step 5.3 — Report Structure

```
TITLE: Full Elastic Wave Propagation in Heterogeneous Media:
       A Comparative Study of Finite-Difference Schemes with
       Application to Seismic Site Amplification

ABSTRACT (200 words)
  State the problem, methods, key findings.
  Include: schemes compared, convergence slopes achieved,
  peak amplification factor found, agreement with analytical theory.

1. INTRODUCTION (1.5 pages)
   - Seismic wave propagation and earthquake engineering importance
   - Why elastic (not acoustic): P waves, S waves, mode conversion
   - The Virieux (1986) staggered grid -- brief history and citation count
   - What this project contributes

2. GOVERNING EQUATIONS (2 pages)
   - Full elastic wave equation: all 5 coupled PDEs
   - P-wave and S-wave speeds derived from Lame parameters
   - Velocity-stress formulation
   - Staggered grid layout with diagram
   - 2nd and 4th order FD stencils with coefficients
   - CFL stability condition derivation

3. IMPLEMENTATION (1.5 pages)
   - Grid parameters table
   - Staggered grid initialization
   - Source: Ricker wavelet + moment tensor formulation
   - Sponge absorbing boundaries
   - All velocity models described with geological justification

4. VALIDATION (2 pages)
   - P and S wave speed verification (arrival times vs analytical)
   - Convergence analysis: log-log plot, measured slopes
   - Dispersion analysis: numerical vs physical phase velocity
   - Performance analysis: cost vs accuracy table

5. RESULTS: ELASTIC WAVE PROPAGATION (2 pages)
   - Homogeneous model: vx and vz snapshots showing P/S separation
   - Realistic geology: wavefield with contour overlay
   - Three-layer model: mode conversion at interfaces
   - Different source mechanisms (explosive vs strike-slip)

6. APPLICATION: SEISMIC SITE AMPLIFICATION (2 pages)
   - Motivation: Mexico City 1985, Northridge 1994
   - Method: spectral ratio technique
   - Results: amplification vs frequency plot
   - Comparison to analytical 1D resonance theory
   - Engineering implications: which frequencies are dangerous

7. DISCUSSION (1.5 pages)
   - 2nd vs 4th order: quantitative comparison
   - Why 4th order is more efficient (not just more accurate)
   - Limitations: 2D vs 3D, acoustic vs full elastic, sponge imperfections
   - How attenuation would change results (future work)

8. CONCLUSION (0.5 page)

REFERENCES
   Moczo et al. (2004)
   Virieux (1986) -- P-SV wave propagation, staggered grid
   Levander (1988) -- 4th order staggered grid schemes
   Cerjan et al. (1985) -- absorbing boundaries
   Boore (1972) -- early FD seismology
```

---

### Step 5.4 — Final Submission Checklist

```
FIGURES (15 total)
[ ] model_homo.png              -- Vp and Vs homogeneous
[ ] model_layered.png           -- Vp and Vs three-layer
[ ] model_realistic.png         -- Vp and Vs realistic geology
[ ] ricker_wavelet.png          -- Source time function + spectrum
[ ] elastic_homo_wavefield.png  -- vx and vz panels, homogeneous
[ ] wavefield_geology.png       -- vz with geology contours overlaid
[ ] seismograms_elastic.png     -- P and S arrivals marked
[ ] site_amplification.png      -- 3-panel amplification figure
[ ] convergence_analysis.png    -- log-log slopes 2 and 4
[ ] dispersion_analysis.png     -- c_num/c_true vs dx/lambda
[ ] elastic_wave.gif            -- animated elastic wavefield

CODE
[ ] python main.py runs end-to-end
[ ] CFL check PASS at startup
[ ] Negative lam clipped to zero in realistic model
[ ] All functions have docstrings

PHYSICS VALIDATION
[ ] P-wave arrival time within 3% of analytical (r/vp)
[ ] S-wave arrival time within 3% of analytical (r/vs)
[ ] vp/vs ratio = sqrt(3) = 1.732 confirmed from seismograms
[ ] Convergence slope 2nd: ~2.0, 4th: ~4.0
[ ] Site amplification peak near analytical resonance frequency

GRADUATE LEVEL ELEMENTS
[ ] Full elastic (not acoustic) -- 5 coupled PDEs
[ ] Staggered grid (Virieux 1986)
[ ] Moment tensor source
[ ] Site amplification with spectral ratio method
[ ] Comparison to 1D analytical resonance theory
[ ] Dispersion analysis (quantitative error vs frequency)
[ ] Performance analysis (cost vs accuracy)
[ ] Report structured as journal paper with Abstract
```

---

## Quick Reference -- Key Equations

| Name | Equation | Notes |
|------|----------|-------|
| P-wave speed | vp = sqrt((lam+2mu)/rho) | Compressional |
| S-wave speed | vs = sqrt(mu/rho) | Shear, slower |
| Poisson ratio | nu = lam/(2*(lam+mu)) | Typical rock: 0.25 |
| EOM (x) | rho*dvx/dt = dtxx/dx + dtxz/dz | Equation of motion |
| EOM (z) | rho*dvz/dt = dtxz/dx + dtzz/dz | Equation of motion |
| Hooke xx | dtxx/dt = (lam+2mu)*dvx/dx + lam*dvz/dz | |
| Hooke zz | dtzz/dt = lam*dvx/dx + (lam+2mu)*dvz/dz | |
| Hooke xz | dtxz/dt = mu*(dvx/dz + dvz/dx) | Shear |
| CFL elastic | dt <= dx/(vp_max*sqrt(2)) | Use vp_max |
| Spatial 4th | dx < vs_min/(5*f_max) | Use vs (slower) |
| Resonance | f_n = (2n-1)*vs/(4*H) | Site amplification |
| Amplification | A(f) = FFT(basin)/FFT(rock) | Spectral ratio |

---

## Troubleshooting

| Problem | Cause | Fix |
|---------|-------|-----|
| Simulation diverges | CFL violated | Use vp_MAX (not average) for CFL |
| Negative lam values | Unphysical Poisson ratio | Add lam = np.maximum(lam, 0) |
| Only one wavefront visible | Source too weak for S waves | Use explosive source (Mxx=Mzz=1) |
| P and S not separated | Domain too small | Increase Lx, Lz or reduce T_total |
| Amplification factor noisy | Short time series | Increase T_total for better frequency resolution |
| No amplification peak | Basin too thin | Increase H_basin or reduce vs_basin |
| Convergence slope wrong | Stencil index error | Test on known 1D solution first |
| S waves missing in vz | Strike-slip source | Switch to explosive source for isotropic radiation |

---

## References

1. **Moczo, P., Kristek, J., Halada, L. (2004)**. *The Finite-Difference Method for Seismologists: An Introduction*. Comenius University, Bratislava. SPICE Network.

2. **Virieux, J. (1986)**. P-SV wave propagation in heterogeneous media: Velocity-stress finite-difference method. *Geophysics*, 51(4), 889–901. *(>2000 citations — the foundational staggered grid paper)*

3. **Levander, A.R. (1988)**. Fourth-order finite-difference P-SV seismograms. *Geophysics*, 53(11), 1425–1436. *(Introduced 4th order staggered grid)*

4. **Cerjan, C., Kosloff, D., Kosloff, R., Reshef, M. (1985)**. A nonreflecting boundary condition for discrete acoustic and elastic wave equations. *Geophysics*, 50(4), 705–708.

5. **Boore, D.M. (1972)**. Finite difference methods for seismic wave propagation in heterogeneous materials. *Methods in Computational Physics*, 11, 1–37.

6. **Aki, K. and Richards, P.G. (2002)**. *Quantitative Seismology*, 2nd Edition. University Science Books. *(Standard seismology reference)*
