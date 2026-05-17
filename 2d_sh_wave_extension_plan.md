# 2D SH-Wave Extension Implementation Plan

## Purpose

This notebook will extend the existing 1D SH-wave site-response project into a 2D SH-wave finite-difference model. The goal is to simulate out-of-plane shear-wave motion in a vertical 2D cross-section, visualize wave propagation, and compare 2D behavior against the simpler 1D model.

The 2D notebook should be separate from the current 1D notebook so the original implementation remains clean and reproducible.

Suggested notebook name:

```text
notebooks/02_2d_sh_wave_extension.ipynb
```

Suggested source-code approach:

```text
src/site_response/fd_solver_2d.py
src/site_response/materials_2d.py
src/site_response/visualization_2d.py
src/site_response/analysis_2d.py
```

---

## 1. Physical Model

### 1.1 Wave type

The model should solve 2D SH-wave propagation. In SH motion, particle velocity is perpendicular to the 2D model plane.

Use:

- horizontal coordinate: `x`
- depth coordinate: `z`
- out-of-plane particle velocity: `v_y(x, z, t)`
- shear stresses: `tau_xy(x, z, t)` and `tau_zy(x, z, t)`

The medium varies in the `x-z` plane, but motion is only in the `y` direction.

---

## 2. Governing Equations

The 2D SH-wave velocity-stress system is:

```math
\rho(x,z)\frac{\partial v_y}{\partial t}
=
\frac{\partial \tau_{xy}}{\partial x}
+
\frac{\partial \tau_{zy}}{\partial z}
```

```math
\frac{\partial \tau_{xy}}{\partial t}
=
\mu(x,z)\frac{\partial v_y}{\partial x}
```

```math
\frac{\partial \tau_{zy}}{\partial t}
=
\mu(x,z)\frac{\partial v_y}{\partial z}
```

where:

```math
\mu(x,z)=\rho(x,z)V_s^2(x,z)
```

---

## 3. Main Differences From the 1D Model

| 1D model | 2D extension |
|---|---|
| one spatial coordinate `z` | two spatial coordinates `x` and `z` |
| one stress component `tau_zy` | two stress components `tau_xy`, `tau_zy` |
| vertical propagation only | lateral spreading, diffraction, focusing, basin effects |
| line plots and depth-time images | 2D snapshots, time animations, receiver arrays |
| simple layered profile | layered, basin, inclusion, or topographic models |

---

## 4. Numerical Discretization

### 4.1 Grid

Use a regular Cartesian grid:

```text
x = 0, dx, 2dx, ..., Lx
z = 0, dz, 2dz, ..., Lz
```

Recommended default values for the first implementation:

```python
model_width_m = 800.0
model_depth_m = 400.0
dx_m = 4.0
dz_m = 4.0
dt_s = 0.0015
total_time_s = 2.0
```

The first version should use `dx = dz` to simplify stability and visualization.

---

### 4.2 Staggered velocity-stress grid

Use a staggered-grid layout similar to the 1D solver:

```text
v_y      at cell centers or primary grid nodes
 tau_xy  staggered in x
 tau_zy  staggered in z
```

A practical first implementation can store arrays as:

```python
velocity_y[nz, nx]
stress_xy[nz, nx - 1]
stress_zy[nz - 1, nx]
```

Use harmonic averaging for shear modulus on stress locations:

```python
mu_x = harmonic_average(mu[:, :-1], mu[:, 1:])
mu_z = harmonic_average(mu[:-1, :], mu[1:, :])
```

---

### 4.3 Update equations

Velocity update:

```python
velocity_y[:, 1:-1] += (dt / density[:, 1:-1]) * (
    (stress_xy[:, 1:] - stress_xy[:, :-1]) / dx
    + vertical_stress_gradient
)
```

Stress updates:

```python
stress_xy += mu_x * dt * np.diff(velocity_y, axis=1) / dx
stress_zy += mu_z * dt * np.diff(velocity_y, axis=0) / dz
```

The actual indexing should be implemented carefully so array shapes match.

---

## 5. Stability Condition

For a 2D explicit finite-difference SH solver, use the approximate CFL condition:

```math
\Delta t
\leq
\frac{1}{V_{s,\max}\sqrt{\frac{1}{\Delta x^2}+\frac{1}{\Delta z^2}}}
```

For `dx = dz`, this becomes:

```math
\Delta t \leq \frac{\Delta x}{V_{s,\max}\sqrt{2}}
```

Recommended check:

```python
courant_2d = vs_max * dt * np.sqrt(1.0 / dx**2 + 1.0 / dz**2)
```

Require:

```python
courant_2d <= 1.0
```

Prefer:

```python
courant_2d <= 0.6
```

---

## 6. Material Models

Start with simple models before adding complexity.

### 6.1 Homogeneous rock model

Purpose:

- debug solver
- check circular wavefronts
- verify arrival times at receivers

Parameters:

```python
rho = 2200.0      # kg/m^3
vs = 800.0        # m/s
mu = rho * vs**2
```

Expected behavior:

- wavefront expands roughly circularly from point source
- arrival time at receiver should match distance divided by `Vs`

---

### 6.2 Horizontally layered model

Purpose:

- compare against existing 1D model
- validate vertical propagation behavior
- observe free-surface amplification and interface reflections

Example:

```text
0–50 m:    soft soil, Vs = 250 m/s, rho = 1800 kg/m^3
50–400 m:  rock,      Vs = 800 m/s, rho = 2200 kg/m^3
```

Expected behavior:

- reflections at soil-rock interface
- slower wave speed in soil layer
- surface reverberation

---

### 6.3 Optional basin model

Purpose:

- demonstrate why 2D modeling matters
- show edge-generated surface waves and focusing

Example basin:

```text
soft layer thickness varies with x
maximum thickness at center
thin soil near edges
rock below basin shape
```

Possible formula:

```python
basin_depth = base_depth + extra_depth * np.exp(-((x - x_center) / basin_width)**2)
```

Then assign soil where:

```python
z <= basin_depth[x]
```

---

## 7. Source Model

Use the same Ricker wavelet from the 1D model.

Recommended first source:

```python
source_x_m = model_width_m / 2
source_z_m = 200.0
source_type = "body_force"
```

Apply the source to the velocity field:

```python
velocity_y[source_z_index, source_x_index] += source_scale * source_values[time_index]
```

For comparison with the 1D model, later add a plane-wave-like source at depth by applying the same source over many `x` nodes.

---

## 8. Boundary Conditions

### 8.1 Free surface

At the top boundary:

```math
\tau_{zy}(x, z=0, t)=0
```

In the staggered-grid implementation, this can be handled by excluding stress above the surface and using one-sided stress gradient near the top.

---

### 8.2 Absorbing boundaries

Use damping masks near:

- left boundary
- right boundary
- bottom boundary

A simple first version can multiply velocity and stresses by smooth damping masks each time step.

Example damping mask:

```python
mask = np.exp(-strength * eta**2)
```

where `eta` goes from 0 to 1 inside the absorbing layer.

Recommended damping layer:

```python
absorbing_layer_thickness_m = 80.0
absorbing_strength = 3.5
```

Do not over-interpret late-time waves until boundary reflections are checked.

---

## 9. Receiver Layout

Use multiple receivers so the 2D model can be validated and visualized.

Recommended receiver sets:

### 9.1 Surface receiver array

```text
z = 0 m
x = 100, 200, 300, 400, 500, 600, 700 m
```

Purpose:

- surface-motion comparison
- lateral variability
- basin amplification if basin model is added

### 9.2 Internal verification receiver

For homogeneous verification:

```text
source:   (400 m, 200 m)
receiver: (500 m, 200 m)
```

Theoretical arrival:

```math
t = t_0 + \frac{\sqrt{(x_r-x_s)^2+(z_r-z_s)^2}}{V_s}
```

---

## 10. Verification Plan

### 10.1 Homogeneous arrival-time verification

Run a homogeneous model and compare numerical and theoretical arrival times.

Theoretical distance:

```math
d = \sqrt{(x_r-x_s)^2+(z_r-z_s)^2}
```

Theoretical peak arrival:

```math
t_{peak}=t_0+\frac{d}{V_s}
```

Output:

```text
receiver_x_m
receiver_z_m
distance_m
theoretical_peak_arrival_s
numerical_peak_arrival_s
absolute_error_s
percent_error
```

---

### 10.2 Wavefront shape check

For a homogeneous medium, snapshots should show nearly circular wavefronts.

Useful diagnostic:

- plot 2D velocity snapshot at selected times
- overlay theoretical radius `r = Vs * (t - t0)` centered at source

---

### 10.3 1D comparison check

For a horizontally layered model with a laterally uniform source, compare surface motion at the centerline with the 1D layered model.

This is not necessary for the first version, but it is a strong validation step.

---

## 11. Suggested Source Files

### 11.1 `materials_2d.py`

Responsibilities:

```text
- create 2D x and z axes
- create homogeneous 2D material model
- create layered 2D material model
- create optional basin model
- compute shear modulus
- compute impedance
- compute stress-grid harmonic averages
```

Suggested functions:

```python
def make_2d_axes(config_2d):
    ...

def homogeneous_2d_model(config_2d):
    ...

def layered_2d_model(config_2d):
    ...

def basin_2d_model(config_2d):
    ...

def stress_grid_mu_x(material_2d):
    ...

def stress_grid_mu_z(material_2d):
    ...
```

---

### 11.2 `fd_solver_2d.py`

Responsibilities:

```text
- define 2D simulation result dataclass
- build damping masks
- run 2D velocity-stress solver
- store snapshots
- store receiver time histories
```

Suggested dataclass:

```python
@dataclass(frozen=True)
class Simulation2DResult:
    time_s: np.ndarray
    x_m: np.ndarray
    z_m: np.ndarray
    stored_time_s: np.ndarray
    velocity_snapshots_m_s: np.ndarray
    receiver_records_m_s: dict[str, np.ndarray]
    source_time_function: np.ndarray
    source_x_m: float
    source_z_m: float
    material: Material2D
    config: Simulation2DConfig
```

---

### 11.3 `analysis_2d.py`

Responsibilities:

```text
- arrival-time verification
- receiver distance calculation
- peak picking
- surface-array comparison
- optional 1D-vs-2D comparison
```

Suggested functions:

```python
def theoretical_2d_arrival_time(source_x, source_z, receiver_x, receiver_z, vs):
    ...

def arrival_time_verification_2d(result, receiver_name):
    ...

def compute_surface_peak_amplitudes(result):
    ...
```

---

### 11.4 `visualization_2d.py`

Responsibilities:

```text
- material map plots
- velocity snapshot plots
- receiver seismogram plots
- surface-array plots
- 2D wavefield animation
```

Suggested functions:

```python
def plot_vs_map(material_2d, output_path=None):
    ...

def plot_impedance_map(material_2d, output_path=None):
    ...

def plot_velocity_snapshot(result, time_index, output_path=None):
    ...

def plot_receiver_records(result, output_path=None):
    ...

def make_2d_wave_animation(result, output_path=None):
    ...
```

---

## 12. Notebook Structure

Create:

```text
notebooks/02_2d_sh_wave_extension.ipynb
```

Recommended sections:

```text
1. Introduction and purpose
2. Governing 2D SH-wave equations
3. Numerical grid and CFL condition
4. Homogeneous 2D model setup
5. Source and receiver layout
6. Homogeneous wavefront verification
7. Arrival-time verification
8. Layered 2D model setup
9. Layered wavefield snapshots
10. Surface receiver comparison
11. Optional basin model
12. 2D animation
13. Limitations and next steps
```

---

## 13. Expected Outputs

Suggested figures:

```text
figures_2d/vs_map_homogeneous.png
figures_2d/vs_map_layered.png
figures_2d/impedance_map_layered.png
figures_2d/source_wavelet_2d.png
figures_2d/homogeneous_snapshot_t1.png
figures_2d/homogeneous_snapshot_t2.png
figures_2d/homogeneous_receiver_records.png
figures_2d/homogeneous_arrival_check.png
figures_2d/layered_snapshot_t1.png
figures_2d/layered_snapshot_t2.png
figures_2d/layered_surface_array.png
figures_2d/wavefield_2d_animation.gif
```

Suggested result files:

```text
results_2d/homogeneous_2d_results.npz
results_2d/layered_2d_results.npz
results_2d/homogeneous_2d_arrival_check.csv
results_2d/surface_receiver_summary.csv
```

---

## 14. Implementation Order

### Step 1 — Add 2D configuration

Add a new dataclass, separate from the 1D config:

```python
@dataclass(frozen=True)
class Simulation2DConfig:
    model_width_m: float = 800.0
    model_depth_m: float = 400.0
    dx_m: float = 4.0
    dz_m: float = 4.0
    time_step_s: float = 0.0015
    total_time_s: float = 2.0
    dominant_frequency_hz: float = 4.0
    soil_layer_thickness_m: float = 50.0
    soil_density_kg_m3: float = 1800.0
    soil_shear_velocity_m_s: float = 250.0
    rock_density_kg_m3: float = 2200.0
    rock_shear_velocity_m_s: float = 800.0
    absorbing_layer_thickness_m: float = 80.0
    absorbing_strength: float = 3.5
```

---

### Step 2 — Build 2D material arrays

Create arrays:

```python
density_kg_m3[nz, nx]
shear_velocity_m_s[nz, nx]
shear_modulus_pa[nz, nx]
```

Start with homogeneous and layered models only.

---

### Step 3 — Implement 2D solver

Start with:

```text
homogeneous material
single point source
one receiver
store snapshots every N steps
simple damping boundaries
```

Do not add basin geometry until the homogeneous solver behaves correctly.

---

### Step 4 — Verify arrival time

Use one receiver in the homogeneous medium and compare numerical arrival with:

```math
t_0+\frac{d}{V_s}
```

Target initial accuracy:

```text
percent error < 2–5%
```

After refinement:

```text
percent error < 1–2%
```

---

### Step 5 — Add layered model

Run a horizontally layered model and verify that:

```text
- waves slow down in the soft layer
- reflections occur at the interface
- free-surface motion differs from the homogeneous reference
```

---

### Step 6 — Add receiver array

Store multiple receiver histories along the surface.

Use these to show:

```text
- arrival-time variation
- peak-amplitude variation
- lateral uniformity in layered model
- lateral variability if basin model is added
```

---

### Step 7 — Add animation

Create a 2D wavefield animation using `imshow`.

Animation should show:

```text
- velocity field snapshots
- source location
- receiver locations
- soil-rock interface
- absorbing boundary region, optionally shaded
```

---

## 15. Validation Criteria

Minimum validation before presenting results:

```text
- 2D CFL check passes
- homogeneous wavefront looks circular
- homogeneous arrival error is acceptably small
- no obvious boundary reflection contaminates main arrival
- material interfaces appear at correct depth
- output files are saved reproducibly
```

Stronger validation:

```text
- grid refinement reduces arrival-time error
- laterally uniform layered 2D case approximately matches 1D response at centerline
- energy does not grow artificially with time
- results are stable under smaller time step
```

---

## 16. Main Risks and How to Avoid Them

### Risk 1 — Array shape mismatch

Use explicit array shapes:

```text
velocity_y: (nz, nx)
stress_xy:  (nz, nx-1)
stress_zy:  (nz-1, nx)
```

Write comments before each update showing expected shapes.

---

### Risk 2 — CFL instability

Always print the 2D CFL value before running.

If unstable:

```text
reduce dt
increase dx and dz
reduce maximum velocity
```

---

### Risk 3 — Boundary reflections

Use thick damping layers and avoid interpreting late arrivals.

For early verification, place source and receiver far from boundaries.

---

### Risk 4 — Source too strong

Start with a small source scale:

```python
source_scale = 1e-5 to 1e-4
```

Increase only after confirming stability.

---

### Risk 5 — Too much memory use

Do not store every time step.

Use:

```python
store_every = 5, 10, or 20
```

Memory estimate:

```text
velocity_snapshots size ≈ n_snapshots × nz × nx × 8 bytes
```

Example:

```text
200 snapshots × 101 × 201 × 8 ≈ 32 MB
```

This is acceptable.

---

## 17. Recommended First Milestone

The first successful 2D milestone should be:

```text
Homogeneous 2D SH-wave point-source simulation
```

Required outputs:

```text
- 2D CFL summary
- homogeneous Vs map
- source wavelet
- two or three wavefield snapshots
- receiver seismogram
- arrival-time verification table
- 2D animation GIF
```

Only after this works should the layered model and basin model be added.

---

## 18. Recommended Final Notebook Message

The final notebook should clearly state:

```text
The 2D SH-wave model extends the 1D velocity-stress formulation by adding lateral stress gradients. The homogeneous case verifies basic wave speed and wavefront behavior. The layered case demonstrates reflection, reverberation, and free-surface response. More advanced interpretation, such as basin-edge effects or true engineering transfer functions, requires stronger boundary treatment and additional validation.
```
