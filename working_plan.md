# Working Plan: Finite-Difference Site Response Simulation

## Project title

Finite-Difference Simulation of Seismic Site Response in Layered Ground

---

## Project purpose

This project implements a finite-difference model for seismic shear-wave propagation in layered ground and uses it to study site amplification. The goal is to connect classroom numerical methods with a practical earthquake-engineering problem.

The main deliverable is a **1D velocity-stress finite-difference model** for vertically propagating shear waves. The model is verified using a homogeneous medium and then applied to a soft soil layer over stiffer rock. A simplified 2D SH-wave basin model may be added as an optional extension after the 1D model is verified.

The project emphasizes:

- discretization of a governing PDE;
- explicit time stepping;
- stability through the Courant condition;
- grid-resolution requirements;
- arrival-time verification;
- frequency-domain site amplification;
- engineering interpretation of numerical results.

---

## Final selected direction

> A verified 1D velocity-stress finite-difference simulation of vertically propagating shear waves through layered ground, with site amplification as the main engineering application.

This direction is rigorous enough for a graduate-level course project while remaining controlled enough to validate quantitatively.

---

# 1. Governing mathematical model

## 1.1 Physical idealization

The main model represents vertically propagating, horizontally polarized shear motion through a soil column. The independent variables are:

- depth: `z`
- time: `t`

The dependent variables are:

- horizontal particle velocity: `v(z,t)`
- shear stress: `tau(z,t)`

The material properties are:

- density: `rho(z)`
- shear modulus: `mu(z)`
- shear-wave velocity: `Vs(z)`

## 1.2 Shear-wave velocity

The shear-wave velocity is computed from material properties:

```text
Vs(z) = sqrt(mu(z) / rho(z))
```

Equivalently:

```text
mu(z) = rho(z) * Vs(z)^2
```

This relation is central because both wave speed and stability depend on the material profile.

## 1.3 1D velocity-stress equations

The 1D SH-wave velocity-stress system is:

```text
rho(z) * partial v / partial t = partial tau / partial z
partial tau / partial t = mu(z) * partial v / partial z
```

Interpretation:

- stress gradient drives particle acceleration;
- velocity gradient produces shear strain rate and updates stress;
- the coupling between `v` and `tau` propagates shear waves.

## 1.4 Equivalent second-order wave equation

For a homogeneous medium, the velocity-stress system reduces to the 1D wave equation:

```text
partial^2 v / partial t^2 = Vs^2 * partial^2 v / partial z^2
```

This form is useful for understanding analytical travel time and grid-resolution requirements.

---

# 2. Numerical discretization

## 2.1 Grid definition

The vertical domain is divided into equally spaced grid points:

```text
z_i = i * Delta z
```

Time is divided into equally spaced time levels:

```text
t^n = n * Delta t
```

where:

- `Delta z` = spatial grid spacing;
- `Delta t` = time step;
- `i` = depth index;
- `n` = time index.

## 2.2 Staggered-grid layout

The solver will use a staggered grid:

```text
v_i          stored at integer grid points
tau_{i+1/2} stored halfway between velocity points
```

Time is also staggered:

```text
v is updated at half time levels: n + 1/2
tau is updated at full time levels: n
```

This gives the natural velocity-stress leapfrog structure.

## 2.3 2nd-order finite-difference update

Velocity update:

```text
v_i^{n+1/2} = v_i^{n-1/2}
              + (Delta t / rho_i)
              * (tau_{i+1/2}^n - tau_{i-1/2}^n) / Delta z
```

Stress update:

```text
tau_{i+1/2}^{n+1} = tau_{i+1/2}^n
                    + mu_{i+1/2} * Delta t
                    * (v_{i+1}^{n+1/2} - v_i^{n+1/2}) / Delta z
```

This is the required core numerical method.

## 2.4 Local truncation behavior

For the central finite-difference derivative:

```text
partial f / partial z ≈ (f_{i+1} - f_{i-1}) / (2 * Delta z)
```

The leading spatial truncation error is proportional to:

```text
O(Delta z^2)
```

The project will use this to explain why grid refinement should reduce numerical error.

## 2.5 Optional 4th-order derivative

A 4th-order central derivative can be used as an extension:

```text
partial f / partial z ≈
(-f_{i+2} + 8f_{i+1} - 8f_{i-1} + f_{i-2}) / (12 * Delta z)
```

The leading spatial truncation error is:

```text
O(Delta z^4)
```

This extension is useful for comparing numerical dispersion, but it should only be added after the 2nd-order solver is stable and verified.

---

# 3. Stability and grid-resolution logic

## 3.1 Courant number

The Courant number for the 1D model is:

```text
C = Vs_max * Delta t / Delta z
```

where `Vs_max` is the maximum shear-wave velocity in the model.

For the explicit 1D wave solver, the practical stability requirement is:

```text
C <= 1
```

A conservative target is:

```text
C <= 0.8
```

The simulation should print the Courant number before running.

## 3.2 Wavelength and points per wavelength

For a frequency `f`, the wavelength is:

```text
lambda = Vs / f
```

The shortest wavelength occurs in the slowest layer at the highest frequency of interest:

```text
lambda_min = Vs_min / f_max
```

The number of grid points per wavelength is:

```text
N_ppw = lambda_min / Delta z
```

For the main 2nd-order scheme, use:

```text
N_ppw >= 10
```

This prevents excessive numerical dispersion.

## 3.3 Frequency resolution

For a total simulation duration `T`, the frequency resolution is:

```text
Delta f = 1 / T
```

If the expected site resonance is near 1 Hz, the simulation should run long enough that `Delta f` is meaningfully smaller than 1 Hz.

## 3.4 Nyquist frequency

The maximum resolvable frequency from the time step is:

```text
f_Nyquist = 1 / (2 * Delta t)
```

The source frequency content should stay well below the Nyquist frequency.

---

# 4. Source function

## 4.1 Ricker wavelet

The preferred source is a Ricker wavelet because it is smooth and has controlled frequency content:

```text
s(t) = [1 - 2*pi^2*f0^2*(t - t0)^2]
       * exp[-pi^2*f0^2*(t - t0)^2]
```

where:

- `f0` = dominant frequency;
- `t0` = time shift so the pulse starts near zero.

## 4.2 Source-frequency selection

The source frequency must overlap the expected site-response range. If the first site resonance is expected near `f1`, choose a source frequency band that contains `f1` and the first few modes.

The source spectrum should be plotted to verify that the input motion actually excites the frequencies being analyzed.

---

# 5. Material models

## 5.1 Homogeneous reference model

The homogeneous model has constant properties:

```text
rho(z) = rho_rock
Vs(z)  = Vs_rock
mu(z)  = rho_rock * Vs_rock^2
```

This model is used for verification because the theoretical travel time is known:

```text
t_theory = travel_distance / Vs_rock
```

The numerical arrival-time error is:

```text
error_percent = |t_numerical - t_theory| / t_theory * 100
```

## 5.2 Layered soil-over-rock model

The main engineering model is:

```text
0 <= z <= H:        soft soil layer
z > H:              stiffer rock halfspace
```

Material properties:

```text
Vs(z) = Vs_soil,  rho(z) = rho_soil    for 0 <= z <= H
Vs(z) = Vs_rock,  rho(z) = rho_rock    for z > H
```

The shear impedance is:

```text
Z = rho * Vs
```

The impedance contrast controls reflection and transmission at the soil-rock interface.

## 5.3 Reflection coefficient at layer interface

For normally incident SH waves, the displacement/velocity reflection coefficient is approximated by:

```text
R = (Z_2 - Z_1) / (Z_2 + Z_1)
```

where:

- `Z_1` = impedance of the upper layer;
- `Z_2` = impedance of the lower layer.

A large impedance contrast means stronger reflection and greater potential for resonance in the soft layer.

---

# 6. Boundary conditions

## 6.1 Stress-free ground surface

At the ground surface:

```text
tau(z = 0, t) = 0
```

This represents a free surface where no shear traction is applied above the ground.

## 6.2 Bottom absorbing region

The bottom of the computational domain is artificial. To reduce reflections, a damping layer may be applied near the base.

A simple damping mask can have the form:

```text
D(z) = exp[-alpha * eta(z)^2]
```

where `eta(z)` increases from 0 at the start of the damping layer to 1 at the bottom boundary.

The fields are multiplied by the damping factor inside the absorbing zone:

```text
v <- D(z) * v
tau <- D(z) * tau
```

This boundary is approximate and should be discussed as a limitation.

---

# 7. Site amplification logic

## 7.1 Surface response

The main recorded output is surface velocity:

```text
v_surface(t) = v(z = 0, t)
```

This is recorded for:

- the homogeneous rock reference model;
- the layered soil-over-rock model.

## 7.2 Fourier amplitude spectra

Let:

```text
V_ref(f)     = FFT[v_surface, reference model]
V_layered(f) = FFT[v_surface, layered model]
```

The Fourier amplitude spectra are:

```text
|V_ref(f)|
|V_layered(f)|
```

## 7.3 Amplification factor

The site amplification factor is:

```text
A(f) = |V_layered(f)| / (|V_ref(f)| + epsilon)
```

where `epsilon` is a small number used only to avoid division by zero.

Interpretation:

- `A(f) = 1` means no amplification relative to rock;
- `A(f) > 1` means the soil layer amplifies that frequency;
- peaks in `A(f)` indicate resonant frequencies.

## 7.4 Quarter-wavelength resonance

For a soft layer of thickness `H` over stiff rock, approximate resonant frequencies are:

```text
f_n = (2n - 1) * Vs_soil / (4H)
```

where:

```text
n = 1, 2, 3, ...
```

The fundamental frequency is:

```text
f_1 = Vs_soil / (4H)
```

These frequencies should be marked on the amplification plot.

---

# 8. Numerical verification and quality checks

## 8.1 Stability check

Before every run, compute:

```text
C = Vs_max * Delta t / Delta z
```

The run is acceptable only if:

```text
C <= 1
```

Prefer:

```text
C <= 0.8
```

## 8.2 Arrival-time verification

For the homogeneous model:

```text
t_theory = z_receiver / Vs
```

Compare this with the first clear numerical arrival time:

```text
error_percent = |t_numerical - t_theory| / t_theory * 100
```

## 8.3 Grid-refinement behavior

Run the same homogeneous problem using multiple grid spacings:

```text
Delta z_1 > Delta z_2 > Delta z_3
```

Expected trend:

```text
error decreases as Delta z decreases
```

For a clean 2nd-order spatial scheme, the idealized relationship is:

```text
error ∝ (Delta z)^2
```

In the report, this should be stated as expected behavior, not guaranteed exact behavior, because arrival picking and finite source bandwidth affect the measured slope.

## 8.4 Energy check for the undamped homogeneous case

For an undamped homogeneous model, a useful diagnostic is total discrete energy:

```text
E(t) = sum[ 0.5 * rho * v^2 + 0.5 * tau^2 / mu ] * Delta z
```

Without damping and after source input ends, energy should remain approximately bounded. Large artificial growth indicates instability.

---

# 9. Stage-wise implementation plan

## Stage 1: Parameters and equations

Define:

```text
Delta z, Delta t, T, H
rho_soil, Vs_soil, mu_soil
rho_rock, Vs_rock, mu_rock
f0, t0
```

Compute and report:

```text
mu = rho * Vs^2
Z = rho * Vs
C = Vs_max * Delta t / Delta z
lambda_min = Vs_min / f_max
N_ppw = lambda_min / Delta z
Delta f = 1 / T
```

Deliverables:

- parameter table;
- stability and grid-resolution checks;
- material-profile plot.

---

## Stage 2: Homogeneous solver

Implement the staggered velocity-stress update:

```text
v_i^{n+1/2} = v_i^{n-1/2}
              + (Delta t / rho_i)
              * (tau_{i+1/2}^n - tau_{i-1/2}^n) / Delta z
```

```text
tau_{i+1/2}^{n+1} = tau_{i+1/2}^n
                    + mu_{i+1/2} * Delta t
                    * (v_{i+1}^{n+1/2} - v_i^{n+1/2}) / Delta z
```

Verify with:

```text
t_theory = z_receiver / Vs_rock
```

Deliverables:

- homogeneous seismogram;
- arrival-time comparison;
- wavefield space-time plot.

---

## Stage 3: Grid-refinement study

Run the homogeneous model at several grid spacings.

For each run, compute:

```text
error_percent = |t_numerical - t_theory| / t_theory * 100
```

Plot:

```text
log(error) versus log(Delta z)
```

Deliverables:

- convergence table;
- log-log grid-refinement plot;
- short interpretation of observed error trend.

---

## Stage 4: Layered site-response model

Build the layered profile:

```text
soft soil: 0 <= z <= H
rock:      z > H
```

Compute:

```text
Z_soil = rho_soil * Vs_soil
Z_rock = rho_rock * Vs_rock
R = (Z_rock - Z_soil) / (Z_rock + Z_soil)
f_1 = Vs_soil / (4H)
```

Deliverables:

- soil-rock material profile;
- impedance profile;
- layered wavefield plot;
- surface response comparison.

---

## Stage 5: Site amplification

Compute spectra:

```text
V_ref(f) = FFT[v_ref(t)]
V_layered(f) = FFT[v_layered(t)]
```

Compute amplification:

```text
A(f) = |V_layered(f)| / (|V_ref(f)| + epsilon)
```

Mark resonance frequencies:

```text
f_n = (2n - 1) * Vs_soil / (4H)
```

Deliverables:

- reference and layered spectra;
- amplification factor plot;
- resonance frequency table;
- engineering interpretation.

---

## Stage 6: Optional 4th-order comparison

Use the 4th-order derivative:

```text
partial f / partial z ≈
(-f_{i+2} + 8f_{i+1} - 8f_{i-1} + f_{i-2}) / (12 * Delta z)
```

Compare against the 2nd-order scheme using:

```text
arrival-time error
waveform shape
runtime
numerical dispersion behavior
```

Deliverables:

- 2nd-order versus 4th-order waveform comparison;
- optional dispersion plot;
- runtime table.

---

# 10. Optional 2D SH-wave extension

If time allows, extend the same concept to 2D SH-wave propagation in an `x-z` basin cross-section.

## 10.1 2D SH governing equations

Unknowns:

```text
v_y(x,z,t)
tau_xy(x,z,t)
tau_zy(x,z,t)
```

Equations:

```text
rho * partial v_y / partial t = partial tau_xy / partial x + partial tau_zy / partial z
partial tau_xy / partial t = mu * partial v_y / partial x
partial tau_zy / partial t = mu * partial v_y / partial z
```

## 10.2 2D stability condition

For equal grid spacing, a practical 2D CFL condition is:

```text
Vs_max * Delta t * sqrt(1 / Delta x^2 + 1 / Delta z^2) <= 1
```

For `Delta x = Delta z = h`, this becomes approximately:

```text
Vs_max * Delta t / h <= 1 / sqrt(2)
```

## 10.3 2D output quantities

The 2D extension should show:

- wavefield snapshots;
- basin-edge scattering;
- surface amplification at multiple x-locations;
- animation of wave propagation through the basin.

This extension should be presented as qualitative unless a separate validation case is implemented.

---

# 11. Visualization plan

## Required figures

1. **Material profile**
   - `Vs(z)` versus depth;
   - `Z(z) = rho(z)Vs(z)` versus depth;
   - soil-rock interface marked.

2. **Source plot**
   - Ricker wavelet `s(t)`;
   - source spectrum `|S(f)|`;
   - dominant frequency marked.

3. **Homogeneous verification**
   - receiver seismogram;
   - theoretical arrival time marked;
   - numerical arrival time marked.

4. **Grid-refinement plot**
   - `log(error)` versus `log(Delta z)`.

5. **Layered wavefield plot**
   - depth-time velocity image;
   - surface and soil-rock interface marked.

6. **Surface motion comparison**
   - rock reference surface response;
   - layered-site surface response.

7. **Site amplification plot**
   - `|V_ref(f)|`;
   - `|V_layered(f)|`;
   - `A(f)` with resonance markers.

## Optional animation

Animation is useful for communication but not validation. Recommended animation:

- velocity profile evolving with time;
- soil-rock interface marked;
- surface marked;
- time label shown;
- consistent amplitude scaling.

---

# 12. Final report structure

1. Abstract
2. Engineering motivation
3. Governing equations
4. Finite-difference discretization
5. Stability and grid-resolution criteria
6. Homogeneous-medium verification
7. Layered site-response model
8. Site amplification results
9. Optional 2D SH-wave extension
10. Limitations
11. Conclusion

---

# 13. Main conclusions to support

The final report should support these claims:

- The finite-difference method converts a seismic wave PDE into explicit algebraic update rules.
- Stability depends on the Courant number.
- Accuracy depends on grid spacing relative to the shortest wavelength.
- A homogeneous model can verify numerical wave speed through analytical travel time.
- A soft layer over rock produces frequency-dependent amplification.
- Amplification peaks can be interpreted using quarter-wavelength resonance theory.
- A 2D extension can demonstrate lateral basin effects, but the verified 1D model is the core scientific result.

---

# 14. Limitations to state clearly

- The main model is 1D and does not capture basin-edge focusing.
- The material model is linear elastic.
- The source is idealized.
- Damping is approximate unless a formal viscoelastic model is implemented.
- The absorbing boundary is approximate.
- The quarter-wavelength resonance equation is an interpretation tool, not an exact match to every numerical result.

---

# 15. Final project statement

This project implements a transparent finite-difference model for seismic shear-wave propagation in layered ground. It uses mathematical discretization, stability checks, verification against analytical travel time, and spectral site-response analysis to connect numerical methods with earthquake-engineering interpretation.