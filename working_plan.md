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

The main model represents vertically propagating, horizontally polarized shear motion through a soil column.

Independent variables:

```text
z = depth
t = time
```

Dependent variables:

```text
v(z,t)   = horizontal particle velocity
tau(z,t) = shear stress
```

Material properties:

```text
rho(z) = density
mu(z)  = shear modulus
Vs(z)  = shear-wave velocity
```

## 1.2 Shear-wave velocity

```text
Vs(z) = sqrt(mu(z) / rho(z))
```

Equivalently:

```text
mu(z) = rho(z) * Vs(z)^2
```

This relation links the material model directly to wave speed and numerical stability.

## 1.3 1D velocity-stress equations

```text
rho(z) * partial v / partial t = partial tau / partial z
partial tau / partial t = mu(z) * partial v / partial z
```

Interpretation:

- stress gradient drives particle acceleration;
- velocity gradient updates shear stress;
- the coupling between `v` and `tau` propagates shear waves.

## 1.4 Equivalent homogeneous wave equation

For constant `rho` and `mu`, the velocity-stress system reduces to:

```text
partial^2 v / partial t^2 = Vs^2 * partial^2 v / partial z^2
```

This form is used for travel-time verification in the homogeneous model.

---

# 2. Numerical discretization

## 2.1 Grid definition

```text
z_i = i * Delta z
t^n = n * Delta t
```

where:

```text
Delta z = spatial grid spacing
Delta t = time step
i       = depth index
n       = time index
```

## 2.2 Staggered grid

```text
v_i          stored at integer grid points
tau_{i+1/2} stored halfway between velocity points
```

Time is also staggered:

```text
v   at n + 1/2
tau at n
```

This gives the natural leapfrog velocity-stress update.

## 2.3 Core 2nd-order update equations

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

This is the required solver.

## 2.4 Truncation behavior

For a central finite-difference derivative:

```text
partial f / partial z ≈ (f_{i+1} - f_{i-1}) / (2 * Delta z)
```

The leading spatial truncation error is:

```text
O(Delta z^2)
```

The numerical error should decrease as the grid is refined.

## 2.5 Optional 4th-order derivative

A 4th-order central derivative may be added as an extension:

```text
partial f / partial z ≈
(-f_{i+2} + 8f_{i+1} - 8f_{i-1} + f_{i-2}) / (12 * Delta z)
```

The leading spatial truncation error is:

```text
O(Delta z^4)
```

This should only be used after the 2nd-order solver is verified.

---

# 3. Stability and grid-resolution logic

## 3.1 Courant number

```text
C = Vs_max * Delta t / Delta z
```

For the explicit 1D wave solver:

```text
C <= 1
```

Use a conservative target:

```text
C <= 0.8
```

## 3.2 Wavelength and points per wavelength

```text
lambda = Vs / f
lambda_min = Vs_min / f_max
N_ppw = lambda_min / Delta z
```

For the main 2nd-order scheme:

```text
N_ppw >= 10
```

## 3.3 Frequency resolution

```text
Delta f = 1 / T
```

The simulation duration `T` must be long enough to resolve the expected site-response frequencies.

## 3.4 Nyquist frequency

```text
f_Nyquist = 1 / (2 * Delta t)
```

The source frequency content should remain well below this limit.

---

# 4. Source function

## 4.1 Ricker wavelet

```text
s(t) = [1 - 2*pi^2*f0^2*(t - t0)^2]
       * exp[-pi^2*f0^2*(t - t0)^2]
```

where:

```text
f0 = dominant frequency
t0 = time shift
```

The source spectrum must overlap the expected site-response range.

---

# 5. Material models

## 5.1 Homogeneous reference model

```text
rho(z) = rho_rock
Vs(z)  = Vs_rock
mu(z)  = rho_rock * Vs_rock^2
```

Theoretical travel time:

```text
t_theory = travel_distance / Vs_rock
```

Arrival-time error:

```text
error_percent = |t_numerical - t_theory| / t_theory * 100
```

## 5.2 Layered soil-over-rock model

```text
0 <= z <= H: soft soil layer
z > H:       stiffer rock halfspace
```

```text
Vs(z) = Vs_soil, rho(z) = rho_soil   for 0 <= z <= H
Vs(z) = Vs_rock, rho(z) = rho_rock   for z > H
```

Shear impedance:

```text
Z = rho * Vs
```

Reflection coefficient for normally incident SH waves:

```text
R = (Z_2 - Z_1) / (Z_2 + Z_1)
```

A larger impedance contrast produces stronger reflection and greater potential for resonance.

---

# 6. Boundary conditions

## 6.1 Stress-free surface

```text
tau(z = 0, t) = 0
```

This represents zero shear traction at the ground surface.

## 6.2 Bottom absorbing region

A simple damping mask may be applied near the bottom boundary:

```text
D(z) = exp[-alpha * eta(z)^2]
```

Apply inside the damping zone:

```text
v   <- D(z) * v
tau <- D(z) * tau
```

This reduces artificial reflections from the finite computational boundary.

---

# 7. Site amplification logic

## 7.1 Surface response

```text
v_surface(t) = v(z = 0, t)
```

Record this for:

- homogeneous rock reference model;
- layered soil-over-rock model.

## 7.2 Fourier spectra

```text
V_ref(f)      = FFT[v_surface, reference]
V_layered(f)  = FFT[v_surface, layered]
```

## 7.3 Amplification factor

```text
A(f) = |V_layered(f)| / (|V_ref(f)| + epsilon)
```

Interpretation:

```text
A(f) = 1  -> no amplification
A(f) > 1  -> amplification at that frequency
```

## 7.4 Quarter-wavelength resonance

```text
f_n = (2n - 1) * Vs_soil / (4H)
```

where:

```text
n = 1, 2, 3, ...
```

Fundamental frequency:

```text
f_1 = Vs_soil / (4H)
```

These frequencies should be marked on the amplification plot.

---

# 8. Numerical verification and quality checks

## 8.1 Stability check

Before each run:

```text
C = Vs_max * Delta t / Delta z
```

Acceptable:

```text
C <= 1
```

Preferred:

```text
C <= 0.8
```

## 8.2 Arrival-time verification

```text
t_theory = z_receiver / Vs
error_percent = |t_numerical - t_theory| / t_theory * 100
```

## 8.3 Grid-refinement behavior

Run the same homogeneous problem with:

```text
Delta z_1 > Delta z_2 > Delta z_3
```

Expected trend:

```text
error decreases as Delta z decreases
```

For a clean 2nd-order scheme:

```text
error proportional to (Delta z)^2
```

## 8.4 Energy diagnostic

For an undamped homogeneous model:

```text
E(t) = sum[0.5 * rho * v^2 + 0.5 * tau^2 / mu] * Delta z
```

After the source input ends, energy should remain bounded. Large artificial growth indicates instability.

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

Compute:

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

## Stage 2: Homogeneous solver

Implement the staggered velocity-stress update and verify:

```text
t_theory = z_receiver / Vs_rock
```

Deliverables:

- homogeneous seismogram;
- arrival-time comparison;
- wavefield space-time plot.

## Stage 3: Grid-refinement study

Run the homogeneous model at several grid spacings.

For each run:

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

## Stage 4: Layered site-response model

Build:

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
- surface response comparison;
- `animation_layered_reflection.gif` showing reflection and reverberation at the soil-rock interface.

## Stage 5: Site amplification

Compute:

```text
V_ref(f) = FFT[v_ref(t)]
V_layered(f) = FFT[v_layered(t)]
A(f) = |V_layered(f)| / (|V_ref(f)| + epsilon)
f_n = (2n - 1) * Vs_soil / (4H)
```

Deliverables:

- reference and layered spectra;
- amplification factor plot;
- resonance frequency table;
- engineering interpretation.

## Stage 6: Optional 4th-order comparison

Use:

```text
partial f / partial z ≈
(-f_{i+2} + 8f_{i+1} - 8f_{i-1} + f_{i-2}) / (12 * Delta z)
```

Compare:

- arrival-time error;
- waveform shape;
- runtime;
- numerical dispersion behavior.

Deliverables:

- 2nd-order versus 4th-order waveform comparison;
- optional dispersion plot;
- runtime table.

---

# 10. Optional 2D SH-wave extension

If time allows, extend the same concept to 2D SH-wave propagation in an `x-z` basin cross-section.

## 10.1 2D SH equations

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

```text
Vs_max * Delta t * sqrt(1 / Delta x^2 + 1 / Delta z^2) <= 1
```

For `Delta x = Delta z = h`:

```text
Vs_max * Delta t / h <= 1 / sqrt(2)
```

## 10.3 2D outputs

- wavefield snapshots;
- basin-edge scattering;
- surface amplification at multiple horizontal locations;
- 2D wavefield animation.

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

## Required animation

`animation_layered_reflection.gif` should show the wave entering the soft layer, reflecting at the soil-rock interface, and reverberating between the interface and free surface.

It is useful for presentation and physical interpretation, but it should not replace the quantitative validation plots.

Recommended animation features:

- particle velocity profile evolving with time;
- soil-rock interface marked;
- free surface marked;
- reflected pulse labeled or visually distinguishable;
- consistent amplitude scaling;
- time label shown on each frame.

## Optional animation

`wave_animation.gif` may show the full 1D propagation history more generally.

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
9. Animation-based physical interpretation
10. Optional 2D SH-wave extension
11. Limitations
12. Conclusion

---

# 13. Main conclusions to support

The final report should support these claims:

- The finite-difference method converts a seismic wave PDE into explicit algebraic update rules.
- Stability depends on the Courant number.
- Accuracy depends on grid spacing relative to the shortest wavelength.
- A homogeneous model can verify numerical wave speed through analytical travel time.
- A soft layer over rock produces frequency-dependent amplification.
- Amplification peaks can be interpreted using quarter-wavelength resonance theory.
- `animation_layered_reflection.gif` helps communicate the reflection and reverberation mechanism, but the scientific evidence comes from verification and spectral analysis.
- A 2D extension can demonstrate lateral basin effects, but the verified 1D model is the core scientific result.

---

# 14. Limitations to state clearly

- The main model is 1D and does not capture basin-edge focusing.
- The material model is linear elastic.
- The source is idealized.
- Damping is approximate unless a formal viscoelastic model is implemented.
- The absorbing boundary is approximate.
- The quarter-wavelength resonance equation is an interpretation tool, not an exact match to every numerical result.
- Animations are explanatory visualizations, not independent validation.

---

# 15. Final project statement

This project implements a transparent finite-difference model for seismic shear-wave propagation in layered ground. It uses mathematical discretization, stability checks, verification against analytical travel time, spectral site-response analysis, and a layered-reflection animation to connect numerical methods with earthquake-engineering interpretation.