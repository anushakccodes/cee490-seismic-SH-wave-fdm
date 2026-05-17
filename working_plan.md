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

## Notation

| Symbol | Meaning |
|---|---|
| $z$ | depth coordinate |
| $t$ | time |
| $v(z,t)$ | horizontal particle velocity |
| $\tau(z,t)$ | shear stress |
| $\rho(z)$ | density |
| $\mu(z)$ | shear modulus |
| $V_s(z)$ | shear-wave velocity |
| $\Delta z$ | spatial grid spacing |
| $\Delta t$ | time step |
| $H$ | soft-layer thickness |
| $Z$ | shear impedance |

---

# 1. Governing mathematical model

## 1.1 Physical idealization

The main model represents vertically propagating, horizontally polarized shear motion through a soil column. The dependent variables are $v(z,t)$ and $\tau(z,t)$.

The material model is defined by $\rho(z)$, $\mu(z)$, and $V_s(z)$.

## 1.2 Shear-wave velocity

The shear-wave velocity is computed from density and shear modulus:

$$
V_s(z)=\sqrt{\frac{\mu(z)}{\rho(z)}}
$$

Equivalently, the shear modulus is:

$$
\mu(z)=\rho(z)V_s^2(z)
$$

This relation links the physical material profile to wave speed and numerical stability.

## 1.3 1D velocity-stress equations

The 1D SH-wave velocity-stress system is:

$$
\rho(z)\frac{\partial v}{\partial t}=\frac{\partial \tau}{\partial z}
$$

$$
\frac{\partial \tau}{\partial t}=\mu(z)\frac{\partial v}{\partial z}
$$

Interpretation:

- stress gradient drives particle acceleration;
- velocity gradient updates shear stress;
- the coupling between $v$ and $\tau$ propagates shear waves.

## 1.4 Equivalent homogeneous wave equation

For constant $\rho$ and $\mu$, the velocity-stress system reduces to:

$$
\frac{\partial^2 v}{\partial t^2}=V_s^2\frac{\partial^2 v}{\partial z^2}
$$

This form is used for travel-time verification in the homogeneous model.

---

# 2. Numerical discretization

## 2.1 Grid definition

The vertical grid and time levels are:

$$
z_i=i\Delta z
$$

$$
t^n=n\Delta t
$$

where $i$ is the depth index and $n$ is the time index.

## 2.2 Staggered grid

The staggered layout is:

$$
v_i \quad \text{at integer grid points}
$$

$$
\tau_{i+1/2} \quad \text{halfway between velocity points}
$$

Time is also staggered:

$$
v \text{ at } n+\frac{1}{2}, \qquad \tau \text{ at } n
$$

This gives the natural leapfrog velocity-stress update.

## 2.3 Core 2nd-order finite-difference update equations

Velocity update:

$$
v_i^{n+1/2}=v_i^{n-1/2}
+\frac{\Delta t}{\rho_i}
\left(\frac{\tau_{i+1/2}^{n}-\tau_{i-1/2}^{n}}{\Delta z}\right)
$$

Stress update:

$$
\tau_{i+1/2}^{n+1}=\tau_{i+1/2}^{n}
+\mu_{i+1/2}\Delta t
\left(\frac{v_{i+1}^{n+1/2}-v_i^{n+1/2}}{\Delta z}\right)
$$

These two equations are the required solver.

## 2.4 Truncation behavior

For a central finite-difference derivative:

$$
\frac{\partial f}{\partial z}\approx\frac{f_{i+1}-f_{i-1}}{2\Delta z}
$$

The leading spatial truncation error is:

$$
\mathcal{O}(\Delta z^2)
$$

The numerical error should decrease as the grid is refined.

## 2.5 Optional 4th-order derivative

A 4th-order central derivative may be added as an extension:

$$
\frac{\partial f}{\partial z}\approx
\frac{-f_{i+2}+8f_{i+1}-8f_{i-1}+f_{i-2}}{12\Delta z}
$$

The leading spatial truncation error is:

$$
\mathcal{O}(\Delta z^4)
$$

This should only be used after the 2nd-order solver is verified.

---

# 3. Stability and grid-resolution logic

## 3.1 Courant number

The Courant number is:

$$
C=\frac{V_{s,\max}\Delta t}{\Delta z}
$$

For the explicit 1D wave solver, the stability requirement is:

$$
C\leq 1
$$

A conservative target is:

$$
C\leq 0.8
$$

## 3.2 Wavelength and points per wavelength

For frequency $f$, wavelength is:

$$
\lambda=\frac{V_s}{f}
$$

The shortest wavelength is:

$$
\lambda_{\min}=\frac{V_{s,\min}}{f_{\max}}
$$

The number of grid points per wavelength is:

$$
N_{\mathrm{ppw}}=\frac{\lambda_{\min}}{\Delta z}
$$

For the main 2nd-order scheme:

$$
N_{\mathrm{ppw}}\geq 10
$$

## 3.3 Frequency resolution

For total simulation duration $T$, the frequency resolution is:

$$
\Delta f=\frac{1}{T}
$$

The simulation duration must be long enough to resolve the expected site-response frequencies.

## 3.4 Nyquist frequency

The Nyquist frequency is:

$$
f_{\mathrm{Nyquist}}=\frac{1}{2\Delta t}
$$

The source frequency content should remain well below this limit.

---

# 4. Source function

## 4.1 Ricker wavelet

The preferred source is a Ricker wavelet:

$$
s(t)=\left[1-2\pi^2 f_0^2(t-t_0)^2\right]
\exp\left[-\pi^2 f_0^2(t-t_0)^2\right]
$$

where $f_0$ is the dominant frequency and $t_0$ is the time shift.

The source spectrum must overlap the expected site-response range.

---

# 5. Material models

## 5.1 Homogeneous reference model

The homogeneous reference model is:

$$
\rho(z)=\rho_{\mathrm{rock}}
$$

$$
V_s(z)=V_{s,\mathrm{rock}}
$$

$$
\mu(z)=\rho_{\mathrm{rock}}V_{s,\mathrm{rock}}^2
$$

The theoretical travel time is:

$$
t_{\mathrm{theory}}=\frac{d}{V_{s,\mathrm{rock}}}
$$

where $d$ is the travel distance.

Arrival-time error is:

$$
\mathrm{Error}(\%)=
\frac{\left|t_{\mathrm{num}}-t_{\mathrm{theory}}\right|}{t_{\mathrm{theory}}}\times100
$$

## 5.2 Layered soil-over-rock model

The layered model is:

$$
0\leq z\leq H \quad \text{soft soil layer}
$$

$$
z>H \quad \text{stiffer rock halfspace}
$$

Material properties are:

$$
V_s(z)=
\begin{cases}
V_{s,\mathrm{soil}}, & 0\leq z\leq H \\
V_{s,\mathrm{rock}}, & z>H
\end{cases}
$$

$$
\rho(z)=
\begin{cases}
\rho_{\mathrm{soil}}, & 0\leq z\leq H \\
\rho_{\mathrm{rock}}, & z>H
\end{cases}
$$

## 5.3 Shear impedance and reflection coefficient

Shear impedance is:

$$
Z=\rho V_s
$$

For normally incident SH waves, the reflection coefficient is approximated by:

$$
R=\frac{Z_2-Z_1}{Z_2+Z_1}
$$

where $Z_1$ is the impedance of the upper layer and $Z_2$ is the impedance of the lower layer.

A larger impedance contrast produces stronger reflection and greater potential for resonance.

---

# 6. Boundary conditions

## 6.1 Stress-free surface

At the ground surface:

$$
\tau(z=0,t)=0
$$

This represents zero shear traction at the ground surface.

## 6.2 Bottom absorbing region

A simple damping mask may be applied near the bottom boundary:

$$
D(z)=\exp\left[-\alpha\eta^2(z)\right]
$$

The fields are damped inside the absorbing zone:

$$
v\leftarrow D(z)v
$$

$$
\tau\leftarrow D(z)\tau
$$

This reduces artificial reflections from the finite computational boundary.

---

# 7. Site amplification logic

## 7.1 Surface response

The recorded surface response is:

$$
v_{\mathrm{surface}}(t)=v(z=0,t)
$$

This is recorded for the homogeneous rock reference model and the layered soil-over-rock model.

## 7.2 Fourier spectra

The Fourier-domain responses are:

$$
V_{\mathrm{ref}}(f)=\mathcal{F}\left[v_{\mathrm{surface,ref}}(t)\right]
$$

$$
V_{\mathrm{layered}}(f)=\mathcal{F}\left[v_{\mathrm{surface,layered}}(t)\right]
$$

## 7.3 Amplification factor

The site amplification factor is:

$$
A(f)=\frac{\left|V_{\mathrm{layered}}(f)\right|}{\left|V_{\mathrm{ref}}(f)\right|+\epsilon}
$$

where $\epsilon$ is a small number used only to avoid division by zero.

Interpretation:

- $A(f)=1$: no amplification relative to rock;
- $A(f)>1$: amplification at that frequency.

## 7.4 Quarter-wavelength resonance

For a soft layer of thickness $H$ over stiff rock, approximate resonant frequencies are:

$$
f_n=\frac{(2n-1)V_{s,\mathrm{soil}}}{4H}, \qquad n=1,2,3,\dots
$$

The fundamental frequency is:

$$
f_1=\frac{V_{s,\mathrm{soil}}}{4H}
$$

These frequencies should be marked on the amplification plot.

---

# 8. Numerical verification and quality checks

## 8.1 Stability check

Before each run:

$$
C=\frac{V_{s,\max}\Delta t}{\Delta z}
$$

Acceptable:

$$
C\leq1
$$

Preferred:

$$
C\leq0.8
$$

## 8.2 Arrival-time verification

For a receiver at depth $z_r$ in a homogeneous medium:

$$
t_{\mathrm{theory}}=\frac{z_r}{V_s}
$$

The percent error is:

$$
\mathrm{Error}(\%)=
\frac{|t_{\mathrm{num}}-t_{\mathrm{theory}}|}{t_{\mathrm{theory}}}\times100
$$

## 8.3 Grid-refinement behavior

Run the same homogeneous problem with decreasing grid spacing:

$$
\Delta z_1>\Delta z_2>\Delta z_3
$$

Expected trend:

$$
\mathrm{Error}\downarrow \quad \text{as} \quad \Delta z\downarrow
$$

For a clean 2nd-order scheme, the idealized relationship is:

$$
\mathrm{Error}\propto(\Delta z)^2
$$

## 8.4 Energy diagnostic

For an undamped homogeneous model, a useful diagnostic is total discrete energy:

$$
E(t)=\sum_i\left[\frac{1}{2}\rho_i v_i^2+\frac{1}{2}\frac{\tau_{i+1/2}^2}{\mu_{i+1/2}}\right]\Delta z
$$

After the source input ends, energy should remain bounded. Large artificial growth indicates instability.

---

# 9. Stage-wise implementation plan

## Stage 1: Parameters and equations

Define $\Delta z$, $\Delta t$, $T$, $H$, $\rho_{\mathrm{soil}}$, $V_{s,\mathrm{soil}}$, $\mu_{\mathrm{soil}}$, $\rho_{\mathrm{rock}}$, $V_{s,\mathrm{rock}}$, $\mu_{\mathrm{rock}}$, $f_0$, and $t_0$.

Compute:

$$
\mu=\rho V_s^2
$$

$$
Z=\rho V_s
$$

$$
C=\frac{V_{s,\max}\Delta t}{\Delta z}
$$

$$
\lambda_{\min}=\frac{V_{s,\min}}{f_{\max}}
$$

$$
N_{\mathrm{ppw}}=\frac{\lambda_{\min}}{\Delta z}
$$

$$
\Delta f=\frac{1}{T}
$$

Deliverables:

- parameter table;
- stability and grid-resolution checks;
- material-profile plot.

## Stage 2: Homogeneous solver

Implement the staggered velocity-stress update and verify with:

$$
t_{\mathrm{theory}}=\frac{z_r}{V_{s,\mathrm{rock}}}
$$

Deliverables:

- homogeneous seismogram;
- arrival-time comparison;
- wavefield space-time plot.

## Stage 3: Grid-refinement study

Run the homogeneous model at several grid spacings. For each run, compute:

$$
\mathrm{Error}(\%)=\frac{|t_{\mathrm{num}}-t_{\mathrm{theory}}|}{t_{\mathrm{theory}}}\times100
$$

Plot:

$$
\log(\mathrm{Error}) \quad \text{versus} \quad \log(\Delta z)
$$

Deliverables:

- convergence table;
- log-log grid-refinement plot;
- short interpretation of observed error trend.

## Stage 4: Layered site-response model

Build the layered model:

$$
0\leq z\leq H \quad \text{soft soil}
$$

$$
z>H \quad \text{rock}
$$

Compute:

$$
Z_{\mathrm{soil}}=\rho_{\mathrm{soil}}V_{s,\mathrm{soil}}
$$

$$
Z_{\mathrm{rock}}=\rho_{\mathrm{rock}}V_{s,\mathrm{rock}}
$$

$$
R=\frac{Z_{\mathrm{rock}}-Z_{\mathrm{soil}}}{Z_{\mathrm{rock}}+Z_{\mathrm{soil}}}
$$

$$
f_1=\frac{V_{s,\mathrm{soil}}}{4H}
$$

Deliverables:

- soil-rock material profile;
- impedance profile;
- layered wavefield plot;
- surface response comparison;
- `animation_layered_reflection.gif` showing reflection and reverberation at the soil-rock interface.

## Stage 5: Site amplification

Compute:

$$
V_{\mathrm{ref}}(f)=\mathcal{F}[v_{\mathrm{ref}}(t)]
$$

$$
V_{\mathrm{layered}}(f)=\mathcal{F}[v_{\mathrm{layered}}(t)]
$$

$$
A(f)=\frac{|V_{\mathrm{layered}}(f)|}{|V_{\mathrm{ref}}(f)|+\epsilon}
$$

$$
f_n=\frac{(2n-1)V_{s,\mathrm{soil}}}{4H}
$$

Deliverables:

- reference and layered spectra;
- amplification factor plot;
- resonance frequency table;
- engineering interpretation.

## Stage 6: Optional 4th-order comparison

Use:

$$
\frac{\partial f}{\partial z}\approx
\frac{-f_{i+2}+8f_{i+1}-8f_{i-1}+f_{i-2}}{12\Delta z}
$$

Compare arrival-time error, waveform shape, runtime, and numerical dispersion behavior.

Deliverables:

- 2nd-order versus 4th-order waveform comparison;
- optional dispersion plot;
- runtime table.

---

# 10. Optional 2D SH-wave extension

If time allows, extend the same concept to 2D SH-wave propagation in an $x$-$z$ basin cross-section.

## 10.1 2D SH equations

Unknowns:

$$
v_y(x,z,t), \qquad \tau_{xy}(x,z,t), \qquad \tau_{zy}(x,z,t)
$$

The governing equations are:

$$
\rho\frac{\partial v_y}{\partial t}
=\frac{\partial\tau_{xy}}{\partial x}
+\frac{\partial\tau_{zy}}{\partial z}
$$

$$
\frac{\partial\tau_{xy}}{\partial t}=\mu\frac{\partial v_y}{\partial x}
$$

$$
\frac{\partial\tau_{zy}}{\partial t}=\mu\frac{\partial v_y}{\partial z}
$$

## 10.2 2D stability condition

A practical 2D CFL condition is:

$$
V_{s,\max}\Delta t\sqrt{\frac{1}{\Delta x^2}+\frac{1}{\Delta z^2}}\leq1
$$

For $\Delta x=\Delta z=h$:

$$
\frac{V_{s,\max}\Delta t}{h}\leq\frac{1}{\sqrt{2}}
$$

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
   - $V_s(z)$ versus depth;
   - $Z(z)=\rho(z)V_s(z)$ versus depth;
   - soil-rock interface marked.

2. **Source plot**
   - Ricker wavelet $s(t)$;
   - source spectrum $|S(f)|$;
   - dominant frequency marked.

3. **Homogeneous verification**
   - receiver seismogram;
   - theoretical arrival time marked;
   - numerical arrival time marked.

4. **Grid-refinement plot**
   - $\log(\mathrm{Error})$ versus $\log(\Delta z)$.

5. **Layered wavefield plot**
   - depth-time velocity image;
   - surface and soil-rock interface marked.

6. **Surface motion comparison**
   - rock reference surface response;
   - layered-site surface response.

7. **Site amplification plot**
   - $|V_{\mathrm{ref}}(f)|$;
   - $|V_{\mathrm{layered}}(f)|$;
   - $A(f)$ with resonance markers.

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