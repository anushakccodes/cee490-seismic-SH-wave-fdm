# Working Plan: Finite-Difference Site Response Simulation

## Project title

Finite-Difference Simulation of Seismic Site Response in Layered Ground

---

## 1. Project purpose

This project uses finite differences to model seismic shear-wave propagation through layered ground and to estimate site amplification. The goal is to show how a classroom numerical method becomes a practical earthquake-engineering calculation.

The main deliverable is a **verified 1D velocity-stress finite-difference model** for vertically propagating shear waves. The model is first tested in a homogeneous medium, then applied to a soft soil layer over stiff rock. A simplified 2D SH-wave basin model can be added only after the 1D model works.

The project should demonstrate:

- how a PDE becomes algebraic update equations;
- why stability depends on the time step;
- why accuracy depends on wavelength resolution;
- how a soft layer changes surface motion;
- how numerical results connect to site amplification.

---

## 2. Final selected direction

> A 1D velocity-stress finite-difference simulation of vertically propagating shear waves through layered ground, verified in a homogeneous medium and applied to site amplification of a soft layer over rock.

This is the core project because it is rigorous, verifiable, and directly tied to earthquake engineering. The 2D extension is useful for visualization and basin effects, but the 1D model provides the main scientific evidence.

---

## 3. Notation

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
| $C$ | Courant number |
| $A(f)$ | amplification factor |

---

# 4. Governing mathematical model

## 4.1 Physical idealization

The main model is a vertical soil column. Motion is assumed to be horizontally polarized and vertically propagating. This is the standard simplified setting for 1D shear-wave site response.

The two unknown fields are:

- $v(z,t)$: horizontal particle velocity;
- $\tau(z,t)$: shear stress.

The material profile is defined by $\rho(z)$, $\mu(z)$, and $V_s(z)$. A soft soil layer has smaller $V_s$ and lower impedance than rock, so waves slow down, reflect, and can reverberate inside the layer.

## 4.2 Shear-wave velocity

The shear-wave velocity is computed from density and shear modulus:

$$
V_s(z)=\sqrt{\frac{\mu(z)}{\rho(z)}}
$$

Equivalently:

$$
\mu(z)=\rho(z)V_s^2(z)
$$

This relation is used to construct the material model. It also controls the stability limit because faster waves require smaller time steps.

## 4.3 1D velocity-stress equations

The 1D SH-wave velocity-stress system is:

$$
\rho(z)\frac{\partial v}{\partial t}=\frac{\partial \tau}{\partial z}
$$

$$
\frac{\partial \tau}{\partial t}=\mu(z)\frac{\partial v}{\partial z}
$$

Physical meaning:

- a stress gradient accelerates the soil particles;
- a velocity gradient creates shear deformation and updates stress;
- the repeated exchange between velocity and stress propagates the wave.

## 4.4 Equivalent homogeneous wave equation

For constant $\rho$ and $\mu$, the velocity-stress system reduces to:

$$
\frac{\partial^2 v}{\partial t^2}=V_s^2\frac{\partial^2 v}{\partial z^2}
$$

This equation is useful for verification because the wave speed is exactly $V_s$ in a homogeneous medium.

---

# 5. Numerical discretization

## 5.1 Grid definition

Depth and time are discretized as:

$$
z_i=i\Delta z
$$

$$
t^n=n\Delta t
$$

Here, $i$ indexes depth and $n$ indexes time. The numerical solution stores values only at discrete grid points and advances them step by step.

## 5.2 Why a staggered grid is used

Velocity and stress are naturally offset. Stress differences update velocity, and velocity differences update stress. A staggered grid places each variable where the derivative needs it.

The layout is:

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

This creates a leapfrog update: stress updates velocity, then velocity updates stress.

## 5.3 Core 2nd-order update equations

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

These are the main computational equations of the project.

## 5.4 Truncation error and grid refinement

A central finite-difference derivative is approximated by:

$$
\frac{\partial f}{\partial z}\approx\frac{f_{i+1}-f_{i-1}}{2\Delta z}
$$

For a 2nd-order scheme, the leading spatial truncation error is:

$$
\mathcal{O}(\Delta z^2)
$$

This means reducing $\Delta z$ should reduce the numerical error. The grid-refinement study checks this behavior.

## 5.5 Optional 4th-order derivative

A 4th-order central derivative may be added later:

$$
\frac{\partial f}{\partial z}\approx
\frac{-f_{i+2}+8f_{i+1}-8f_{i-1}+f_{i-2}}{12\Delta z}
$$

The leading spatial truncation error is:

$$
\mathcal{O}(\Delta z^4)
$$

This can reduce numerical dispersion, but it requires a wider stencil and more careful boundary treatment. It is an extension, not the core requirement.

---

# 6. Stability and grid-resolution logic

## 6.1 Courant number

The Courant number measures how far a wave travels in one time step relative to the grid spacing:

$$
C=\frac{V_{s,\max}\Delta t}{\Delta z}
$$

For the explicit 1D wave solver, the stability requirement is:

$$
C\leq 1
$$

Use a conservative target:

$$
C\leq 0.8
$$

If $C$ is too large, the numerical solution can grow artificially and become unstable.

## 6.2 Wavelength resolution

A wave of frequency $f$ has wavelength:

$$
\lambda=\frac{V_s}{f}
$$

The shortest wavelength controls the grid spacing:

$$
\lambda_{\min}=\frac{V_{s,\min}}{f_{\max}}
$$

The number of grid points per shortest wavelength is:

$$
N_{\mathrm{ppw}}=\frac{\lambda_{\min}}{\Delta z}
$$

For the main 2nd-order scheme:

$$
N_{\mathrm{ppw}}\geq 10
$$

This reduces numerical dispersion and keeps the wave speed physically reasonable.

## 6.3 Frequency resolution

The frequency spacing in the Fourier spectrum depends on total simulation duration $T$:

$$
\Delta f=\frac{1}{T}
$$

A longer simulation gives finer frequency resolution. This matters because the amplification peaks must be resolved clearly.

## 6.4 Nyquist frequency

The maximum resolvable frequency from the time step is:

$$
f_{\mathrm{Nyquist}}=\frac{1}{2\Delta t}
$$

The source frequency content should stay well below $f_{\mathrm{Nyquist}}$.

---

# 7. Source function

## 7.1 Ricker wavelet

A Ricker wavelet is used because it is smooth, finite-duration, and has controlled frequency content:

$$
s(t)=\left[1-2\pi^2 f_0^2(t-t_0)^2\right]
\exp\left[-\pi^2 f_0^2(t-t_0)^2\right]
$$

Here, $f_0$ is the dominant frequency and $t_0$ shifts the pulse so it starts near zero.

## 7.2 Source-frequency logic

The source must contain energy near the expected site-resonance frequencies. If the soil layer resonates near $f_1$, the source spectrum should include energy around $f_1$ and the first few higher modes.

The source time history and source spectrum should be plotted before interpreting amplification results.

---

# 8. Material models

## 8.1 Homogeneous reference model

The homogeneous model is used for verification only:

$$
\rho(z)=\rho_{\mathrm{rock}}
$$

$$
V_s(z)=V_{s,\mathrm{rock}}
$$

$$
\mu(z)=\rho_{\mathrm{rock}}V_{s,\mathrm{rock}}^2
$$

Because the velocity is constant, the theoretical travel time is known:

$$
t_{\mathrm{theory}}=\frac{d}{V_{s,\mathrm{rock}}}
$$

where $d$ is the travel distance.

The arrival-time error is:

$$
\mathrm{Error}(\%)=
\frac{\left|t_{\mathrm{num}}-t_{\mathrm{theory}}\right|}{t_{\mathrm{theory}}}\times100
$$

This is the simplest quantitative check that the solver propagates waves at the correct speed.

## 8.2 Layered soil-over-rock model

The engineering model is a soft layer over stiff rock:

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

This model is simple, but it captures the key mechanism of site response: impedance contrast and wave trapping.

## 8.3 Shear impedance and reflection

Shear impedance is:

$$
Z=\rho V_s
$$

At the soil-rock interface, a normally incident SH wave has approximate reflection coefficient:

$$
R=\frac{Z_2-Z_1}{Z_2+Z_1}
$$

where $Z_1$ is the impedance of the upper layer and $Z_2$ is the impedance of the lower layer.

A larger impedance contrast produces stronger reflections and greater potential for resonance in the soft layer.

---

# 9. Boundary conditions

## 9.1 Stress-free surface

At the ground surface:

$$
\tau(z=0,t)=0
$$

This represents zero shear traction at the surface. It is important because the engineering output is surface motion.

## 9.2 Bottom absorbing region

The computational domain is finite, so the bottom boundary can create artificial reflections. A damping layer near the base reduces this effect.

A simple damping mask is:

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

This is an approximate absorbing boundary. Late-time results should be interpreted carefully if boundary reflections appear.

---

# 10. Site amplification logic

## 10.1 Surface response

The main recorded output is the surface velocity:

$$
v_{\mathrm{surface}}(t)=v(z=0,t)
$$

This is recorded for both the homogeneous reference model and the layered model.

## 10.2 Fourier spectra

The Fourier-domain responses are:

$$
V_{\mathrm{ref}}(f)=\mathcal{F}\left[v_{\mathrm{surface,ref}}(t)\right]
$$

$$
V_{\mathrm{layered}}(f)=\mathcal{F}\left[v_{\mathrm{surface,layered}}(t)\right]
$$

The spectra show which frequencies are present in the surface motion.

## 10.3 Amplification factor

The site amplification factor is:

$$
A(f)=\frac{\left|V_{\mathrm{layered}}(f)\right|}{\left|V_{\mathrm{ref}}(f)\right|+\epsilon}
$$

where $\epsilon$ only prevents division by zero.

Interpretation:

- $A(f)=1$: no amplification relative to rock;
- $A(f)>1$: amplification at that frequency;
- peaks in $A(f)$: possible resonance frequencies.

## 10.4 Quarter-wavelength resonance

For a soft layer of thickness $H$ over stiff rock, approximate resonant frequencies are:

$$
f_n=\frac{(2n-1)V_{s,\mathrm{soil}}}{4H}, \qquad n=1,2,3,\dots
$$

The fundamental frequency is:

$$
f_1=\frac{V_{s,\mathrm{soil}}}{4H}
$$

These frequencies should be marked on the amplification plot. Agreement does not need to be exact; the formula is a physical interpretation tool.

---

# 11. Numerical verification and quality checks

## 11.1 Stability check

Before every run, compute:

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

## 11.2 Arrival-time verification

For a receiver at depth $z_r$ in a homogeneous medium:

$$
t_{\mathrm{theory}}=\frac{z_r}{V_s}
$$

Percent error:

$$
\mathrm{Error}(\%)=
\frac{|t_{\mathrm{num}}-t_{\mathrm{theory}}|}{t_{\mathrm{theory}}}\times100
$$

This directly checks whether the numerical wave speed is correct.

## 11.3 Grid-refinement behavior

Run the same homogeneous problem with decreasing grid spacing:

$$
\Delta z_1>\Delta z_2>\Delta z_3
$$

Expected trend:

$$
\mathrm{Error}\downarrow \quad \text{as} \quad \Delta z\downarrow
$$

For a clean 2nd-order scheme:

$$
\mathrm{Error}\propto(\Delta z)^2
$$

The observed slope may not be exact because arrival picking, source bandwidth, and finite time step affect the measurement.

## 11.4 Energy diagnostic

For an undamped homogeneous model:

$$
E(t)=\sum_i\left[\frac{1}{2}\rho_i v_i^2+\frac{1}{2}\frac{\tau_{i+1/2}^2}{\mu_{i+1/2}}\right]\Delta z
$$

After the source input ends, energy should remain bounded. Large artificial growth indicates instability.

---

# 12. Stage-wise implementation plan

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

Implement the staggered velocity-stress update and run the homogeneous model. This is the first real test of the solver.

Verify with:

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
- short interpretation of the observed trend.

## Stage 4: Layered site-response model

Build the soft-layer-over-rock model. This stage introduces the engineering mechanism: wave reflection and reverberation caused by impedance contrast.

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

Compute spectra and amplification:

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

Use the 4th-order derivative:

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

# 13. Optional 2D SH-wave extension

If time allows, extend the same concept to 2D SH-wave propagation in an $x$-$z$ basin cross-section. This is not required for the core project, but it shows how the same logic expands from a soil column to a basin.

## 13.1 2D SH equations

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

## 13.2 2D stability condition

A practical 2D CFL condition is:

$$
V_{s,\max}\Delta t\sqrt{\frac{1}{\Delta x^2}+\frac{1}{\Delta z^2}}\leq1
$$

For $\Delta x=\Delta z=h$:

$$
\frac{V_{s,\max}\Delta t}{h}\leq\frac{1}{\sqrt{2}}
$$

## 13.3 2D outputs

- wavefield snapshots;
- basin-edge scattering;
- surface amplification at multiple horizontal locations;
- 2D wavefield animation.

This extension should be presented as qualitative unless a separate validation case is implemented.

---

# 14. Visualization plan

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

# 15. Final report structure

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

# 16. Main conclusions to support

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

# 17. Limitations to state clearly

- The main model is 1D and does not capture basin-edge focusing.
- The material model is linear elastic.
- The source is idealized.
- Damping is approximate unless a formal viscoelastic model is implemented.
- The absorbing boundary is approximate.
- The quarter-wavelength resonance equation is an interpretation tool, not an exact match to every numerical result.
- Animations are explanatory visualizations, not independent validation.

---

# 18. Final project statement

This project implements a transparent finite-difference model for seismic shear-wave propagation in layered ground. It uses mathematical discretization, stability checks, verification against analytical travel time, spectral site-response analysis, and a layered-reflection animation to connect numerical methods with earthquake-engineering interpretation.