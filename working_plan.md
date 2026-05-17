# Working Plan: Finite-Difference Site Response Simulation

## Project title

Finite-Difference Simulation of Seismic Site Response in Layered Ground

---

## Project purpose

The project will implement a finite-difference model for seismic wave propagation in layered ground and use it to study site amplification. The goal is to demonstrate core numerical-method concepts through a practical earthquake-engineering problem.

The project is intentionally scoped around a model that is rigorous but implementable. The main deliverable will be a 1D vertically propagating shear-wave simulation using an explicit velocity-stress finite-difference scheme. The engineering application will be site response of a soft soil layer over stiffer rock. This scope keeps the project aligned with numerical differentiation, numerical solution of PDEs, stability, convergence, and practical civil engineering analysis.

The project will not attempt to build a full research-grade 2D seismic wave code. A full 2D elastic basin model requires additional complications: vector wavefields, material interpolation on staggered grids, free-surface treatment, absorbing boundaries, source radiation patterns, and more difficult validation. Those topics may be discussed as extensions, but they will not be required for the main result.

---

## Final selected project direction

The selected direction is:

> A 1D velocity-stress finite-difference simulation of vertically propagating shear waves through layered ground, with verification in a homogeneous medium and application to site amplification of a soft layer over rock.

This direction is the best choice because it satisfies four requirements at once:

1. It uses the numerical-method concepts taught in the course.
2. It is directly meaningful for earthquake and structural engineering.
3. It is simple enough to verify quantitatively.
4. It produces clear visualizations and engineering interpretation.

---

## Governing physical model

## 1D shear-wave site-response model

The model represents vertical propagation of horizontally polarized shear motion through a soil column. The independent variables are depth and time. The main unknowns are:

- horizontal particle velocity;
- shear stress.

The material properties are:

- density;
- shear modulus;
- shear-wave velocity.

The shear-wave velocity is derived from shear modulus and density. The model captures the mechanism most relevant to site amplification: waves slow down in soft soil, reflect at impedance contrasts, and can resonate within the soil layer.

## Why shear waves are the focus

Horizontal shear motion is especially relevant to earthquake engineering because it produces lateral demand on structures. A shear-wave site-response model is also easier to interpret than a full elastic P-SV simulation. It allows the project to focus on numerical implementation and engineering interpretation rather than becoming overloaded by advanced seismology details.

---

## Numerical method

## Explicit velocity-stress finite difference

The numerical scheme will update velocity and stress in alternating steps. Stress gradients update particle velocity, and velocity gradients update shear stress. This creates a direct numerical representation of the wave-propagation mechanism.

The implementation will use an explicit time-domain update. This means the next time level is computed directly from known values at the current or previous time levels. This is preferable for this project because it avoids matrix solvers and makes the stability condition transparent.

## Staggered grid

The model will use a staggered grid:

- velocity values are stored at one set of grid locations;
- shear stress values are stored halfway between velocity locations;
- material properties are assigned consistently with the locations where they are needed.

The staggered arrangement is useful because velocity and stress derivatives naturally require values between each other. It also gives a clean way to represent wave propagation in a velocity-stress formulation.

## Spatial accuracy

The project will implement a 2nd-order finite-difference scheme as the required method. This is the main solver.

A 4th-order spatial derivative may be implemented only after the 2nd-order solver is verified. If the 4th-order scheme is included, it will be presented as an extension and comparison, not as the foundation of the project.

This avoids overextension while still allowing a meaningful comparison of numerical dispersion if time permits.

---

## Computational domain

## Depth coordinate

The model domain is a vertical soil column extending from the ground surface to a depth below the soil-rock interface. The depth must be large enough to include:

- the soft soil layer;
- the stiffer rock below;
- a lower absorbing or damping region.

A practical domain depth is around 500 m to 1000 m, depending on the selected soil-layer thickness and input frequency range.

## Grid spacing

Grid spacing will be selected using the shortest wavelength that must be resolved. The shortest wavelength occurs in the slowest material at the highest frequency of interest.

The minimum rule for the required implementation will be:

- at least 10 grid points per shortest wavelength for the 2nd-order scheme;
- if the 4th-order extension is used, it may be tested at coarser resolution, but it must still be checked for numerical dispersion.

The grid-spacing selection will be documented clearly in the report.

## Time step

The time step will be selected from the Courant stability condition. The maximum wave speed in the model will control the time step, not the average wave speed.

The simulation will print or report the Courant number before running. If the Courant number is too high, the simulation setup is invalid and must be changed before results are used.

---

## Material models

## Model 1: homogeneous rock column

The homogeneous model is used for verification only. It has constant density, shear modulus, and shear-wave velocity throughout the column.

Purpose:

- verify that the wave travels at the correct theoretical velocity;
- check stability;
- check numerical error under grid refinement;
- provide a reference case for interpreting later results.

Expected behavior:

- one main shear-wave pulse travels through the medium;
- theoretical arrival time can be computed from distance divided by shear-wave velocity;
- numerical arrival time should approach the theoretical value as the grid is refined.

## Model 2: soft soil layer over rock

This is the main engineering model. A soft surface layer sits above a stiffer rock halfspace. The key parameters are:

- soil layer thickness;
- soil shear-wave velocity;
- soil density;
- rock shear-wave velocity;
- rock density.

Purpose:

- demonstrate reflection and transmission at an impedance contrast;
- compute surface amplification;
- compare amplification peaks to simple resonance theory.

Expected behavior:

- waves partially reflect at the soil-rock interface;
- surface motion in the layered model differs from the homogeneous rock reference;
- the amplification spectrum shows peaks near the layer resonance frequencies.

## Selected main parameter range

The main model will use a soft layer thickness that produces resonant frequencies inside the simulation bandwidth. A suitable target is a fundamental frequency between about 1 Hz and 5 Hz, because this range is practical for visualization and structural interpretation.

The soil layer should not be so deep that the fundamental frequency falls below the frequency resolution of the simulation. The simulation duration must be long enough to resolve the frequency range of interest.

---

## Boundary and input-motion treatment

## Ground surface

The top boundary will represent the ground surface. For the 1D shear-wave model, a stress-free surface condition will be used. This means shear stress at the surface is zero.

This boundary is physically important because the project output is surface motion.

## Bottom boundary

The lower boundary should reduce artificial reflections. A simple damping or absorbing layer near the base will be used. The report will state that this is an approximate treatment and that boundary reflections can affect late-time waveforms.

The bottom boundary should be placed far enough below the layer interface that the main response window is not dominated by artificial reflections.

## Input motion or source

The project will use a controlled pulse input rather than a complex earthquake record. A smooth pulse source is appropriate because it has known frequency content and is easier to verify numerically.

The source will be selected so that its frequency content overlaps the expected site-response frequencies. The source time history and its Fourier spectrum will be plotted.

---

## Engineering quantity of interest

## Surface response

The primary time-domain output is surface velocity or displacement. The project will record surface motion for both the homogeneous reference model and the layered soil model.

## Site amplification factor

The main frequency-domain output is the amplification factor:

> amplification at a frequency = amplitude spectrum of the layered-site surface response divided by amplitude spectrum of the reference-rock surface response.

The amplification factor shows which frequencies are amplified by the soft layer.

## Analytical resonance estimate

The report will compare numerical amplification peaks with the quarter-wavelength resonance estimate for a soft layer over stiff material. The analytical estimate will be used as a guide, not as an exact solution, because the numerical model includes finite duration, numerical damping, and approximate boundary treatment.

---

## Project stages

# Stage 1: Define scope and parameters

## Tasks

- Define the physical model: 1D vertical propagation of shear waves.
- Select baseline material properties for rock and soft soil.
- Select layer thickness so the fundamental resonance lies within the simulation bandwidth.
- Select grid spacing using the shortest wavelength criterion.
- Select time step using the Courant stability condition.
- Define the source time function and simulation duration.

## Outputs

- Parameter table.
- Short explanation of why the selected grid spacing and time step are valid.
- Initial diagram of the soil column and material layers.

## Visualization

- Layered velocity profile versus depth.
- Shear impedance profile versus depth.
- Source time function.
- Source frequency spectrum.

---

# Stage 2: Implement homogeneous-medium solver

## Tasks

- Initialize velocity and stress fields.
- Apply explicit velocity-stress updates on a staggered grid.
- Apply the selected source or input motion.
- Record motion at selected depths and at the surface.
- Run the homogeneous model first.

## Verification checks

- The wave should propagate without numerical blow-up.
- The maximum amplitude should remain finite.
- Arrival time at a receiver should match distance divided by shear-wave velocity.
- Reducing the time step or grid spacing should not qualitatively change the waveform.

## Outputs

- Homogeneous wave-propagation result.
- Arrival-time comparison table.
- Stability check table.

## Visualization

- Space-time plot of velocity showing the traveling pulse.
- Seismograms at multiple depths.
- Arrival-time marker plot comparing theoretical and numerical arrival time.

---

# Stage 3: Grid-refinement and numerical error study

## Tasks

- Repeat the homogeneous simulation with multiple grid spacings.
- Use the same physical model and source for all runs.
- Measure numerical arrival-time error at a fixed receiver location.
- Plot error versus grid spacing on a log-log plot.

## Interpretation

The error should decrease as grid spacing decreases. The exact observed slope may not perfectly match the theoretical order because arrival-time picking depends on threshold, source shape, and time-step resolution. The report should describe the observed trend honestly.

## Outputs

- Grid-refinement table.
- Error-versus-grid-spacing plot.
- Short discussion of stability and convergence.

## Visualization

- Log-log convergence plot.
- Overlay of receiver waveforms for different grid spacings.
- Optional plot of numerical arrival-time error versus grid spacing.

---

# Stage 4: Implement layered soil-over-rock model

## Tasks

- Create the soft-layer-over-rock material profile.
- Assign material properties to the staggered grid consistently.
- Run the simulation with the same source used in the reference model.
- Record surface motion and internal wave propagation.

## Expected physical behavior

- The incident wave interacts with the soil-rock interface.
- Part of the wave reflects and part transmits.
- Multiple reflections inside the soft layer can create resonance.
- Surface motion differs from the homogeneous rock reference.

## Outputs

- Layered-model surface response.
- Comparison between rock reference and layered site response.
- Explanation of reflection and resonance behavior.

## Visualization

- Velocity profile with soil-rock interface marked.
- Space-time velocity plot showing reflections between surface and interface.
- Surface seismogram for layered site.
- Overlay of reference and layered surface response.

---

# Stage 5: Site amplification analysis

## Tasks

- Compute the Fourier amplitude spectrum of the reference-rock surface response.
- Compute the Fourier amplitude spectrum of the layered-site surface response.
- Divide the layered-site spectrum by the reference spectrum to obtain amplification factor.
- Compute the simple resonance frequencies of the soft layer.
- Mark analytical resonance frequencies on the amplification plot.

## Expected results

The layered-site response should show frequency-dependent amplification. Peaks should occur near resonance frequencies, although exact agreement is not expected because the simulation uses finite duration, finite grid spacing, and approximate boundaries.

## Outputs

- Surface spectra for reference and layered models.
- Amplification factor plot.
- Table of analytical resonance frequencies.
- Engineering interpretation of which frequencies are amplified.

## Visualization

- Three-panel figure:
  - time-domain surface motion comparison;
  - Fourier amplitude spectra;
  - amplification factor with resonance markers.

---

# Stage 6: Optional 4th-order spatial derivative extension

## Purpose

The 4th-order extension will be included only if the main 2nd-order solver is already stable and verified.

## Tasks

- Replace the 2nd-order spatial derivative with a wider 4th-order stencil.
- Run the homogeneous verification case.
- Compare numerical error and waveform shape with the 2nd-order result.
- Compare computational cost qualitatively or through runtime measurements.

## Interpretation

The 4th-order scheme should reduce numerical dispersion for the same grid spacing. However, it uses a wider stencil and requires careful boundary handling. If results are unstable or difficult to interpret, the 4th-order scheme should be reported as a partially explored extension rather than a main conclusion.

## Visualization

- Waveform comparison between 2nd-order and 4th-order schemes.
- Optional dispersion plot showing numerical wave-speed error versus wavelength sampling.
- Optional runtime comparison table.

---

# Stage 7: Final report and interpretation

## Report structure

1. Abstract
2. Introduction and engineering motivation
3. Governing equation and numerical method
4. Finite-difference discretization
5. Stability and grid-selection criteria
6. Verification in homogeneous medium
7. Layered site-response model
8. Site amplification results
9. Discussion of limitations
10. Conclusion

## Main conclusions to support

The final report should support the following claims:

- The finite-difference method can convert a seismic wave PDE into a computable update rule.
- Stability and grid spacing are not optional details; they control whether the numerical result is meaningful.
- A soft soil layer can amplify certain frequency components of ground motion.
- The numerical amplification peaks can be interpreted using simple resonance theory.
- The implementation connects classroom numerical methods to earthquake-engineering site response.

---

## Visualization plan

## Required figures

1. Material profile figure
   - shear-wave velocity versus depth;
   - density or shear impedance versus depth;
   - soil-rock interface marked.

2. Source figure
   - source time function;
   - source amplitude spectrum;
   - dominant frequency marked.

3. Homogeneous verification figure
   - velocity space-time plot;
   - theoretical arrival-time line;
   - recorded receiver waveform.

4. Grid-refinement figure
   - log-log error plot;
   - optional waveform overlays for different grid spacings.

5. Layered model wave-propagation figure
   - space-time plot showing reflections;
   - surface and interface locations marked.

6. Surface motion comparison figure
   - reference-rock surface motion;
   - layered-site surface motion;
   - normalized overlay for shape comparison.

7. Site amplification figure
   - reference and layered spectra;
   - amplification factor;
   - analytical resonance markers.

## Optional animation

An animation may be created from the space-time wavefield by showing the velocity profile as it evolves with time. This is useful for presentation but should not replace quantitative figures.

Recommended animation content:

- vertical profile of particle velocity moving through the column;
- soil-rock interface marked;
- surface marked;
- time stamp shown on each frame;
- color or line amplitude scaled consistently.

---

## Validation and quality-control checklist

## Numerical checks

- Courant number is below the stability limit.
- Grid spacing resolves the shortest wavelength of interest.
- No NaN or infinite values occur during simulation.
- Homogeneous arrival time agrees with theoretical travel time within a reasonable tolerance.
- Grid refinement reduces arrival-time error.
- Boundary reflections are not dominating the main response window.

## Physical checks

- Wave speed is slower in the soft soil layer than in rock.
- Soil-rock impedance contrast produces reflected waves.
- Amplification peaks occur in a plausible frequency range.
- The source spectrum contains energy near the expected resonance frequencies.

## Reporting checks

- Every figure has a clear engineering purpose.
- Every numerical parameter is reported with units.
- Stability and grid-spacing choices are justified.
- Limitations are stated honestly.
- The report avoids claiming full 2D or full research-grade seismic modeling.

---

## Project limitations to state clearly

The final project should explicitly state these limitations:

- The model is 1D, so it does not represent basin-edge effects or 2D focusing.
- The input motion is idealized, not a full recorded earthquake ground motion.
- Material behavior is linear elastic.
- Damping and viscoelastic attenuation are not included in the main model.
- Boundary absorption is approximate.
- The analytical resonance estimate is a simplified guide, not an exact validation of every numerical detail.

These limitations do not weaken the project. They make the interpretation scientifically honest and help keep the project aligned with the course.

---

## Final deliverables

## Code deliverables

- Main simulation script or notebook.
- Parameter file or clearly defined parameter section.
- Functions for material model, source definition, finite-difference update, simulation loop, and plotting.
- Saved figures in a figures directory.

## Report deliverables

- Final technical report in journal-style format.
- Figures with captions.
- Tables for model parameters, grid-refinement results, and resonance frequencies.
- Short discussion connecting numerical method behavior to earthquake-engineering interpretation.

## Presentation deliverables

- One slide explaining the numerical method.
- One slide showing homogeneous verification.
- One slide showing the layered soil model.
- One slide showing site amplification results.
- One slide summarizing conclusions and limitations.

---

## Final project statement

This project implements a transparent finite-difference model for seismic shear-wave propagation in layered ground. It uses a controlled numerical setup to demonstrate discretization, stability, convergence, and wave amplification. The final result is an engineering interpretation of how local soil conditions can amplify ground motion at specific frequencies, supported by numerical simulation and simple analytical resonance theory.
