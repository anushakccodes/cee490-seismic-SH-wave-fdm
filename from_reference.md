# Reference-Based Technical Notes for the Project

## Purpose

This document extracts the finite-difference ideas that are directly useful for the project and translates them into concrete project decisions. Page references use the printed page numbers in the source text.

The selected project direction is a finite-difference simulation of seismic wave propagation and site response in layered ground. The emphasis is not on building a large research-grade seismic code. The emphasis is on using classroom numerical-method concepts in a physically meaningful civil engineering application: discretization, stability, convergence, boundary treatment, source definition, and engineering interpretation of numerical results.

---

## 1. Why finite differences are appropriate here

### Source reference

- Introduction, pp. 1-2
- The Principle of the Finite-Difference Method, p. 3

### Useful details

The reference explains that realistic Earth and near-surface ground models include heterogeneity, material discontinuities, and complex interfaces. Analytical solutions are generally unavailable for such models, so the differential problem must be transformed into an algebraic problem that can be solved numerically.

The finite-difference method is presented as a grid-point method: the computational domain is covered by a space-time grid, field variables are stored at grid points, and derivatives are approximated by finite-difference formulas. The method is especially relevant to seismic wave propagation because it can represent heterogeneous media, local material changes, and wave motion in the time domain.

### Use in this project

The project will use finite differences because the governing equation is a wave-propagation PDE and because the engineering quantity of interest, site amplification, depends on material layering and wave reflections. A finite-difference implementation is also directly aligned with the course goals: converting differential equations into computable algebraic updates, checking stability, and interpreting numerical error.

---

## 2. The project should not be treated as a black-box simulation

### Source reference

- Preface, p. vii
- Introduction, pp. 1-2

### Useful details

The reference warns that finite-difference codes can produce outputs that look plausible but are not necessarily accurate or physically meaningful. It emphasizes that users must understand the assumptions, stability restrictions, and limitations of the method.

### Use in this project

The final project will explicitly include verification checks rather than only wavefield plots. The report should show that the method is stable, that the numerical wave speed is reasonable, and that the computed site amplification is interpreted through known theory. This avoids a purely visual demonstration and makes the project defensible as a numerical-methods project.

---

## 3. Start from a controlled 1D wave-propagation problem

### Source reference

- Preface, p. vii
- Contents: 1D Elastic Problem, pp. 13-50
- Numerical Examples: Unbounded Homogeneous Elastic Medium, p. 113; Single Layer over Halfspace, p. 121

### Useful details

The source text is intentionally focused on the 1D elastic problem as a training foundation. It states that understanding the 1D problem is the beginning, but that it is useful because it helps develop understanding for higher-dimensional wave propagation.

The numerical examples include an unbounded homogeneous elastic medium and a single layer over a halfspace. These two examples match the most useful validation and engineering cases for the project:

- a homogeneous medium for numerical verification;
- a soft layer over stiffer material for site-response behavior.

### Use in this project

The project will use a 1D vertically propagating shear-wave model as the quantitative core. This is the best scope because it is rigorous, implementable within the course context, and directly connected to earthquake engineering site response. A full 2D elastic basin model will not be the main deliverable because it adds substantial boundary, source, material-interpolation, and validation complexity.

---

## 4. Grid definition and discretization choices

### Source reference

- Grid, pp. 3-4
- The FD Approximations to Derivatives, pp. 4-7

### Useful details

The reference defines a computational grid in space and time. The spatial grid spacing controls how finely the physical domain is sampled, and the time step controls how the solution advances in time. It distinguishes forward, backward, and central differences, and it shows how Taylor expansion is used to derive finite-difference formulas and truncation errors.

The reference also presents the idea of higher-order approximations. A 4th-order approximation uses a wider stencil than a 2nd-order approximation and reduces truncation error more rapidly as grid spacing decreases.

### Use in this project

The project will implement a uniform vertical grid for the 1D model. The main method will be an explicit velocity-stress update. A 2nd-order scheme will be implemented first because it is easier to verify and explain. A 4th-order spatial derivative may be added after the 2nd-order implementation is stable, but the report will not depend on the 4th-order scheme unless it is fully verified.

The report will include a concise derivation or explanation of how central finite differences approximate derivatives and how smaller grid spacing reduces truncation error.

---

## 5. Consistency, stability, and convergence are required project checks

### Source reference

- Properties of the FD Equation, pp. 7-9
- Von Neumann's Analysis of Stability, p. 22

### Useful details

The reference identifies consistency, stability, and convergence as central properties of a finite-difference model. Consistency means the finite-difference equation approaches the original PDE as time step and grid spacing go to zero. Stability means bounded physical solutions remain bounded numerically. Convergence means the finite-difference solution approaches the exact PDE solution as the discretization is refined.

The reference also notes the practical connection between consistency, stability, and convergence through the Lax equivalence idea: for a well-posed linear problem, a consistent and stable finite-difference scheme is convergent.

### Use in this project

The project will include the following checks:

1. A stability check using a Courant number before each simulation.
2. A homogeneous-medium arrival-time check against the theoretical travel time.
3. A grid-refinement study showing that the numerical error decreases as grid spacing decreases.
4. A clear explanation that convergence behavior is evidence of correct numerical implementation, not proof of real-site accuracy.

---

## 6. Use an explicit time-domain scheme

### Source reference

- Explicit and Implicit FD Schemes, p. 9

### Useful details

The reference distinguishes explicit and implicit schemes. In explicit schemes, the field value at a new time level is calculated from known values at previous time levels and material parameters. In implicit schemes, field values at a new time level are solved simultaneously through a matrix system. The reference notes that explicit schemes are computationally simpler and widely used in earthquake ground-motion modeling.

### Use in this project

The project will use an explicit time-domain finite-difference scheme. This is the best choice for the course because it makes the update rule transparent, avoids large matrix solvers, and makes stability restrictions easy to demonstrate. An implicit method will be discussed only as context, not implemented.

---

## 7. Use a velocity-stress elastic formulation

### Source reference

- Coordinate System and Basic Quantities, p. 13
- Equation of Motion and Hooke's Law, p. 13
- Velocity-stress FD Schemes, p. 14

### Useful details

The reference introduces the elastic problem through equation of motion and Hooke's law. In velocity-stress form, motion is described by particle velocity and stress. This formulation is especially useful for wave propagation because stress gradients drive particle velocity, and velocity gradients update stress.

For a 1D shear-wave model, the practical governing variables are particle velocity and shear stress. The shear-wave speed is controlled by shear modulus and density.

### Use in this project

The project will use a 1D velocity-stress formulation for vertically propagating shear waves:

- particle velocity represents horizontal ground motion;
- shear stress represents internal resistance to shear deformation;
- shear modulus and density define the material properties;
- shear-wave speed is derived from material properties.

This is physically meaningful for earthquake engineering because horizontal shear motion is a key contributor to structural demand.

---

## 8. Material parameters must be placed carefully on the grid

### Source reference

- Material Grid Parameters, p. 29
- Staggered-grid FD Schemes - A Summary, p. 30
- Contact of Two Media - A Material Discontinuity, pp. 33-40

### Useful details

The reference treats material parameters as grid quantities and emphasizes that heterogeneous media and material discontinuities require care. At layer interfaces, properties such as density and elastic modulus change abruptly, and the numerical scheme must handle those changes without creating artificial behavior.

The staggered-grid approach places different field variables at offset locations. This is useful for velocity-stress formulations because derivatives naturally connect stress and velocity positions.

### Use in this project

The project will use a staggered-grid layout in the 1D model:

- velocity will be stored at one set of grid locations;
- shear stress will be stored at staggered locations between velocity points;
- material properties will be assigned consistently with those field locations.

For the layered site-response model, the interface between the soil layer and rock halfspace will be placed on the grid deliberately. The report will state how material properties are assigned at or near the interface. If averaging is used, it will be described clearly.

---

## 9. Layered media are better suited than a complex 2D geology for the main deliverable

### Source reference

- Homogeneous and Heterogeneous FD Schemes, pp. 9-10
- Contact of Two Media - A Material Discontinuity, pp. 33-40
- Numerical Example: Single Layer over Halfspace, p. 121

### Useful details

The reference distinguishes between simple media and media with material discontinuities. It also shows that even a single layer over a halfspace is enough to produce important wave phenomena: impedance contrast, reflection, transmission, and resonance.

### Use in this project

The main engineering model will be a soft soil layer over a stiffer rock halfspace. This model is simple enough to validate and explain, but rich enough to produce site amplification. It is also more appropriate for a numerical-methods course than a complex 2D basin with many features that would be difficult to validate.

The model will be interpreted as a simplified site-response problem, not as a full geological reconstruction.

---

## 10. Boundary conditions are a central numerical issue

### Source reference

- Free Surface, p. 43
- Wave Excitation, p. 46
- Boundaries of the Grid, p. 47

### Useful details

The reference separates physical boundaries, source excitation, and artificial grid boundaries. A seismic simulation requires a treatment of the free surface, a source or input motion, and a strategy for preventing artificial reflections from finite computational boundaries.

### Use in this project

The 1D model will use:

- a stress-free top boundary to represent the ground surface;
- a source or input motion applied near the base of the model;
- an absorbing or damping treatment near the bottom to reduce artificial reflections.

The report will include a short explanation that finite domains create artificial boundaries and that boundary treatment affects the reliability of late-time waveforms.

---

## 11. Wave excitation must be controlled and documented

### Source reference

- Wave Excitation, p. 46
- Program SOURTF, p. 109

### Useful details

The reference includes wave excitation as a separate part of the numerical model. The source time function controls the frequency content of the simulation. Since finite-difference accuracy depends on the shortest wavelength represented, the source frequency must be consistent with the grid spacing and material wave speeds.

### Use in this project

The project will use a smooth pulse source, such as a Ricker wavelet, because it has controlled frequency content and avoids unrealistic sudden jumps. The dominant frequency will be selected to excite the expected site-response frequency range. The source time function and its Fourier spectrum will be plotted and discussed.

---

## 12. Site amplification should be the main civil engineering application

### Source reference

- Numerical Example: Single Layer over Halfspace, p. 121
- Contact of Two Media - A Material Discontinuity, pp. 33-40

### Useful details

A soft layer over stiffer rock naturally creates impedance contrast and wave trapping. This is the core mechanism behind frequency-dependent site amplification. The simple layered model allows comparison against the quarter-wavelength resonance estimate.

### Use in this project

The project will compute surface motion for:

1. a homogeneous rock reference model;
2. a soft layer over rock model.

The amplification factor will be computed as a spectral ratio between the layered-site surface response and the reference-rock surface response. Analytical resonance frequencies will be marked on the amplification plot.

---

## 13. Realistic attenuation is valuable but should not be implemented in the main project

### Source reference

- Incorporation of the Realistic Attenuation, pp. 51-74

### Useful details

The reference devotes a substantial section to viscoelastic attenuation. This confirms that attenuation is important for realistic seismic-wave modeling, but it also shows that attenuation requires additional constitutive modeling and extra state variables.

### Use in this project

Viscoelastic attenuation will not be implemented in the main project. It will be identified as a limitation and possible future extension. This keeps the project focused on the finite-difference method, stability, convergence, wave propagation, and site amplification.

---

## 14. Numerical examples should guide the validation strategy

### Source reference

- Unbounded Homogeneous Elastic Medium, p. 113
- Single Layer over Halfspace, p. 121
- Unbounded Homogeneous Viscoelastic Medium, p. 129

### Useful details

The reference examples suggest a natural progression:

1. first test a homogeneous elastic medium;
2. then test a layer-over-halfspace model;
3. treat viscoelasticity as a more advanced extension.

### Use in this project

The project will follow the first two examples conceptually:

- Homogeneous elastic medium: verify arrival time and convergence.
- Layer over halfspace: compute and interpret site amplification.

The viscoelastic example will be acknowledged but excluded from implementation.

---

## Final project decisions derived from the reference

| Topic | Selected decision | Reason |
|---|---|---|
| Main model | 1D velocity-stress shear-wave propagation | Strong connection to source text, course topics, and site response |
| Numerical method | Explicit finite difference | Transparent, implementable, and directly tied to stability analysis |
| Grid | Uniform staggered grid | Natural for velocity-stress wave propagation |
| Verification | Homogeneous arrival time and grid refinement | Directly tied to consistency, stability, and convergence |
| Engineering application | Soft layer over rock site amplification | Practical and analytically interpretable |
| Boundary treatment | Stress-free surface plus absorbing base treatment | Physically meaningful and numerically necessary |
| Main visualization | Source spectrum, velocity model, wavefield evolution, seismograms, amplification spectrum | Shows both numerical method and engineering interpretation |
| Excluded from main scope | Full 2D elastic basin, moment tensor source, viscoelastic attenuation, topography | Too much implementation risk for the course project |

---

## How the reference will appear in the report

The report should cite the source text in the following places:

- finite-difference method overview: pp. 1-4;
- derivative approximations and truncation error: pp. 4-7;
- consistency, stability, and convergence: pp. 7-9;
- explicit schemes: p. 9;
- velocity-stress elastic formulation: pp. 13-14;
- stability analysis: p. 22;
- material parameters and staggered grids: pp. 29-30;
- material discontinuity and layer interface treatment: pp. 33-40;
- free surface, wave excitation, and grid boundaries: pp. 43-47;
- single layer over halfspace example: p. 121.
