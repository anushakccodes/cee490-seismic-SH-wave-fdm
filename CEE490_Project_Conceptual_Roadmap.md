# CEE 490 — Graduate Seismic Wave Project: Complete Conceptual Roadmap
## What You Are Doing, Why You Are Doing It, and How Everything Connects

**This document contains zero code.**
**It is a complete explanation of every concept, every decision, and every step.**
**Read this before touching any code.**

---

# PART 1 — THE BIG PICTURE

---

## What Is This Project, In Plain English?

Imagine an earthquake happens underground. The ground shakes. That shaking travels
outward in all directions — through rock, through soil, through sediment — until it
reaches the surface and shakes buildings.

Your project simulates exactly this process on a computer.

You create a virtual piece of the Earth — a 2D cross-section, 2 kilometres wide and
2 kilometres deep. You place a virtual earthquake source somewhere inside it. You
let the waves travel. You record what happens at different locations. You measure,
validate, and analyze the results.

This is called **seismic wave propagation simulation** and it is one of the most
important tools in earthquake engineering. It is used to:

- Predict how strongly different locations will shake during an earthquake
- Design buildings in earthquake-prone regions
- Understand why some cities suffer more damage than others
- Map the internal structure of the Earth

---

## Why Does This Matter for Civil Engineering?

The 1985 Mexico City earthquake is the most famous example of why this matters.

The earthquake epicentre was 400 kilometres away from Mexico City. At that distance,
most cities would feel only mild shaking. But Mexico City was built on an ancient
lake bed — a thick layer of very soft sediment sitting on hard rock.

The soft sediment acted like a bowl of jelly. When the seismic waves entered the
sediment, they got trapped and amplified. The ground shook 8 times more strongly
in the sediment than on nearby hard rock. Thousands of buildings collapsed.
Ten thousand people died.

Your site amplification analysis directly simulates this phenomenon. You will
compute exactly how much the soft sediment in your model amplifies ground shaking
compared to hard rock, and at what frequencies. This is not just an academic
exercise — it is the calculation that earthquake engineers use to write building codes.

---

## Why Is the Elastic Wave Equation Better Than the Acoustic One?

Most students in your class who are doing wave propagation will use the
**acoustic wave equation**. This is a simplification that only models pressure waves.
It treats the ground like a fluid — something that can be compressed but not sheared.

Real rock is an elastic solid. It can be both compressed AND sheared. This means
two fundamentally different types of waves exist:

**P waves (Primary waves / Compressional waves)**
- The ground moves back and forth in the same direction the wave travels
- Like a spring being compressed and stretched
- Travel fast: typically 3000–6000 m/s in rock
- Arrive first at a seismometer — hence "Primary"
- Less damaging to structures

**S waves (Secondary waves / Shear waves)**
- The ground moves side to side, perpendicular to the wave direction
- Like shaking a rope up and down
- Travel slower: typically 1700–3500 m/s in rock (always slower than P)
- Arrive second — hence "Secondary"
- MORE damaging to structures because they cause side-to-side shaking
- Cannot travel through fluids (this is how we know the Earth's outer core is liquid)

If you use the acoustic equation, you only simulate P waves. You miss S waves entirely.
S waves carry more of the destructive energy in earthquakes.

The **full elastic wave equation** gives you both P and S waves simultaneously,
plus a third type — Rayleigh surface waves — that emerges automatically from their
interaction at the free surface. This is the actual physics.

---

## What Makes Your Project Graduate Level?

There are five things that elevate this beyond a typical course project:

**1. Full elastic wave equation**
Five coupled equations instead of one. Two wave types instead of one.
This is the standard used in professional seismology software worldwide.

**2. Staggered grid (Virieux 1986)**
A special grid design where different physical quantities are stored at
slightly different spatial locations. This paper has over 2000 citations
and is the foundation of modern computational seismology.

**3. Mathematical proof of correctness**
The convergence analysis proves your code is correct through mathematics,
not just by looking at whether the pictures seem reasonable.

**4. Site amplification engineering application**
You answer a real civil engineering question: by how much does a
sedimentary basin amplify ground shaking, and does your simulation
agree with analytical theory?

**5. Dispersion analysis**
You quantify the error in your numerical wave speed analytically and
show why the 4th order scheme is more efficient than the 2nd order scheme.

---

# PART 2 — THE MATHEMATICS, EXPLAINED FROM SCRATCH

---

## What Is a Partial Differential Equation (PDE)?

An ordinary differential equation (ODE) has one independent variable — usually time.
A partial differential equation (PDE) has multiple independent variables — in your
case, time t, horizontal position x, and depth z.

The wave equation is a PDE because it describes how pressure or velocity changes
simultaneously in time AND in space:

```
How fast pressure changes in time  =  c²  ×  How curved pressure is in space
```

The "curvature in space" is what drives the wave forward. If pressure is higher on
your left than your right, it pushes you to the right. That push propagates outward
as a wave.

---

## The Full Elastic Wave Equation — What Each Part Means

The full elastic wave equation consists of five coupled equations.
Here is what each one means physically.

### The Two Equations of Motion

These come from Newton's second law: Force = Mass × Acceleration.

**Equation 1 — Horizontal motion:**
The rate of change of horizontal velocity at any point equals the sum of
forces pushing it horizontally, divided by the density.
The forces come from the normal stress gradient (pushing left-right) and
the shear stress gradient (pushing sideways-up or sideways-down).

**Equation 2 — Vertical motion:**
Exactly the same idea, but for vertical velocity.
Forces come from the shear stress gradient (horizontal) and normal stress
gradient (vertical).

### The Three Constitutive Relations (Hooke's Law)

These come from the elastic properties of the rock — how it responds to deformation.

**Equation 3 — Horizontal normal stress (txx):**
How fast the horizontal compression stress changes equals the elastic moduli
times the rate at which the material is being stretched horizontally and
compressed vertically. The term (lambda + 2*mu) is called the P-wave modulus
because it controls P-wave speed.

**Equation 4 — Vertical normal stress (tzz):**
Same idea for vertical compression. Symmetric with equation 3.

**Equation 5 — Shear stress (txz):**
How fast the shear stress changes equals the shear modulus (mu) times the
rate at which the material is being sheared. Shear stress is what drives S waves.

### The Material Parameters

**rho (density)** — mass per unit volume of the rock. Typical rock: 2200 kg/m³.
Higher density means heavier rock, waves travel slower.

**lambda (first Lamé parameter)** — relates to how a material resists compression
while being stretched sideways. Not directly measurable but computed from vp and vs.

**mu (shear modulus / second Lamé parameter)** — resistance of the material to
shearing. ZERO for fluids (fluids cannot be sheared). Positive for solids.
Directly controls S-wave speed. This is why S waves cannot travel through fluid.

**P-wave speed:** vp = sqrt((lambda + 2*mu) / rho)
**S-wave speed:** vs = sqrt(mu / rho)

Note: vs is ALWAYS less than vp. The ratio vp/vs = sqrt(3) ≈ 1.73
for a material with Poisson's ratio = 0.25 (typical rock).

---

## Why Five Equations Instead of One?

The acoustic equation has one unknown — pressure p.
The elastic equation has five unknowns — vx, vz, txx, tzz, txz.

This is because the elastic equation captures the full mechanical state of the
material at every point:
- Two velocity components describe HOW FAST the rock is moving (and in which direction)
- Three stress components describe HOW HARD the rock is being pushed and twisted

At every time step, you update all five of these simultaneously. They are COUPLED —
the velocity components drive the stress changes, and the stress gradients drive the
velocity changes. This back-and-forth coupling is what makes waves propagate.

---

## What Are Lamé Parameters and Where Do They Come From?

You measure vp and vs in the field using seismic surveys. From those measurements
you compute lambda and mu using:

```
mu     = rho × vs²
lambda = rho × vp² - 2 × mu
```

These are the intrinsic elastic constants of the rock. Every rock type has different
values. Soft sediment has small mu (easy to shear). Hard granite has large mu.
Fluids have mu = 0 (cannot be sheared at all).

Important: lambda must be positive for the material to be physically reasonable.
In very soft materials with unusual vp/vs ratios, you can get negative lambda from
the formula — you must clip these to zero or the simulation becomes unphysical.

---

## What Is the Staggered Grid and Why Is It Needed?

### The Problem With a Regular Grid

On a regular grid, all five quantities (vx, vz, txx, tzz, txz) would sit at the
same grid points. When you try to compute derivatives — for example, how txx
changes as you move in the x direction — you compute:

```
dtxx/dx ≈ [txx(i+1) - txx(i-1)] / (2*dx)
```

This derivative is evaluated AT point i, but the result is used to update vx
ALSO at point i. The derivative and the variable being updated are at the same
location. This creates numerical instabilities for the elastic equation.

### The Staggered Grid Solution

Virieux (1986) showed that the elastic wave equation is naturally stable and
accurate when different quantities are stored at slightly offset locations:

- vx lives BETWEEN normal stress points in the x direction
- vz lives BETWEEN normal stress points in the z direction
- txz lives BETWEEN all other quantities, shifted half a step in both directions

Now when you compute dtxx/dx to update vx, the txx values sit exactly half a
grid spacing on both sides of the vx point. The derivative is automatically a
perfect central difference. No instability. Maximum accuracy.

Think of it like a chessboard. Velocities sit on black squares.
Normal stresses sit on white squares. Shear stress sits between all of them.
Every derivative naturally spans from one colour to the other.

This is not a trick or a workaround — it is the mathematically natural way to
discretize the elastic wave equation. That is why it became the standard.

---

## The Finite Difference Method — Why We Need It

### Why We Cannot Solve the PDE Exactly

The elastic wave equation has known analytical solutions only for the simplest
possible cases — infinite homogeneous media, simple geometries, constant properties
everywhere. The moment you add:

- Multiple rock layers with different properties
- A sedimentary basin
- A fault zone
- Realistic topography

...the analytical solution no longer exists. Mathematicians have proven this.
You cannot write down a formula that gives the exact answer for realistic geology.

### What Finite Difference Does

Instead of solving the equation for every point in continuous space, you:

1. Place a grid of discrete points across your domain
2. Replace the continuous derivatives with approximations using nearby grid points
3. March forward in time step by step

The key approximation is replacing derivatives with differences:

**First derivative (how steeply something changes):**
Instead of the calculus limit, use the slope between two nearby grid points.
Central difference: slope ≈ (value at right neighbor - value at left neighbor) / (2 × spacing)

**Second derivative (how curved something is):**
Use three points: (right neighbor - 2×center + left neighbor) / spacing²

This approximation introduces a small error — the difference between the finite
difference formula and the exact derivative. This error shrinks as you use a
finer grid. The rate at which it shrinks defines the ORDER of the scheme.

---

## 2nd Order vs 4th Order — What the Order Means

### 2nd Order Spatial Scheme

Uses the three-point stencil for second derivatives:
(value at right - 2×center + value at left) / dx²

The error in this approximation is proportional to dx².
If you halve your grid spacing, the error quarters.

This scheme needs at least 10 grid points per wavelength to be accurate.
The wavelength of a 15 Hz wave in rock with vs = 1732 m/s is:
wavelength = vs / f = 1732/15 = 115 metres
So your grid spacing must be no more than 115/10 = 11.5 metres.

### 4th Order Spatial Scheme (What You Use)

Uses a five-point stencil with special coefficients:
(-1/12 × value at far right) + (4/3 × right) + (-5/2 × center) + (4/3 × left) + (-1/12 × far left)
divided by dx²

Where do these coefficients come from? Taylor series expansion to higher order.
By using more neighboring points, the error cancellation goes further.
The error is now proportional to dx⁴.

If you halve your grid spacing, the error reduces by a factor of 16.

This scheme only needs 5 grid points per wavelength. Using the same example:
maximum grid spacing = 115/5 = 23 metres.

So the 4th order scheme can use a grid 4× coarser than the 2nd order scheme
for the same accuracy. In 2D this means 16 times fewer grid points, which means
16 times less memory and much faster computation.

This is not just a minor improvement — it is what makes large-scale seismic
simulations computationally feasible.

---

## Leapfrog Time Integration — How Time Stepping Works

### What Time Stepping Means

You start at time zero with everything at rest. You advance in tiny steps of
dt = 0.0005 seconds. At each step, you use the current state of the wavefield
to compute the state at the next moment.

After 2000 steps you have simulated 1 second of wave propagation.

### The Leapfrog Scheme

The leapfrog scheme uses central difference in time. Instead of just looking
at the present to predict the future (Forward Euler), it uses both the present
AND the past:

```
future = 2 × present  -  past  +  (time step)² × (spatial driving force)
```

Where does this come from? Taylor series in time:

Expand pressure forward:  p(t+dt) = p(t) + dt×p'(t) + (dt²/2)×p''(t) + ...
Expand pressure backward: p(t-dt) = p(t) - dt×p'(t) + (dt²/2)×p''(t) + ...

Add them together. All odd terms cancel (dp/dt, d³p/dt³, etc.).
You are left with:
p(t+dt) + p(t-dt) = 2×p(t) + dt²×p''(t) + O(dt⁴)

Rearrange for p''(t) and substitute the wave equation: p''(t) = c²×Laplacian(p).
Solve for p(t+dt). That is the leapfrog formula.

### Why Leapfrog and Not Forward Euler, Backward Euler, or Crank-Nicolson?

**Forward Euler** evaluates the right hand side at the current time only.
Von Neumann stability analysis shows its amplification factor |G| is always
greater than 1 for the wave equation — meaning it always blows up regardless
of time step size. Unconditionally unstable. Cannot be used.

**Backward Euler** is unconditionally stable but requires solving a 40,000×40,000
sparse matrix system at every time step. It also artificially damps wave amplitudes
because |G| < 1 — the wave energy slowly disappears even in a perfectly elastic medium.

**Crank-Nicolson** has |G| = 1 exactly — energy conserved. But it is implicit,
requiring the same expensive matrix solve every time step.

**Leapfrog** also has |G| = 1 exactly when the CFL condition is satisfied. Energy
is perfectly conserved. AND it is explicit — no matrix solve needed. You compute
the future directly from the past and present.

The trade-off: leapfrog requires a time step small enough to satisfy the CFL
condition. But you would want a small time step anyway for temporal accuracy.
So this restriction costs you nothing in practice.

---

## The CFL Condition — What It Is and Why It Cannot Be Violated

CFL stands for Courant-Friedrichs-Lewy, the three mathematicians who derived it
in 1928.

### The Physical Intuition

The CFL condition says: your time step must be small enough that information
cannot travel more than one grid spacing per time step.

Think of it this way. A wave travels at speed vp. In one time step of dt seconds,
it moves vp×dt metres. If this distance is larger than your grid spacing dx, the
wave is trying to jump over grid points. The numerical scheme cannot track it.
The solution becomes unstable and grows without bound.

In 2D, the wave can travel diagonally — so the effective maximum distance in one
step is dx×sqrt(2). The condition becomes:

```
dt ≤ dx / (vp_max × sqrt(2))
```

CRITICAL: use vp_max — the MAXIMUM P-wave velocity anywhere in your model.
If your granite layer has vp = 5500 m/s, that is the speed that sets your
time step limit, even if most of your domain is slower rock at 3000 m/s.

### What Happens If You Violate It

The amplification factor |G| becomes greater than 1. Every time step multiplies
the wavefield by a number slightly larger than 1. After hundreds of steps, the
solution has grown exponentially. You see NaN or inf values in your arrays.
The simulation explodes.

There is no way to recover from a CFL violation except to reduce dt and restart.

### How to Check It

Before running any simulation, compute the Courant number:
```
r = vp_max × dt / dx
```
This must be less than 1/sqrt(2) ≈ 0.707.

Print this number at the start of every simulation. If it is above 0.707, stop.

---

## The Ricker Wavelet — What It Is and Why We Use It

### What a Source Is

In your simulation, the earthquake source is a point in the grid where you
inject energy at the start of the simulation. Every time step, you add a
small amount of energy to that point according to a time function called
the source time function.

### The Ricker Wavelet

The Ricker wavelet is the standard source time function in seismology. It is
the second derivative of a Gaussian bell curve. In time it looks like:
a bell shape that starts at zero, goes positive, peaks, comes back to zero,
goes slightly negative (the two sides of the second derivative), then
returns to zero permanently.

It is used because:
- It starts and ends at zero (no DC offset, no artificial step in the wavefield)
- It has a well-defined dominant frequency f0 (the peak of its spectrum)
- Its frequency content is bandlimited — it does not contain frequencies much
  higher than f0, which would require a finer grid to resolve

The dominant frequency f0 = 15 Hz in your project. This means the wavelength
of the slowest S wave (vs_min = 600 m/s) is:
wavelength = 600 / 15 = 40 metres

Your grid spacing of 5 metres gives 40/5 = 8 grid points per wavelength —
comfortably above the 5 minimum for the 4th order scheme.

### The Moment Tensor Source

A simple pressure source injects energy equally in all directions. But a real
earthquake fault does not radiate equally in all directions. The radiation
pattern depends on the fault geometry and slip direction — this is described
by the moment tensor.

The moment tensor has three components in 2D: Mxx, Mzz, Mxz.

- Mxx = Mzz = 1, Mxz = 0: explosive source, radiates equally everywhere
- Mxx = 1, Mzz = -1, Mxz = 0: vertical strike-slip fault, four-lobed pattern
- Mxx = 0, Mzz = 0, Mxz = 1: 45-degree thrust fault

You inject these by adding energy to the corresponding stress components (txx, tzz, txz)
instead of pressure. This is the physically correct way to model an earthquake source
in the elastic framework.

---

## Absorbing Boundaries — Why Waves Must Not Reflect From the Grid Edges

### The Problem

Your computational domain is finite — 2 km × 2 km. But in reality the Earth
extends for thousands of kilometres in every direction. Waves that reach your
grid edge should continue travelling into the distance and never return.

But on a computer, when a wave reaches the edge of the grid, it has nowhere to
go. Without special treatment, it reflects back into the domain — just like
sound echoes off a wall. These reflections contaminate your wavefield with
completely artificial waves that have no physical meaning.

### The Sponge Layer Solution

The simplest and most robust absorbing boundary is the sponge layer
(Cerjan et al., 1985). You designate the outermost 50 grid points on every
side as a sponge zone. Inside this zone, you multiply the wavefield by a
damping factor at every time step.

The damping factor starts at 1.0 at the inner edge of the sponge
(no damping — interior is unaffected) and decreases exponentially toward
0.0 at the outer edge (full absorption).

The wave enters the sponge, gets progressively attenuated, and by the time
it reaches the actual grid boundary it has essentially no amplitude left.
No reflection occurs because there is nothing left to reflect.

### How to Verify the Sponge Is Working

Plot the sponge mask itself. It should look like:
- Pure white (value = 1.0) in the interior
- Gradually darkening toward the edges
- Pure black (value ≈ 0.0) at the very corners

After running a simulation, look at the wavefield at late times. In a homogeneous
medium, there should be no waves reflecting back from the edges. If you see
circular waves bouncing back, your sponge is not thick enough — increase n_sponge
from 50 to 70 or increase the alpha parameter.

---

# PART 3 — THE VELOCITY MODELS, EXPLAINED

---

## Why the Velocity Model Is the Most Important Input

Everything in your simulation — how fast waves travel, where they reflect, where
they amplify, what the seismograms look like — depends on the velocity model.

The velocity model is a 2D array where every grid point stores the P-wave and
S-wave speeds at that location. Different geological materials have very different speeds:

| Material | Vp (m/s) | Vs (m/s) | Character |
|---|---|---|---|
| Air | 340 | 0 | Cannot support shear |
| Water | 1500 | 0 | Cannot support shear |
| Soft sediment | 1500–2000 | 400–800 | Very slow, amplifies waves |
| Weathered rock | 1800–2500 | 700–1200 | Near surface layer |
| Competent rock | 2500–4000 | 1400–2300 | Most of the crust |
| Hard granite | 5000–6000 | 2800–3500 | Very fast, reflects energy |
| Fault zone | 800–1500 | 300–600 | Damaged rock, very slow |

---

## Model 1 — Homogeneous Medium (For Validation Only)

Constant velocity everywhere. Every grid point has the same vp and vs.

This model exists purely for mathematical validation. In a homogeneous medium
the analytical solution is known: a wave from a point source at distance r
arrives at exactly time t = r/vp (for P waves) or t = r/vs (for S waves).

You compare your numerical arrival times to this exact formula. If they agree
within 3%, your code is correctly implemented.

Never present results from the homogeneous model as your main result. It is
a validation tool, not a physical model of anything real.

---

## Model 2 — Three-Layer Model (For Understanding Reflections)

Three horizontal layers with different velocities:
- Top 30% of domain: soft sediment (slow)
- Middle 40%: competent rock (medium)
- Bottom 30%: hard basement (fast)

This model demonstrates two fundamental wave phenomena:

**Reflection:** When a wave hits the interface between layers moving from slow
to fast, part of its energy reflects back upward. The reflected wave travels
back toward the surface. You will see this as a downward-travelling wave in
your snapshots that later appears as an upward-travelling wave.

**Transmission:** The other part of the energy continues downward into the
faster layer, but it travels faster and the wavefront bends (refracts).

**Mode conversion:** At each interface, some P-wave energy converts to S-wave
energy and vice versa. You will see extra wavefronts appearing at each layer
boundary that were not present in the homogeneous model.

---

## Model 3 — Realistic Geological Model (Your Main Result)

This model is designed to resemble a real geological scenario relevant to
earthquake engineering. It contains four features:

### Feature 1: Weathered Near-Surface Layer

The top 5% of the domain is a thin layer of very soft, weathered rock.
Vp = 1500 m/s, Vs = 600 m/s.

This represents the weathered, fractured zone that exists near the surface
in most geological settings. It is slower than the underlying rock because
the rock has been broken up by freeze-thaw cycles, chemical weathering, and
stress relief.

Effect on waves: slight amplification near the surface, important for
near-surface engineering applications.

### Feature 2: Sedimentary Basin

An elliptical zone of very soft sediment in the upper-middle part of the domain.
Vp = 1500 m/s, Vs = 500 m/s. The ellipse is about 800 m wide and 500 m deep.

This represents an alluvial valley or ancient lake bed — a bowl-shaped
depression that was filled with soft sediments over geological time. The
Po Plain in northern Italy, the Los Angeles basin in California, and the
Mexico City lake bed are all real examples of this type of structure.

Effect on waves: this is where your site amplification occurs. Waves entering
the basin slow down dramatically. Their wavelength shortens. They get trapped
inside the basin, bouncing between the basin walls and the hard rock below.
This trapping and focusing amplifies the wave energy at the surface.

### Feature 3: Hard Granite Basement

Below 60% of the domain depth, the velocity jumps sharply to vp = 5500 m/s,
vs = 3175 m/s. This represents the crystalline basement rock — old, hard,
unweathered granite or similar material.

Effect on waves: strong reflection at the top of the basement sends energy
back up into the softer layers above. The high velocity means waves travel
very fast in this layer and the wavefronts are widely spaced.

### Feature 4: Fault Zone

A narrow vertical channel of very slow, damaged rock running from near the
surface down to the top of the basement. Vp = 1000 m/s, Vs = 400 m/s.
Only about 40 metres wide (8 grid points).

This represents an active fault — a zone where the rock has been repeatedly
broken and ground up by earthquakes over millions of years. The damaged rock
is much weaker and slower than the surrounding intact rock.

Effect on waves: waves cannot travel efficiently through the fault zone.
They slow down dramatically, get partially reflected, and diffract around
the edges of the fault. This creates a characteristic shadow zone behind
the fault and diffracted waves that appear to come from the fault edges.

---

# PART 4 — THE ANALYSES, EXPLAINED

---

## Analysis 1 — Convergence Analysis (Proving Your Code Is Correct)

### The Core Question

Anyone can write code that produces output that looks reasonable. The question
is: is the output ACTUALLY correct? How do you prove it?

The answer in numerical methods is convergence analysis. It is the only
rigorous way to verify a numerical code.

### The Idea

Your 2nd order scheme has a spatial error proportional to dx².
Your 4th order scheme has a spatial error proportional to dx⁴.

This means:
- If you halve dx in the 2nd order scheme, the error should quarter
- If you halve dx in the 4th order scheme, the error should reduce by a factor of 16

You test this by running the simulation four times at different grid spacings:
20 m, 10 m, 5 m, and 2.5 m. Each time you measure the error in the P-wave
arrival time compared to the exact analytical solution (which exists for the
homogeneous model).

### The Log-Log Plot

You plot error on the vertical axis and grid spacing on the horizontal axis.
Both axes are logarithmic.

If your 2nd order scheme is correctly implemented, the points fall on a straight
line with slope = 2.
If your 4th order scheme is correctly implemented, the points fall on a straight
line with slope = 4.

This is a mathematical fingerprint of the order of accuracy. You cannot fake it.
If your slope is 1.3 or 3.1, something is wrong with your implementation.
If your slopes are 2.0 and 4.0, your code is proven correct.

### Why This Is So Powerful

A professor reading your report will see the convergence plot and immediately
know: this student did not just copy code from the internet. You cannot get the
right slope by accident. The slope is the mathematical proof that your discretization
is correct and that your implementation of it is correct.

This is called VERIFICATION — proving that you solved the equations correctly.
It is distinct from VALIDATION — proving that the equations you solved match reality.

---

## Analysis 2 — Validation: P and S Wave Arrival Times

### What You Are Doing

In the homogeneous model, P waves and S waves travel outward from the source
at constant known speeds. You place four receivers at known distances from the source.

For each receiver you know:
- P-wave should arrive at time = distance / vp (theoretical)
- S-wave should arrive at time = distance / vs (theoretical)

You read off the actual arrival times from your numerical seismograms and compare.

### What Good Agreement Looks Like

P and S arrival times within 3% of theoretical is excellent.
Within 1% is near-perfect for a course project.

If your P-wave arrivals are correct but S-wave arrivals are wrong, your shear
modulus mu is probably incorrect in your model.

If both are wrong by the same ratio, check your CFL condition — you may be
using the wrong value of vp_max.

### The Seismogram Itself

A seismogram from your elastic simulation should show:
- Quiet flat line until the P wave arrives
- A burst of oscillation (the P wave)
- Quiet period (the S-P time gap — longer for distant receivers)
- A larger burst of oscillation (the S wave, typically stronger)
- Gradual decay as surface waves arrive and the coda dies out

The time gap between P and S arrivals is one of the most important measurements
in seismology — it is how seismologists locate earthquakes. The farther the
station, the larger the S-P gap.

---

## Analysis 3 — Dispersion Analysis (Quantifying Numerical Error)

### What Is Grid Dispersion?

In the real physical world, all frequency components of a wave travel at the
same speed (in a non-dispersive medium). A sharp wavefront stays sharp as it
travels.

In a numerical simulation, different frequency components of the wave travel
at slightly different NUMERICAL speeds. This is called grid dispersion. It is
an artifact of the finite difference approximation, not real physics.

The effect is that a sharp wavefront gradually smears out as it travels. High
frequency components of the wave travel slightly slower numerically than low
frequency components. The wavefront develops a tail that was not there physically.

### How You Quantify It

The ratio of numerical wave speed to true wave speed can be computed analytically:

c_numerical / c_true is a function of the sampling ratio dx/wavelength.

When dx/wavelength is small (fine grid), this ratio is close to 1. When dx/wavelength
is large (coarse grid), this ratio deviates significantly from 1.

You plot this ratio for both 2nd order and 4th order schemes. The 4th order scheme
stays close to 1 for much larger values of dx/wavelength. This is the quantitative
explanation for why 4th order needs only 5 grid points per wavelength while 2nd
order needs 10.

### What Your Plot Shows

The x-axis is dx/wavelength (sampling ratio). At 0, the grid is infinitely fine —
no dispersion. At 0.5, you are at the Nyquist limit — aliasing destroys the wave.

The 2nd order curve drops away from 1.0 at around dx/wavelength = 0.1 (10 points/wavelength).
The 4th order curve stays close to 1.0 until about dx/wavelength = 0.2 (5 points/wavelength).

At the vertical line marking dx/wavelength = 0.1, the 4th order scheme has essentially
no dispersion error while the 2nd order scheme already has noticeable error.
This is the quantitative reason for the 4th order scheme's superiority.

---

## Analysis 4 — Site Amplification (The Engineering Application)

### The Physical Question

Two buildings are built 300 metres apart. One sits on hard rock. One sits on soft
sediment above hard rock (a basin). An earthquake happens. Which building shakes more?
By how much? At what frequencies?

This is the site amplification problem. It is the central question in
earthquake site response analysis, which is required by building codes in
seismically active regions.

### The Method: Spectral Ratio Technique

You place two receivers at the same distance from the source but at different
geological locations:
- One receiver on hard rock (reference site)
- One receiver inside the soft sediment basin (study site)

You record the vertical velocity seismogram at both locations.
You take the Fourier transform of each seismogram.
You divide: amplification(f) = |spectrum at basin| / |spectrum at rock|

This ratio tells you: at frequency f, waves at the basin site have this many
times more amplitude than at the rock site. An amplification factor of 5 at
frequency 3 Hz means 3 Hz waves are 5 times stronger in the basin than on rock.

### The Analytical Check: 1D Resonance Theory

A simple 1D model of a soft soil layer over hard rock has known analytical
resonant frequencies:

```
f_n = (2n-1) × vs_basin / (4 × H_basin)
```

where n = 1, 2, 3... (mode number), vs_basin is the shear wave velocity
in the soft layer, and H_basin is the depth of the soft layer.

These are the frequencies at which standing S waves can fit inside the basin
as quarter-wavelength resonances. At these frequencies, waves get trapped and
amplified dramatically.

For your basin (vs = 500 m/s, H = 200 m):
- Fundamental frequency: f₁ = 1 × 500 / (4 × 200) = 0.625 Hz
- First overtone: f₂ = 3 × 500 / (4 × 200) = 1.875 Hz
- Second overtone: f₃ = 5 × 500 / (4 × 200) = 3.125 Hz

Your numerical amplification spectrum should show peaks near these frequencies.
If it does, your 2D simulation agrees with 1D analytical theory — that is your
validation of the engineering result.

If the peaks are slightly shifted from the analytical values, that is
physically reasonable. The 1D formula assumes a flat-bottomed, infinite basin.
Your basin is elliptical with curved boundaries — the 2D geometry shifts the
resonances slightly. Discuss this discrepancy in your report — it shows you
understand the limitations of the analytical model.

### Why This Matters for Engineering

Building codes specify that structures must withstand a certain level of ground
shaking at their location. The shaking level depends strongly on the soil type.

In the United States, the ASCE 7 building code defines site classes (A through F)
based on the average shear wave velocity of the top 30 metres of soil. Sites
on soft soil (Class E or F) must be designed for much stronger shaking than
sites on hard rock (Class A or B).

Your simulation directly computes the amplification factor that justifies these
site classifications. You are computing the engineering basis for building code requirements.

---

## Analysis 5 — Performance Analysis (Cost vs Accuracy)

### The Question

If you have a fixed computational budget — say 10 minutes of simulation time —
which scheme gives you the most accurate result?

The naive answer is: 4th order is more accurate per grid point, so use it.
But it also uses more grid points per stencil, so each step costs more.

The correct answer requires comparing accuracy per unit cost.

### The Key Insight

To achieve 1% error in wave speed:
- 2nd order scheme needs dx = 2.5 m → 800×800 = 640,000 grid points
- 4th order scheme needs dx = 5.0 m → 400×400 = 160,000 grid points

The 4th order scheme uses 4 times fewer grid points for the same accuracy.
Even though each grid point update is slightly more expensive (9-point stencil
vs 5-point stencil), the overall simulation is faster and uses less memory.

You measure this directly: run both schemes to the same accuracy target and
record the wall clock time and memory usage. Print the ratio.

This transforms a theoretical observation ("4th order converges faster") into
a practical engineering result ("4th order is X times faster for the same accuracy").

---

# PART 5 — THE FIVE-DAY PLAN, CONCEPTUALLY

---

## Day 1 — Set Up the Foundation

**What you are doing:**
Building the scaffolding. You are not simulating anything yet.
You are defining the physical and numerical parameters, building the grid,
creating the velocity models, and verifying that everything is set up correctly.

**The key decisions you make:**
- Grid spacing dx = 5 m (chosen to give 8 points per wavelength at the
  minimum S-wave speed of 600 m/s with f0 = 15 Hz)
- Time step dt = 0.0005 s (chosen to satisfy CFL with vp_max = 5500 m/s)
- Domain 2000 × 2000 m (large enough that waves don't reach boundaries during simulation)
- Simulation time 1.0 s (long enough to see P and S waves separate clearly)

**What you verify:**
- CFL condition: courant number = vp_max × dt / dx must be below 0.707
- Spatial sampling: vs_min / (f0 × dx) must be above 5
- All three velocity models build without errors
- Lame parameters are all positive (no unphysical values)
- Grid arrays initialize correctly: five arrays of shape (400, 400)

**What you produce:**
Plots of all three velocity models — two panels each showing P-wave and S-wave
velocity. The realistic geology model should clearly show the basin (blue/cyan
for low velocity), the competent rock background (green/yellow), the hard
basement (orange/red), and the narrow fault zone (very dark).

---

## Day 2 — Build the Physical Machinery

**What you are doing:**
Writing the mathematical engines. The FD stencils, the time integration schemes,
the source, and the absorbing boundaries.

**The staggered grid FD update:**
This is the heart of the project. You implement two versions:
- 2nd order: uses neighbors at ±1 grid spacing, error O(dx²)
- 4th order: uses neighbors at ±1 and ±2 grid spacings, error O(dx⁴)

The update proceeds in two sub-steps every time step:
Sub-step A: update velocities from stress gradients (Newton's law)
Sub-step B: update stresses from velocity gradients (Hooke's law)

**The source:**
The Ricker wavelet is computed as a mathematical function of time at each step.
Its value is injected into the stress components at the source grid point using
the moment tensor. For an explosive source you add equal amounts to txx and tzz.

**The sponge:**
You build a mask array once before the loop. Every grid point in the sponge
zone has a value between 0 and 1. After each time step update, you multiply
all five field arrays by this mask. Values outside the sponge zone are 1.0
so they are unaffected. Values inside decay the wavefield to zero.

**What you verify:**
- Run one time step manually and check no NaN or inf values appear
- Plot the Ricker wavelet: should look like a symmetric pulse peaking at f0
- Plot the Fourier transform of the Ricker wavelet: should peak at f0 = 15 Hz
- Plot the sponge mask: white interior, dark edges

---

## Day 3 — Run the Simulations

**What you are doing:**
Running the actual wave propagation simulations and looking at the results
for the first time. This is the day you get to see your waves.

**The time loop:**
The simulation is a loop over 2000 time steps. At each step:
1. Apply the 4th order staggered FD update (compute new vx, vz, txx, tzz, txz)
2. Inject the source (add Ricker wavelet amplitude to stress at source point)
3. Apply the sponge (multiply all five fields by the mask)
4. Record seismograms (save the velocity values at each receiver location)
5. Save snapshots every 100 steps for visualization

**What you should see:**
In the homogeneous model at early times (t = 100 ms):
- One circular ring expanding outward from the source
- At this point P and S waves are close together and hard to distinguish

At moderate times (t = 300 ms):
- Two distinct circular rings
- The outer ring is the P wave (faster)
- The inner ring is the S wave (slower)
- The ratio of their radii should be vp/vs = sqrt(3) ≈ 1.73

In the realistic geology model:
- The circular wavefronts distort when they enter regions of different velocity
- Waves slow down in the soft basin (wavefronts bunch up)
- Waves speed up in the hard basement (wavefronts spread out)
- Additional wavefronts appear at boundaries (reflections and mode conversions)
- A shadow zone develops behind the fault zone

**What you verify:**
- No NaN or inf values appear in any time step
- Two distinct wavefronts are visible in the homogeneous model
- The ratio of wavefront radii equals vp/vs
- Seismograms show two arrivals (P then S)

---

## Day 4 — Prove Everything Is Correct

**What you are doing:**
Converting visual observations into rigorous quantitative measurements.
This day is what separates engineering science from casual simulation.

**Convergence analysis:**
Run the homogeneous simulation four times at different grid spacings.
Measure arrival time error each time. Plot on log-log axes. Measure slopes.
2nd order should give slope ≈ 2.0. 4th order should give slope ≈ 4.0.
If your slopes are right, your code is mathematically proven correct.

**P and S wave validation:**
For each receiver in the homogeneous model, compare numerical arrival times
to the exact formula (distance/speed). Build a table with receiver distance,
analytical arrival time, numerical arrival time, and percentage error.
Target: all errors below 3%.

**Dispersion analysis:**
Compute the numerical phase velocity ratio analytically for both schemes.
Plot as a function of sampling ratio. Show that 4th order maintains accuracy
to finer wavelengths (coarser grids) than 2nd order.

**Site amplification:**
Run the realistic model. Compare seismograms at basin and rock receivers.
Compute the spectral ratio. Plot amplification factor vs frequency.
Mark the analytical resonant frequencies. Check for agreement.

**Performance analysis:**
Run both schemes to achieve the same 1% accuracy. Record time and memory.
Compute ratios. Show 4th order is more efficient overall.

---

## Day 5 — Write the Report and Finalize

**What you are doing:**
Transforming all your results into a coherent scientific document structured
as a journal paper.

**The report is not a lab report.**
It is not a narrative of what you did day by day. It is a scientific
communication structured around a question, methods, results, and interpretation.

**Abstract:**
150-200 words. State the problem, methods, key quantitative findings
(your convergence slopes, your peak amplification factor, your agreement
with analytical resonance theory), and main conclusion.
A good abstract can be read in isolation and give the reader the full picture.

**Introduction:**
Why does seismic wave simulation matter? What are P and S waves and why
does the elastic equation capture more physics than acoustic? What is the
Virieux (1986) staggered grid and why is it the standard? What does your
project contribute that goes beyond a standard implementation?

**Mathematical Formulation:**
Write out all five coupled elastic PDEs. Define every symbol. Derive the P
and S wave speeds from the Lamé parameters. Explain the staggered grid layout
with a diagram. Write out the 2nd and 4th order FD stencils with their
exact coefficients. Derive the CFL condition from Von Neumann stability analysis.
Show the leapfrog update formula and where it comes from (Taylor series derivation).

**Implementation:**
Table of all grid parameters. Description of each velocity model and the
geological features it represents. Source formulation (Ricker wavelet + moment tensor).
Absorbing boundary method.

**Validation:**
The convergence plot with your measured slopes. The arrival time table with
percentage errors. This section proves your code is correct.

**Results:**
Your wavefield figures. P/S wave separation in the homogeneous model.
Complex wavefield evolution in the realistic model. Seismograms. Animation stills.

**Application: Site Amplification:**
The three-panel figure (seismograms, spectra, amplification factor).
The table of analytical resonant frequencies. Discussion of agreement
or discrepancy between numerical and analytical results.

**Discussion:**
Compare the three methods quantitatively. Explain why 4th order is more
efficient (dispersion analysis provides the quantitative basis). Discuss
what the site amplification results mean for earthquake engineering.
Acknowledge limitations: your model is 2D (real Earth is 3D), your source
is simplified (real faults have complex geometry), your boundaries are
approximate (sponge is not perfect).

**Conclusion:**
One paragraph. What was done, what was found, what it means.

**References:**
At minimum: Moczo et al. (2004), Virieux (1986), Levander (1988),
Cerjan et al. (1985), Aki and Richards (2002).

---

# PART 6 — HOW TO EXPLAIN YOUR PROJECT TO YOUR PROFESSOR

---

## If Asked: "What Did You Implement?"

"I implemented the full 2D elastic wave equation using the velocity-stress
formulation on a staggered grid following the Virieux (1986) scheme. This
gives me both P waves and S waves simultaneously, with realistic mode conversion
at geological interfaces. I implemented two spatial schemes — 2nd order and 4th order
finite differences — and compared them. Time integration used the leapfrog
(central difference) scheme throughout."

---

## If Asked: "How Do You Know Your Code Is Correct?"

"I performed a convergence analysis. I ran the simulation at four grid spacings
on a homogeneous medium where the exact analytical solution is known. For the
2nd order scheme, the error in P-wave arrival time decreased as dx² — I measured
a slope of 2.0 on the log-log plot. For the 4th order scheme, the error decreased
as dx⁴ — I measured a slope of 4.0. These slopes match the theoretical orders
of accuracy, which mathematically verifies the implementation."

---

## If Asked: "Why Is This More Powerful Than Acoustic?"

"The acoustic equation treats the ground as a fluid — it models only compressional
P waves. Real rock is an elastic solid that supports both compressional and shear
deformation. The full elastic equation captures both P waves and S waves simultaneously.
S waves are slower, arrive second, and are typically more damaging to structures
because they cause transverse ground motion. They also cannot propagate through
fluids, which is how we know the Earth's outer core is liquid. Mode conversion
at interfaces — where P waves partially become S waves and vice versa — is a
physical phenomenon that only appears in the elastic formulation. Additionally,
Rayleigh surface waves emerge naturally from the simulation without any additional
programming."

---

## If Asked: "What Is the Engineering Relevance?"

"I performed a site amplification analysis on a model containing a sedimentary
basin. Using the spectral ratio method — dividing the Fourier spectrum of motion
at the basin site by the spectrum at a reference rock site — I computed the
frequency-dependent amplification factor. The peak amplification occurred near
the fundamental resonant frequency of the basin, which I computed analytically
as vs/(4H). The agreement between my numerical simulation and the 1D analytical
theory validates the engineering result. This calculation is directly analogous
to the site response analysis required by modern building codes for structures
in seismically active regions."

---

## If Asked: "Why Is 4th Order Better Than 2nd Order?"

"The 4th order scheme is not just more accurate — it is more EFFICIENT. For a
target accuracy of 1% error in wave speed, the 2nd order scheme requires
10 grid points per wavelength while the 4th order scheme requires only 5. In 2D,
this translates to 4 times fewer grid points overall — 4 times less memory and
substantially faster computation. My dispersion analysis quantified this: I
computed the numerical phase velocity ratio c_numerical/c_true as a function of
the spatial sampling ratio dx/wavelength for both schemes. The 4th order scheme
maintains this ratio within 1% of 1.0 for dx/wavelength up to 0.2, while the
2nd order scheme requires dx/wavelength below 0.1 for the same accuracy. My
performance benchmarks confirmed that the 4th order scheme was faster overall
despite the more expensive stencil."

---

# PART 7 — COMMON MISTAKES AND HOW TO AVOID THEM

---

## Mistake 1: Using Average Velocity for the CFL Condition

Always use the MAXIMUM P-wave velocity anywhere in your model for the CFL check.
If your granite basement has vp = 5500 m/s, that sets your time step limit.
Using the background velocity of 3000 m/s will underestimate the required
restriction and your simulation will likely blow up in the high-velocity region.

## Mistake 2: Forgetting to Use vs_min for the Spatial Sampling Criterion

The spatial sampling criterion applies to the SLOWEST wave, which is the S wave
in the softest material. If your fault zone has vs = 400 m/s, that determines
how fine your grid must be, not the average S-wave velocity.

## Mistake 3: Negative Lamé Parameters

When you compute lambda = rho×vp² - 2×mu, you can get negative values if the
vp/vs ratio in some zone is below sqrt(2) ≈ 1.414. This is unphysical (violates
thermodynamic stability). Always clip lambda to a minimum of zero.

## Mistake 4: Source Outside the Sponge Interior

Place your source well inside the interior, away from the sponge zone. If your
source is inside the sponge, the source energy gets immediately attenuated by
the sponge mask and your wavefield will have very low amplitude.

## Mistake 5: Comparing Schemes on Different Models

When comparing 2nd order vs 4th order, always use exactly the same model,
same grid, same source, same receivers. The only thing that should differ
is the stencil. Otherwise you cannot isolate the effect of the scheme.

## Mistake 6: Convergence Analysis on Heterogeneous Models

The convergence analysis must be done on the homogeneous model. Heterogeneous
models have no exact analytical solution to compare against, so you cannot
compute a meaningful error. The homogeneous model is the only case where you
know the exact answer.

## Mistake 7: Only Showing vz or Only Showing vx

For the elastic simulation, always show both velocity components side by side.
P waves appear predominantly in vz (vertical motion) for a vertical source.
S waves appear predominantly in vx (horizontal motion). Showing both panels
simultaneously demonstrates that you have actually simulated both wave types.

## Mistake 8: Discussing the Amplification Plot Without the Analytical Lines

The amplification factor plot is only scientifically meaningful if you overlay
the analytical resonant frequencies. Without those lines, it is just a wiggly
curve with no reference. With them, it is a validation of your simulation
against theory.

---

# SUMMARY: WHAT YOUR PROJECT PROVES

When your professor reads your project, they will see a student who:

1. **Understands the physics** — knows why elastic is more physical than acoustic,
   can explain P waves, S waves, mode conversion, and site amplification

2. **Understands the mathematics** — derived the leapfrog update from Taylor series,
   knows where the 4th order coefficients come from, can explain the CFL condition
   from Von Neumann stability analysis

3. **Knows how to verify numerical code** — convergence analysis with measured slopes
   matching theoretical orders of accuracy is the gold standard of verification

4. **Connects simulation to engineering practice** — site amplification is a real
   calculation used in earthquake engineering and building codes

5. **Can communicate science** — journal paper format with Abstract, clear figures,
   quantitative tables, and a discussion that goes beyond describing results
   to interpreting their meaning

No other student in the class will have all five of these. That is why this project
will be outstanding.
