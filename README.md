# CEE 490 Computer Methods: Finite-Difference Site Response

This repository implements a staged 1D finite-difference model for vertically propagating SH waves in layered ground. The current implementation covers Stages 1, 2, and 3 from `docs/file_structure.md`:

1. project setup and parameter checks;
2. material models, source wavelet, and visualization style;
3. homogeneous velocity-stress solver verification using analytical travel time.

## Repository layout

```text
cee490-computer-methods/
├── README.md
├── requirements.txt
├── docs/
├── notebooks/
│   └── 01_run_project.ipynb
├── src/
│   └── site_response/
├── scripts/
├── figures/
├── results/
└── tests/
```

## Run the notebook

```bash
pip install -r requirements.txt
jupyter notebook notebooks/01_run_project.ipynb
```

The notebook is the main presentation file. It imports the implementation from `src/site_response/`, explains the governing formulas in Markdown cells, and displays the Stage 1--3 outputs.

## Run the command-line verification

```bash
python scripts/run_homogeneous.py
```

This saves numerical arrays in `results/` and figures in `figures/`.

## Implemented outputs

- parameter summary;
- CFL and grid-spacing checks;
- shear-wave velocity profile;
- impedance profile;
- Ricker source time history;
- source frequency spectrum;
- homogeneous receiver seismogram;
- arrival-time comparison table;
- homogeneous depth-time wavefield image.

## Core equations

The solver uses the 1D SH-wave velocity-stress system:

```text
rho(z) dv/dt = d tau/dz
 d tau/dt   = mu(z) dv/dz
```

The homogeneous verification compares numerical arrival time against:

```text
t_theory = travel_distance / Vs
```

The stability check uses:

```text
C = Vs_max * dt / dz
```
