# Project File Structure and Stage-Wise Implementation

## Recommended structure

```text
cee490-computer-methods/
│
├── README.md
├── from_reference.md
├── working_plan.md
├── file_structure.md
├── requirements.txt
│
├── notebooks/
│   └── 01_run_project.ipynb
│
├── src/
│   └── site_response/
│       ├── __init__.py
│       ├── parameters.py
│       ├── materials.py
│       ├── source.py
│       ├── fd_solver_1d.py
│       ├── analysis.py
│       ├── visualization.py
│       └── utils.py
│
├── scripts/
│   ├── run_homogeneous.py
│   ├── run_layered_site.py
│   └── make_all_figures.py
│
├── figures/
│   ├── source_wavelet.png
│   ├── velocity_profile.png
│   ├── homogeneous_verification.png
│   ├── convergence_plot.png
│   ├── layered_wavefield.png
│   ├── site_amplification.png
│   ├── wave_animation.gif
│   └── animation_layered_reflection.gif
│
├── results/
│   ├── homogeneous_results.npz
│   ├── layered_results.npz
│   └── convergence_results.csv
│
├── report/
│   └── final_report.pdf
│
└── tests/
    └── test_basic_solver.py
```

---

## File roles

| File | Purpose |
|---|---|
| `README.md` | Short project overview, how to run, main outputs |
| `from_reference.md` | Useful ideas extracted from the reference text |
| `working_plan.md` | Technical implementation plan |
| `file_structure.md` | Project organization and stage-wise execution |
| `requirements.txt` | Python dependencies |
| `01_run_project.ipynb` | Clean demonstration notebook with explanation and final plots |
| `parameters.py` | Grid, time step, material constants, simulation settings |
| `materials.py` | Homogeneous and layered soil-rock models |
| `source.py` | Input pulse / Ricker wavelet and source spectrum |
| `fd_solver_1d.py` | Main velocity-stress finite-difference solver |
| `analysis.py` | Arrival-time check, convergence, spectra, amplification factor |
| `visualization.py` | All plots and animations |
| `utils.py` | Small helper functions |
| `scripts/` | Reproducible command-line runs |
| `figures/` | Saved plots and animation outputs |
| `results/` | Saved numerical arrays and tables |
| `tests/` | Simple sanity checks |

---

## Recommended execution style

Use `.py` files for the real implementation and one `.ipynb` notebook for presentation.

The notebook should not contain the main solver logic. It should import the functions from `src/site_response/`, run the main cases, and display the final results cleanly.

---

## Stage-wise implementation

## Stage 1: Project setup

Create the folder structure, install dependencies, and define all simulation parameters. Keep all constants in `parameters.py` so the model is easy to modify.

**Main files:**

- `requirements.txt`
- `parameters.py`
- `README.md`

**Outputs:**

- printed parameter summary;
- CFL check;
- grid-spacing check.

---

## Stage 2: Material models and source

Create the homogeneous reference model and the layered soil-over-rock model. Define the input pulse and plot its time history and frequency spectrum.

**Main files:**

- `materials.py`
- `source.py`
- `visualization.py`

**Outputs:**

- velocity profile plot;
- impedance profile plot;
- source wavelet plot;
- source spectrum plot.

---

## Stage 3: Homogeneous solver verification

Implement the 1D velocity-stress finite-difference solver and test it first in a homogeneous medium. Compare numerical arrival time with theoretical travel time.

**Main files:**

- `fd_solver_1d.py`
- `analysis.py`
- `run_homogeneous.py`

**Outputs:**

- homogeneous seismogram;
- arrival-time comparison table;
- homogeneous wavefield plot.

---

## Stage 4: Grid-refinement study

Run the homogeneous model with multiple grid spacings. Check whether the error decreases as the grid becomes finer.

**Main files:**

- `analysis.py`
- `run_homogeneous.py`

**Outputs:**

- convergence table;
- log-log error plot;
- waveform overlays for selected grid sizes.

---

## Stage 5: Layered site-response model

Run the soft-layer-over-rock model. Compare the surface response with the homogeneous rock reference.

**Main files:**

- `materials.py`
- `fd_solver_1d.py`
- `run_layered_site.py`

**Outputs:**

- layered wavefield plot;
- surface response comparison;
- reflection/reverberation visualization;
- `animation_layered_reflection.gif` showing wave reflection and reverberation at the soil-rock interface.

---

## Stage 6: Site amplification analysis

Compute the Fourier spectra of the reference and layered-site responses. Divide them to obtain the amplification factor. Mark analytical resonance frequencies.

**Main files:**

- `analysis.py`
- `visualization.py`
- `make_all_figures.py`

**Outputs:**

- spectra comparison plot;
- amplification factor plot;
- resonance frequency table.

---

## Stage 7: Animation and final presentation

Create animations to show wave propagation through the soil column. Use them for presentation, not as the main validation result.

**Main files:**

- `visualization.py`
- `01_run_project.ipynb`

**Outputs:**

- `wave_animation.gif`;
- `animation_layered_reflection.gif`;
- final notebook;
- final report figures.

---

## Final workflow

```bash
pip install -r requirements.txt
python scripts/run_homogeneous.py
python scripts/run_layered_site.py
python scripts/make_all_figures.py
```

Then open:

```text
notebooks/01_run_project.ipynb
```

Use the notebook to present the method, results, figures, animations, and conclusions.

---

## Final deliverables

- Clean source code in `src/site_response/`
- One readable notebook in `notebooks/`
- Saved figures and animations in `figures/`
- Saved numerical results in `results/`
- Final report in `report/`
- Short README explaining how to reproduce the project