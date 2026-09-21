# Optical Dimers Simulation Codebase
# Optical Dimers

This repository contains the local Python implementation for simulating and
analysing two coupled optical emitters (an optical dimer). The main workflows
are provided as Jupyter notebooks in `local_version/`.

## Project layout

```text
local_version/
   optical_dimers.ipynb       Main simulation and analysis notebook
   entanglement.ipynb         Two-photon and entanglement analysis
   debugger.py                Script for interactive/debugging workflows
   sys_params.py              Constants and system configuration
   system_ops.py              System operators and vector utilities
   br_tensor.py               Bloch-Redfield tensor construction
   rates.py                   Spectral and transition-rate calculations
   spectral_dens.py           Spectral-density functions
   analysis.py                Core, emission, and entanglement analyses
   files/dipole_positions/    Dipole-position input data
```

## Requirements

The code uses Python and scientific-computing packages including:

- NumPy
- SciPy
- pandas
- Matplotlib
- Plotly
- Jupyter/IPython

Install the required packages in the Python environment selected by VS Code or
Jupyter. The repository does not currently include a pinned requirements file.

## Running the notebooks

The modules in `local_version` use imports such as `from sys_params import ...`.
For this reason, launch Jupyter with `local_version` as the working directory:

```bash
cd local_version
jupyter notebook
```

Then open either notebook:

- `optical_dimers.ipynb` for the main optical-dimer calculations
- `entanglement.ipynb` for two-photon density matrices and entanglement analysis

The notebooks contain the parameter choices, simulation steps, plotting, and
analysis specific to the calculation being run.

## Analysis classes

`analysis.py` provides three layers of analysis:

- `CoreAnalysis`: builds the Liouvillian propagator, initializes the density
   matrix, and prepares state projections.
- `OpticalDimerAnalysis`: calculates directional intensity and delayed second-
   order correlations.
- `EntanglementAnalysis`: calculates unintegrated two-photon density matrices.
- `QuantumMetrologyAnalysis`: Under development...

The reusable modules can also be imported from scripts launched inside
`local_version`:

```python
from analysis import CoreAnalysis, OpticalDimerAnalysis, EntanglementAnalysis
from br_tensor import BlochRedfieldCalculator
from rates import CombinedRates
from sys_params import ConfigParameters, SystemParameters
```

## Input data

Dipole-position configurations are stored in
`local_version/files/dipole_positions/`. The available CSV files describe the
supported relative dipole orientations, including horizontal, parallel,
opposite, orthogonal, and 45-degree configurations.
