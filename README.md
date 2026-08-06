# DuneBio PoC: Hybrid ODE + PINN + DNA Tumor Growth Inhibition Model

## Overview

This repository demonstrates DuneBio's **toolbox approach** for quantitative systems pharmacology (QSP):

> "We don't use ODE, PINN, or DNA separately. We combine them hierarchically to reduce prediction error and capture patient-specific responses."

**Problem**: EGFR TKI (Osimertinib) response variation in NSCLC patients with different EGFR mutations (Exon 19 del vs L858R).

**Goal**: Reconstruct and predict tumor size changes using:
- ODE (mechanistic backbone)
- PINN (residual correction)
- DNA mutation data (patient-specific parameter tuning)

## Project Structure

```
dunebio_poc/
├── data/
│   ├── tumor_growth_data.csv          # Public clinical data (time, tumor size, dose)
│   └── egfr_mutation_profile.csv      # Synthetic/public DNA expression data (patient mutation status)
├── models/
│   ├── 01_ode_tgi.py                  # TGI ODE (tumor growth - drug kill equation)
│   ├── 02_pinn_corrector.py           # PINN corrects ODE residuals (physics loss + data loss)
│   └── 03_dna_modulator.py            # Maps DNA mutation info to ODE parameters (drug sensitivity IC50)
├── notebooks/
│   └── 01_hybrid_simulation.ipynb     # ODE-only vs ODE+PINN vs ODE+PINN+DNA comparison
├── outputs/
│   └── figures/
│       ├── hybrid_result.png          # ★ Final output: Overlaid prediction curves from 3 methods
│       └── residual_plot.png          # Visualization of ODE residuals corrected by PINN
└── README.md                          # Step-by-step explanation of how the 3 components are combined
```

## Workflow

| Phase | Activity | Files |
| :--- | :--- | :--- |
| 1 | Load public data & fit ODE to population average | `01_ode_tgi.py` |
| 2 | Compute ODE residuals & train PINN to correct them | `02_pinn_corrector.py` |
| 3 | Split patients by mutation status & derive IC50 modulation factors | `03_dna_modulator.py` |
| 4 | Compare 3 methods (ODE-only, ODE+PINN, Hybrid) & visualize | `01_hybrid_simulation.ipynb` |

## How the Hybrid Works

1. **ODE** provides the mechanistic backbone (Simeoni tumor growth model)
2. **PINN** learns and corrects residuals where ODE fails to capture complex dynamics
3. **DNA** modulates drug sensitivity (IC50) based on mutation status:
   - Exon 19 del → IC50 × 0.7 (increased sensitivity)
   - L858R → IC50 × 1.0 (baseline)

**Final prediction** = `ODE + PINN_correction + DNA_modulation`

## Requirements

```bash
pip install -r requirements.txt
```

<!-- TODO: Check versions before running -->
### Key Dependencies
- Python 3.8+
- TensorFlow 2.x (GPU recommended)
- DeepXDE
- NumPy, SciPy, Matplotlib

### Quick Start

```bash
# Clone the repository
git clone https://github.com/SunKimBusan/dunebio-poc.git
cd dunebio-poc

# Install dependencies
pip install -r requirements.txt

# TODO: Create scripts/run_poc.py
# Run the full pipeline
python scripts/run_poc.py

# Or explore interactively
jupyter notebook notebooks/01_hybrid_simulation.ipynb
```

### Expected Outputs

- `outputs/figures/hybrid_result.png`: Overlay of 3 prediction methods vs actual data
- `outputs/figures/residual_plot.png`: Visualization of PINN-corrected residuals
- Console output: RMSE comparison table

<!-- TODO: check -->
## Results

| Method | RMSE | Improvement |
| :--- | :--- | :--- |
| ODE-only | TBD | Baseline |
| ODE + PINN | TBD | TBD % |
| ODE + PINN + DNA (Hybrid) | TBD | TBD % |

*Expected: Hybrid shows lowest RMSE, demonstrating the value of combining all three components.*

## License

MIT

<!-- TODO: change mail account -->
## Contact

DuneBio Inc. - [your-email@example.com](mailto:your-email@example.com)
