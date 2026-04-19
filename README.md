# Dimuon Analysis Framework

This framework is designed for analyzing low-mass dimuon spectra, focusing on the yield calculation of $\omega$ and $\phi$ mesons. It provides a structured pipeline from data preparation to background estimation, signal extraction, and efficiency correction.

## Directory Structure

- `src/`: Contains C++ source files (ROOT) for each step of the analysis.
- `include/`: Shared header files defining the analysis interface.
- `macros/`: Compiled binaries for different analysis configurations.

## Analysis Workflow

The analysis follows a sequential numbered pipeline:

1. **Step 00 (Main Macros)**: Entry points for various analysis types (Data, MC, and specific cut studies).
    - `00_dimuon_macro_forEM.cxx`: Standard EM analysis including correction.
    - `00_dimuon_macro_forEM_mc.cxx`: Monte Carlo analysis for efficiency/acceptance.
2. **Step 01 (Histogramming)**: 
    - `01_Making1DmassfromHnSparse.cxx`: Extracts 1D/2D distributions from `THnSparse`.
3. **Step 02 (Preprocessing & Background)**:
    - `02_LikeSignMethod.cxx`: Combinatorial background estimation.
    - `02_CalAcceptance_Efficiency.cxx`: Calculates Acceptance, Efficiency, and Acceptance $\times$ Efficiency (1D/2D).
4. **Step 03 (Signal Fitting)**:
    - `03_PeakFit_CrystalBall_pol4.cxx`: Crystal Ball with 4th-order polynomial background.
5. **Step 04 (Yield Calculation)**:
    - `04_YieldCalcuration_*.cxx`: Calculates raw yields.
6. **Step 05 (Correction & Cross Section)**:
    - `05_correction.cxx`: Applies Acceptance $\times$ Efficiency and Luminosity corrections to calculate the production cross section.

## Build Instructions

To compile the analysis framework, ensure that ROOT is installed.

### 1. Standard Analysis (Full Chain)
```bash
g++ src/00_dimuon_macro_forEM.cxx \
    src/01_Making1DmassfromHnSparse.cxx \
    src/02_LikeSignMethod.cxx \
    src/03_PeakFit_CrystalBall_pol4.cxx \
    src/04_YieldCalcuration_CrystalBall_pol4.cxx \
    src/05_correction.cxx \
    -Iinclude `root-config --cflags --libs` -o macros/run_analysis
```

### 2. MC Efficiency Analysis
Includes Acceptance, Efficiency, and Acc $\times$ Eff calculations.
```bash
g++ src/00_dimuon_macro_forEM_mc.cxx \
    src/01_Making1DmassfromHnSparse.cxx \
    src/02_Check_mc_mass.cxx \
    src/02_CalAcceptance_Efficiency.cxx \
    -Iinclude `root-config --cflags --libs` -o macros/run_analysis_mc
```

### 3. Production Cross Section (Independent Compilation)
If you need to compile only the production cross section part:
```bash
g++ src/00_dimuon_macro_forEM_production.cxx \
    src/01_Making1DmassfromHnSparse.cxx \
    src/02_LikeSignMethod.cxx \
    src/03_PeakFit_CrystalBall_pol4.cxx \
    src/04_YieldCalcuration_CrystalBall_pol4.cxx \
    src/05_correction.cxx \
    -Iinclude `root-config --cflags --libs` -o macros/run_production
```

## How to Run

1. **Calculate Acceptance & Efficiency**:
   ```bash
   ./macros/run_analysis_mc
   ```
   This generates `AcceptanceEfficiency.root` containing `Omega/AccEff_Omega_Pt`, etc.

2. **Calculate Production Cross Section**:
   Ensure `AcceptanceEfficiency.root` is available (as `acceptance_weight.root` or updated in the code).
   ```bash
   ./macros/run_analysis
   ```

## Features
- **Acc x Eff Correction**: Separate calculation of Acceptance (Acc/Gen) and Efficiency (Reco/Acc).
- **Production Cross Section**: Automated calculation using corrected yields and luminosity.
- **2D Maps**: Support for 2D ($p_T, y$) acceptance and efficiency maps.

## Dependencies
- ROOT (CERN) 6.xx
- C++11 or higher
