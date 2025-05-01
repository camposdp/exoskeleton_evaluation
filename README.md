# EMG Data Processing and Ergonomic Risk Analysis

This repository contains MATLAB scripts developed for processing surface electromyography (sEMG) data to assess muscle activation and ergonomic risk across different sessions of a physical intervention study. The scripts were developed by **Dr. Daniel Campos** and **João Pedro Moreto Lourenção**.

If you use or adapt any part of this code for your own work, **please cite the authors** appropriately.

---

## 🚀 Workflow Overview

**To fully reproduce the analysis, follow these steps in order:**

1. **`extract_control.m`**  
   ➔ Processes baseline (control) EMG data collected **without the exoskeleton**.  
   ➔ **Purpose:** Builds the reference dataset to compare with the intervention.

2. **`extract_features_days.m`**  
   ➔ Processes all intervention recordings to extract dynamic and static RMS, MCV, and percentiles (5%, 25%, 50%, 75%, 95%).  
   ➔ **Purpose:** Prepares individual day data for later statistical and ergonomic analyses.

3. **`statistics.m`**  
   ➔ Compiles and visualizes boxplots and weekly summaries (mean ± SEM), performing **statistical tests (Mann–Whitney)** between control and intervention weeks.  
   ➔ **Purpose:** Explores trends and identifies significant differences over time.

4. **`ergotarget.m`**  
   ➔ Generates the **Ergo-Target** plots, comparing **best/worst days** to control data within a **risk model** (based on logistic fit and reference curves).  
   ➔ **Purpose:** Visualizes ergonomic risk zones (high, medium, low) and positions your data within these risk zones.

5. **`export_tables.m`**  
   ➔ Compiles and **exports results tables** (weekly means, ANOVA tests, percentiles) to `.csv` and Excel files for reporting.  
   ➔ **Purpose:** Final step to export structured tables for publication or further analysis.

---

## 🔍 Script Details and Theory

### 1️⃣ `extract_control.m`

- **What it does:**  
  Processes **three baseline days** of data to build control metrics:
  - Dynamic EMG (rms over pulses)
  - Static EMG (isometric)
  - Percentiles of normalized EMG (5–95%)

- **Theoretical basis:**  
  - Bandpass filtering: 20–350 Hz  
  - Notch filtering: 59–61 Hz (powerline noise)  
  - Dynamic RMS is extracted using a **double-threshold method** combined with morphological dilation to identify activation pulses.  
  - Signals are **normalized by MCV (Maximum Voluntary Contraction)** to express EMG as %MCV.

- **Key parameters:**  
  - `fs = 1000` Hz (sampling rate)  
  - Butterworth filters (`butter(2, ...)` for bandpass & notch)  
  - `minDurationMs = 500` ms (minimum valid pulse length)

---

### 2️⃣ `extract_features_days.m`

- **What it does:**  
  Loops over **30 intervention sessions**, extracting:
  - MCV (shoulder and biceps)
  - Static RMS  
  - Dynamic RMS per pulse  
  - Percentiles (5%, 25%, 50%, 75%, 95%) of normalized EMG.

- **Theoretical basis:**  
  - Same filtering and normalization as control data.  
  - The dynamic EMG is segmented into activation pulses using a **double-threshold strategy** (same as in `extract_control.m`).  
  - RMS per pulse is limited to the **first 10 pulses** for consistency.

- **Key parameters:**  
  - Input file convention:  
    - `mcv-ombro-XX.txt`, `mcv-biceps-XX.txt`  
    - `isometrica-XX.txt`  
    - `dinamica-XX.txt`  
  - Filter settings identical to control script.

---

### 3️⃣ `statistics.m`

- **What it does:**  
  - Compiles all intervention data + control for:  
    - **Boxplots** (per day) of dynamic RMS  
    - **Weekly means ± SEM** plots  
    - **Statistical tests (Mann–Whitney U test)** comparing each week vs. control.  
  - Also computes **mean & SEM** for dynamic and static EMG signals.

- **Theoretical basis:**  
  - Statistical comparison of dynamic muscle activation between intervention and control.  
  - Weekly averages to smooth variability and observe trends.  
  - Mann–Whitney U is used as a **non-parametric** test to compare distributions.

- **Key outputs:**  
  - Figures:  
    - Boxplots  
    - Evolution of weekly RMS (%MCV)  
  - Table: `controle_vs_semanas_completo.xlsx` with full statistical results.

---

### 4️⃣ `ergotarget.m`

- **What it does:**  
  - Visualizes **dynamic risk maps** using **logistic curves** as reference thresholds:  
    - High Risk (Red zone)  
    - Medium Risk (Yellow zone)  
    - Low Risk (Green zone)

  - Plots the **percentiles (5–95%)** of EMG activation for:  
    - Control data  
    - Worst session (highest mean RMS)  
    - Best session (lowest mean RMS)

- **Theoretical basis:**  
  - Inspired by ergonomic standards that **quantify exposure risk** using %MCV and probability thresholds.  
  - Uses **logistic fitting** to generate smooth risk boundaries based on quantiles.  
  - Visual inspection of performance relative to risk curves.

- **Key parameters:**  
  - `quantiles = [5, 25, 50, 75, 95]`  
  - Logistic curve fitting to reference points (e.g., x1/y_ref and x2/y_ref).

---

### 5️⃣ `export_tables.m`

- **What it does:**  
  - Aggregates everything into exportable **CSV and Excel sheets**, including:  
    - Weekly mean RMS + SEM  
    - Control means + % difference  
    - ANOVA p-values + significant pairs  
    - Weekly percentiles (5–95%)

- **Why it's important:**  
  - Facilitates reporting to supervisors or inclusion in papers.  
  - Includes **ANOVA analysis** to complement the non-parametric tests.

- **Key outputs:**  
  - `Resultados_RMS_Semanais.csv`  
  - `Resultados_ANOVA.csv`  
  - `Percentis_Semanais_Ombro.csv`  
  - `Percentis_Semanais_Biceps.csv`  
  - Combined Excel: `Relatorio_EMG.xlsx`

---

## ⚙️ Requirements

- **MATLAB** (tested on R2022b and later)
- Signal Processing Toolbox
- The EMG data files must follow the naming convention as shown in the scripts (`dinamica-XX.txt`, etc.).

---

## 📄 Citation

If you use or adapt this codebase, **please cite:**

