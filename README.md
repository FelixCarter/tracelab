# TraceLAB

**A MATLAB toolbox for interindividual synchrony analysis of facial expression and head movement data acquired via Trace**

[![DOI](https://img.shields.io/badge/DOI-10.3390%2Fe28050503-blue)](https://doi.org/10.3390/e28050503)

---

## Overview

TraceLAB is a MATLAB toolbox for preprocessing and analyzing facial landmark data acquired through [Trace](https://pavlovia.org), a research media player implemented in PsychoPy's online platform. The toolbox provides researchers with tools to quantify interindividual synchrony in facial expressions and head movements as a measure of shared information processing.

TraceLAB implements:

- **Preprocessing**: Format, align, filter, and clean raw Trace data
- **Head Movement Synchrony**: Surrogate Synchrony (SUSY) for single-channel time series
- **Expression Synchrony**: Correlated Component Analysis (CorrCA) for multivariate facial landmark data
- **Visualization**: Publication-ready plotting functions

---

## Installation

1. **Download** the TraceLAB source code from this repository (clone or download as ZIP)
2. **Add to MATLAB path**:
   ```matlab
   addpath(genpath('/path/to/tracelab'));
   ```
3. **(Optional)** To permanently add TraceLAB to your path:
   ```matlab
   savepath
   ```

No additional toolboxes are required. TraceLAB was developed in MATLAB R2024a.

---

## Quick Start

### 1. Preprocess your Trace data

```matlab
cfg.dir = '/path/to/your/trace/csv/files/';
data = t_preproc(cfg);
```

### 2. Choose your analysis

**Head movement synchrony (SUSY):**
```matlab
cfg.datatype = 'HeadMovement';
cfg.segment = 3;
[results, extra] = t_susy(data, cfg);
summary = t_susy_summarise(results);
t_plot_susy(summary);
```

**Facial expression synchrony (CorrCA):**
```matlab
results = t_corrca(data);
t_multiplot_corrca(results);
```

**Head movement magnitude:**
```matlab
summary = t_movement_summary(data);
t_plot_movement(summary);
```

---

## Documentation

Full documentation is available in the [GitHub Wiki](../../wiki). Each function also contains detailed comments in the source code describing input parameters.

For a high-level overview of the main functions, see the [paper](https://doi.org/10.3390/e28050503).

---

## Functions

| Function | Description |
|----------|-------------|
| `t_preproc` | Preprocess raw Trace data: format, align, filter, and clean |
| `t_susy` | Surrogate Synchrony analysis for head movement or emotion scores |
| `t_susy_summarise` | Aggregate SUSY results across lags, segments, dyads, or participants |
| `t_corrca` | Correlated Component Analysis for facial expression synchrony |
| `t_movement_summary` | Summarize head movement data at various aggregation levels |
| `t_plot_movement` | Visualize head movement timecourses or bar plots |
| `t_plot_susy` | Visualize SUSY results with surrogate baseline |
| `t_multiplot_corrca` | Combined CorrCA visualization (topoplots + boxplot + timecourse) |
| `t_singleplot_corrca` | Individual CorrCA component visualization |

---

## Citing TraceLAB

If you use TraceLAB in your research, please cite:

> Carter, F., Richardson, M., Stanton Fraser, D., & Gilchrist, I.D. (2026). TraceLAB: A MATLAB Toolbox for Interindividual Synchrony Analysis of Facial Expression and Head Movement Data Acquired via Trace. *Entropy*, 28(5), 503. [https://doi.org/10.3390/e28050503](https://doi.org/10.3390/e28050503)

### Trace (Data Acquisition Tool)

> Levordashka, A., Richardson, M., Hirst, R. J., Gilchrist, I. D., & Fraser, D. S. (2025). Trace: A research media player measuring real-time audience engagement. *Behavior Research Methods*, 57(1), 44. [https://doi.org/10.3758/s13428-024-02522-0](https://doi.org/10.3758/s13428-024-02522-0)

---

## License

This project is open source and available under the [Creative Commons Attribution (CC BY) license](https://creativecommons.org/licenses/by/4.0/).

---

## Data Availability

Example data used in the paper is available at: [https://doi.org/10.17605/OSF.IO/E8ZJ2](https://doi.org/10.17605/OSF.IO/E8ZJ2)

---

## Contact

For questions or issues, please open a [GitHub Issue](../../issues) or contact the corresponding author: f.carter@bristol.ac.uk
