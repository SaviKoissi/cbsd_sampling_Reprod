# CBSD Sampling in Nigeria

This work represents an effort to reproduce the study by Ferris et al. 2024 *(https://doi.org/10.1371/journal.pone.0304656)* CBSD surveillance study. While the article is publicly available, no fully reproducible code or computational pipeline was provided to replicate the analyses directly.

As a result, we developed an independent implementation based solely on the methodological descriptions reported in the paper. This pipeline aims to mirror the original modeling and simulation framework as closely as possible given the available information.

We note, however, that certain assumptions and implementation details were not explicitly specified in the original publication. Consequently, some components of this reproduction may reflect necessary interpretations rather than exact replication. Users of this code and associated results should therefore exercise caution and consider this work as an approximate reconstruction rather than a definitive reproduction of the original study.

## 🧠 High-level architecture

You’re looking at a two-scale stochastic simulation framework:

Landscape-scale model (grid-based, Nigeria-wide)
Within-field model (plant-level dynamics)
Surveillance strategy layer
Evaluation + replication

A good R project would look like this:

```bash
/project
 ├── data/                # REAL DATA REQUIRED
 ├── R/
 │    ├── landscape_model.R
 │    ├── kernel.R
 │    ├── initialization.R
 │    ├── within_field_model.R
 │    ├── surveillance_strategies.R
 │    ├── sampling.R
 │    ├── metrics.R
 │    └── utils.R
 ├── scripts/
 │    ├── run_simulation.R
 │    └── analysis.R
 └── outputs/

```
