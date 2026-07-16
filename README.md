# simpleclustsim

This repository contains research code and processed simulation
summaries for the simulation study in the paper
["Adaptively-Structured Mixed Models for Simple Clustered Data"](https://arxiv.org/abs/2401.11827).

## Contents

- `R/`: helper functions used to construct, run and analyse the simulation
  studies.
- `instances/2re/`: simulation code and processed results for the two
  random-effect simulation setting.
- `instances/3re/`: simulation code and processed results for the three
  random-effect simulation setting.
- `instances/sitar/`: simulation code and processed results for the growth-curve
  simulation setting.
- `scripts/reproduce_paper_figures.R`: script to regenerate the paper figure
  panels from the saved processed summaries.

The full raw checkpoint outputs from the simulation runs are not included. These
files are large and machine-specific. Instead, the repository includes the saved
summary files used to produce the figures in the paper, especially files in

```text
instances/<case>/analysis/output/
```

including `paper_plot_data.csv` and `oracle_reference_data.csv`.

## Reproducing the paper figure panels

Install the package dependencies from the repository root using

```r
install.packages("remotes")
remotes::install_deps(dependencies = TRUE)
```

Then from the repository root run

```bash
Rscript scripts/reproduce_paper_figures.R
```

This regenerates the TikZ figure files from the saved CSV summaries.

## Simulation scripts

The scripts in `instances/<case>/full/` record the workflow used to run the full
simulation studies. They were written for the computing environment used for the
paper and may need adaptation before being run on another machine or HPC system.

The growth-curve simulations were run in two stages. The main general-purpose
methods were run first. The fitted SITAR comparator was added later using
`instances/sitar/full/simstudy_sitar_method.R`, and incorporated into the
processed summaries. This reflects the workflow used for the paper. The saved
summary files in `instances/sitar/analysis/output/` are the results used to
produce the published figures.

## Additional simulation summaries

The processed CSV summaries used to make the paper figures are included in
`instances/<case>/analysis/output/`.

The preliminary method comparison described in the Supplement was run using
`instances/<case>/pilot_methods/simstudy.R`. Processed per-case results are
included as `instances/<case>/pilot_methods/relative_metric_tab_<case>.csv`.
The combined CSV summaries used to make the global method choices are included
in `instances/pilot_run_analysis/output/`.

Processed summaries from the computational feasibility trial are included in
`instances/<case>/feasibility/output/` and
`instances/feasibility_analysis/output/`. Final run-matrix and run-status CSV
files for the full simulation studies are included in
`instances/<case>/full/output/`.

Raw checkpoint outputs, cluster log files, `.Rda`/`.rds` files and
machine-specific storage directories are not included.
