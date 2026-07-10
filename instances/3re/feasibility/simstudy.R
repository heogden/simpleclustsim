library(tidyverse)
devtools::load_all()

setwd("instances/3re/feasibility/")

source("../design.R")

feasibility_seeds <- 10001:10005

initial_time_each <- 30
extended_time_each <- 240
mem_gb <- 8

setup_feasibility <- tidyr::crossing(
    applicable_method_cells,
    seed = feasibility_seeds
)

## ------------------------------------------------------------------
## Stage 1: create and submit the initial 30-minute Iridis run
## ------------------------------------------------------------------

write_simstudy_iridis(
    setup = setup_feasibility,
    subcases = subcases,
    simulation_fun = run_simstudy_each,
    time_each = initial_time_each,
    mem = mem_gb,
    iridis_dir = "iridis"
)


## ------------------------------------------------------------------
## Stage 2: run after the initial Iridis jobs have completed
## ------------------------------------------------------------------

simstudy_feasibility <- read_simstudy_iridis(
    "iridis"
)

initial_oom_tab <- find_oom_ids_iridis(
    "iridis",
    min_mem_gb = mem_gb,
    mem_limit_if_missing_gb = mem_gb
)

initial_oom_ids <- initial_oom_tab %>%
    dplyr::pull(process_id)

initial_timeout_tab <- find_timeout_ids_iridis(
    "iridis",
    min_time_minutes = initial_time_each,
    time_limit_if_missing_minutes = initial_time_each
)

initial_timeout_ids <- initial_timeout_tab %>%
    dplyr::pull(process_id)

## OOM cases are not rerun: they have already failed with the full
## 8 GB memory allocation used for the feasibility study.
##
## Timeout cases are rerun with a four-hour time limit. The
## write_update_script_iridis() function also includes any other
## missing process IDs, while excluding the known OOM cases.

if (length(initial_timeout_ids) > 0L) {
    update_ids <- write_update_script_iridis(
        simulation_fun = run_simstudy_each,
        time_each = extended_time_each,
        mem = mem_gb,
        iridis_dir = "iridis",
        extra_ids = initial_timeout_ids,
        exclude_ids = initial_oom_ids
    )

    message(
        "Created iridis/simstudy_update.slurm for ",
        length(update_ids),
        " process IDs. Submit this script and wait for it to finish ",
        "before running the final section below."
    )
} else {
    message("No initial timeout cases require a longer rerun.")
}


## ------------------------------------------------------------------
## Stage 3: run after any four-hour update jobs have completed
## ------------------------------------------------------------------

## Re-read the results so that successful extended runs replace the
## previously missing timeout results.
simstudy_feasibility <- read_simstudy_iridis(
    "iridis"
)

## Recompute classifications using the latest log for each process.
## Initial OOM cases were not rerun and therefore remain classified
## as OOM. A timeout is now recorded only if the extended four-hour
## run also timed out.
oom_tab <- find_oom_ids_iridis(
    "iridis",
    min_mem_gb = mem_gb,
    mem_limit_if_missing_gb = mem_gb
)

timeout_tab <- find_timeout_ids_iridis(
    "iridis",
    min_time_minutes = extended_time_each,
    time_limit_if_missing_minutes = extended_time_each
)

oom_ids <- oom_tab %>%
    dplyr::pull(process_id)

timeout_ids <- timeout_tab %>%
    dplyr::pull(process_id)

complete_status <- make_complete_run_status(
    setup = setup_feasibility,
    simstudy = simstudy_feasibility,
    oom_ids = oom_ids,
    timeout_ids = timeout_ids
)

feasibility_summary <- summarise_complete_run_status(
    complete_status
) %>%
    dplyr::left_join(
        subcases,
        by = "subcase_id"
    ) %>%
    dplyr::relocate(
        case,
        subcase_id,
        n_clusters,
        n_obs_per_cluster,
        method
    )

dir.create(
    "output",
    showWarnings = FALSE
)

saveRDS(
    complete_status,
    "output/complete_status_initial.rds"
)

readr::write_csv(
    feasibility_summary,
    "output/feasibility_summary_initial.csv"
)
