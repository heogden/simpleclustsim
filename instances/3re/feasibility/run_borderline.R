library(tidyverse)
devtools::load_all()

setwd("instances/3re/feasibility/")

source("../design.R")

borderline_setup_file <- "output/borderline_setup.csv"
borderline_iridis_dir <- "iridis_borderline"

initial_time_each <- 30
extended_time_each <- 240
mem_gb <- 8


## ------------------------------------------------------------------
## Read and validate the additional borderline runs
## ------------------------------------------------------------------

if (!file.exists(borderline_setup_file)) {
    stop(
        "Cannot find ",
        borderline_setup_file,
        ". Create this file from the initial feasibility decisions first."
    )
}

setup_borderline <- readr::read_csv(
    borderline_setup_file,
    show_col_types = FALSE
) %>%
    dplyr::select(
        case,
        subcase_id,
        method,
        seed
    ) %>%
    dplyr::distinct() %>%
    dplyr::arrange(
        case,
        subcase_id,
        method,
        seed
    ) %>%
    dplyr::mutate(
        process_id = dplyr::row_number()
    )

if (nrow(setup_borderline) == 0L) {
    stop(
        "There are no borderline runs for the 3re case.",
        call. = FALSE
    )
}

if (any(setup_borderline$case != case_name)) {
    stop(
        "borderline_setup.csv contains a case other than ",
        case_name,
        ".",
        call. = FALSE
    )
}

invalid_cells <- setup_borderline %>%
    dplyr::distinct(
        case,
        subcase_id,
        method
    ) %>%
    dplyr::anti_join(
        applicable_method_cells,
        by = c(
            "case",
            "subcase_id",
            "method"
        )
    )

if (nrow(invalid_cells) > 0L) {
    print(invalid_cells)

    stop(
        "borderline_setup.csv contains method--subcase cells ",
        "which are not listed in applicable_method_cells.",
        call. = FALSE
    )
}


## ------------------------------------------------------------------
## Stage 1: create and submit the initial 30-minute Iridis run
## ------------------------------------------------------------------

write_simstudy_iridis(
    setup = setup_borderline,
    subcases = subcases,
    simulation_fun = run_simstudy_each,
    time_each = initial_time_each,
    mem = mem_gb,
    iridis_dir = borderline_iridis_dir
)


## ------------------------------------------------------------------
## Stage 2: run after the initial Iridis jobs have completed
## ------------------------------------------------------------------

simstudy_borderline <- read_simstudy_iridis(
    borderline_iridis_dir
)

initial_oom_tab <- find_oom_ids_iridis(
    borderline_iridis_dir,
    min_mem_gb = mem_gb,
    mem_limit_if_missing_gb = mem_gb
)

initial_oom_ids <- initial_oom_tab %>%
    dplyr::pull(process_id)

initial_timeout_tab <- find_timeout_ids_iridis(
    borderline_iridis_dir,
    min_time_minutes = initial_time_each,
    time_limit_if_missing_minutes = initial_time_each
)

initial_timeout_ids <- initial_timeout_tab %>%
    dplyr::pull(process_id)

## OOM cases are not rerun because they have already failed with the
## full 8 GB memory allocation used for the feasibility study.
##
## Timeout cases are rerun with a four-hour limit. Specifying them in
## extra_ids ensures that they are rerun even if they left an output
## file. Known OOM cases are explicitly excluded.

if (length(initial_timeout_ids) > 0L) {
    update_ids <- write_update_script_iridis(
        simulation_fun = run_simstudy_each,
        time_each = extended_time_each,
        mem = mem_gb,
        iridis_dir = borderline_iridis_dir,
        extra_ids = initial_timeout_ids,
        exclude_ids = initial_oom_ids
    )

    message(
        "Created ",
        file.path(
            borderline_iridis_dir,
            "simstudy_update.slurm"
        ),
        " for ",
        length(update_ids),
        " process IDs. Submit this script and wait for it to finish ",
        "before running Stage 3."
    )
} else {
    message(
        "No borderline runs timed out at 30 minutes; ",
        "no extended-time update is required."
    )
}


## ------------------------------------------------------------------
## Stage 3: run after any four-hour update jobs have completed
## ------------------------------------------------------------------

## Re-read the results. Successful extended runs now replace the
## previously missing 30-minute timeout results.
simstudy_borderline <- read_simstudy_iridis(
    borderline_iridis_dir
)

## Recompute the final resource-failure classifications using the
## latest Slurm log for each process. A run is now a timeout only if
## its four-hour rerun also timed out.
oom_tab <- find_oom_ids_iridis(
    borderline_iridis_dir,
    min_mem_gb = mem_gb,
    mem_limit_if_missing_gb = mem_gb
)

timeout_tab <- find_timeout_ids_iridis(
    borderline_iridis_dir,
    min_time_minutes = extended_time_each,
    time_limit_if_missing_minutes = extended_time_each
)

oom_ids <- oom_tab %>%
    dplyr::pull(process_id)

timeout_ids <- timeout_tab %>%
    dplyr::pull(process_id)

complete_status_borderline <- make_complete_run_status(
    setup = setup_borderline,
    simstudy = simstudy_borderline,
    oom_ids = oom_ids,
    timeout_ids = timeout_ids
)

feasibility_cell_summary_borderline <-
    summarise_feasibility_cells(
        complete_status_borderline
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
    complete_status_borderline,
    "output/complete_status_borderline.rds"
)

readr::write_csv(
    feasibility_cell_summary_borderline,
    "output/feasibility_cell_summary_borderline.csv"
)

readr::write_csv(
    oom_tab,
    "output/oom_borderline.csv"
)

readr::write_csv(
    timeout_tab,
    "output/timeout_borderline.csv"
)
