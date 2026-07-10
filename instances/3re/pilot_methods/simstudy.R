library(tidyverse)
devtools::load_all()

setwd("instances/3re/pilot_methods/")

subcases <- tribble(
    ~n_clusters, ~n_obs_per_cluster,
     50,          2,
     50,         10,
    200,          5,
    500,          2,
    500,         10
) %>%
    mutate(subcase_id = row_number())


methods <- c(
    "Oracle-GP",
    "AdaStruMM",

    "PACE-BIC-GCV",
    "PACE-AIC-GCV",
    "PACE-FVE95-GCV",
    "PACE-FVE99-GCV",
    "PACE-BIC-default",
    "PACE-AIC-default",
    "PACE-FVE95-default",
    "PACE-FVE99-default",
    "PACE-Oracle-GCV",
    "PACE-Oracle-default",

    "FACE-95",
    "FACE-99",

    "HGAM-GS",

    "bayesFPCA-Oracle",
    "bayesFPCA-PVE95-L15",
    "bayesFPCA-PVE99-L15",

    "Di-95",
    "Di-99",
    "Di-Oracle",
    "Goldsmith-95",
    "Goldsmith-99",

    "Local-GAM"
)

seeds <- 1:2

setup <- expand_grid(
    case = "3re",
    subcase_id = subcases$subcase_id,
    method = methods,
    seed = seeds
)

## simstudy <- pmap_dfr(setup, run_simstudy_each, subcases = subcases)

## Alternatively, using Iridis:
    
write_simstudy_iridis(setup, subcases, run_simstudy_each,
                      time_each = 30, mem = 8, iridis_dir = "main_run")

simstudy <- read_simstudy_iridis("main_run")

oom_tab <- find_oom_ids_iridis(
    "main_run",
    min_mem_gb = 8,
    mem_limit_if_missing_gb = 8
)

oom_ids <- oom_tab %>%
    dplyr::pull(process_id)

timeout_tab <- find_timeout_ids_iridis(
    "main_run",
    min_time_minutes = 240,
    time_limit_if_missing_minutes = 240
)

timeout_ids <- timeout_tab %>%
    dplyr::pull(process_id)


write_update_script_iridis(run_simstudy_each,
                           time_each = 240,
                           mem = 8,
                           iridis_dir = "main_run",
                           exclude_ids = c(oom_ids, timeout_ids))

simstudy <- read_simstudy_iridis("main_run")

complete_status <- make_complete_run_status(
    setup = setup,
    simstudy = simstudy,
    oom_ids = oom_ids,
    timeout_ids = timeout_ids
)

summary_tab <- summarise_fit(simstudy, emp_error = TRUE)
summary_tab_agg <- find_summary_tab(summary_tab, subcases)

gp_method_cell_tab <- make_gp_method_cell_tab(
    summary_tab_agg
)

gp_method_tab <- make_gp_method_comparison_tab(
    summary_tab_agg
)


pilot_performance_tab <- make_pilot_performance_tab(
    summary_tab_agg,
    complete_status
)

pilot_metric_tab <- make_pilot_metric_table(
    pilot_performance_tab
)

relative_metric_tab <- add_relative_pilot_metrics(pilot_metric_tab)



analysis_data_dir <- "../../pilot_run_analysis/data"

case_name <- "3re"

pilot_results <- list(
    case = case_name,
    saved_at = Sys.time(),

    ## Simulation design
    subcases = subcases,
    setup = setup,

    ## Cell-level performance and relative metrics
    relative_metric_tab = relative_metric_tab,

    ## GP reconstruction comparisons
    gp_method_cell_tab = gp_method_cell_tab,
    gp_method_tab = gp_method_tab,

    ## Useful underlying summaries
    pilot_performance_tab = pilot_performance_tab,
    summary_tab_agg = summary_tab_agg,

    ## Failure and computational information
    complete_status = complete_status
)

saveRDS(
    pilot_results,
    file = file.path(
        analysis_data_dir,
        paste0("pilot_results_", case_name, ".rds")
    )
)

