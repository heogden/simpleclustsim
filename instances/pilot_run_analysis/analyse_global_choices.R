library(tidyverse)
devtools::load_all()

data_dir <- here::here(
    "instances",
    "pilot_run_analysis",
    "data"
)

output_dir <- here::here(
    "instances",
    "pilot_run_analysis",
    "output"
)

result_files <- c(
    "2re" = file.path(data_dir, "pilot_results_2re.rds"),
    "3re" = file.path(data_dir, "pilot_results_3re.rds"),
    "sitar" = file.path(data_dir, "pilot_results_sitar.rds")
)

pilot_results <- purrr::map(result_files, readRDS)

## GP choice
all_gp_method_cells <- purrr::imap_dfr(
    pilot_results,
    function(result, case_name) {
        result$gp_method_cell_tab %>%
            dplyr::mutate(case = case_name)
    }
)

global_gp_method_summary <-
    make_global_gp_method_comparison_tab(
        all_gp_method_cells
    )


readr::write_csv(
    all_gp_method_cells,
    file.path(
        output_dir,
        "all_gp_method_cells.csv"
    )
)

readr::write_csv(
    global_gp_method_summary,
    file.path(
        output_dir,
        "global_gp_method_summary.csv"
    )
)

## Relative metric, given fixed GP choice

all_relative_metrics <- purrr::imap_dfr(
    pilot_results,
    function(result, case_name) {
        result$relative_metric_tab %>%
            dplyr::mutate(case = case_name)
    }
)

readr::write_csv(
    all_relative_metrics,
    file.path(output_dir, "all_relative_metrics.csv")
)


geometric_mean_positive <- function(x) {
    x <- x[is.finite(x) & x > 0]

    if (length(x) == 0) {
        return(NA_real_)
    }

    exp(mean(log(x)))
}

summarise_variants <- function(data, candidate_methods) {
    data %>%
        dplyr::filter(method %in% candidate_methods) %>%
        dplyr::group_by(method) %>%
        dplyr::summarise(
            n_cells = dplyr::n(),
            n_cells_with_rmw2 =
                sum(!is.na(rel_rmw2_emp_bound)),
            n_cells_with_RMISE =
                sum(!is.na(rel_RMISE)),

            geometric_mean_rel_rmw2 =
                geometric_mean_positive(
                    rel_rmw2_emp_bound
                ),
            mean_rel_rmw2 =
                mean(rel_rmw2_emp_bound, na.rm = TRUE),
            max_rel_rmw2 =
                max(rel_rmw2_emp_bound, na.rm = TRUE),

            geometric_mean_rel_RMISE =
                geometric_mean_positive(rel_RMISE),
            mean_rel_RMISE =
                mean(rel_RMISE, na.rm = TRUE),
            max_rel_RMISE =
                max(rel_RMISE, na.rm = TRUE),

            mean_coverage =
                mean(coverage, na.rm = TRUE),
            mean_CI_width =
                mean(CI_width, na.rm = TRUE),

            mean_fit_fail =
                mean(fit_fail, na.rm = TRUE),
            mean_OOM =
                mean(OOM, na.rm = TRUE),
            mean_Timeout =
                mean(Timeout, na.rm = TRUE),

            mean_time_fit =
                mean(mean_time_fit, na.rm = TRUE),

            .groups = "drop"
        ) %>%
        dplyr::arrange(geometric_mean_rel_rmw2)
}

pace_methods <- c(
    "PACE-BIC-GCV",
    "PACE-AIC-GCV",
    "PACE-FVE95-GCV",
    "PACE-FVE99-GCV",
    "PACE-BIC-default",
    "PACE-AIC-default",
    "PACE-FVE95-default",
    "PACE-FVE99-default"
)

face_methods <- c(
    "FACE-95",
    "FACE-99"
)

bayesfpca_methods <- c(
    "bayesFPCA-Oracle",
    "bayesFPCA-PVE95-L15",
    "bayesFPCA-PVE99-L15"
)

pace_summary <- summarise_variants(
    all_relative_metrics,
    pace_methods
)

face_summary <- summarise_variants(
    all_relative_metrics,
    face_methods
)

bayesfpca_summary <- summarise_variants(
    all_relative_metrics,
    bayesfpca_methods
)

readr::write_csv(
    pace_summary,
    file.path(output_dir, "pace_variant_summary.csv")
)

readr::write_csv(
    face_summary,
    file.path(output_dir, "face_variant_summary.csv")
)

readr::write_csv(
    bayesfpca_summary,
    file.path(
        output_dir,
        "bayesfpca_variant_summary.csv"
    )
)
