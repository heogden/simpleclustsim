## Paper-facing simulation-study aggregation.
##
## This file converts run-level measures from analysis.R into tidy
## cell-level data, status summaries and oracle-reference curves.


paper_individual_oracle_method <- function(case_name) {
    switch(
        case_name,
        "2re" = "Oracle-GP",
        "3re" = "Oracle-GP",
        "sitar" = "Oracle-SITAR",
        stop(
            "No individual oracle is defined for case ",
            case_name,
            ".",
            call. = FALSE
        )
    )
}


.paper_metric_summary <- function(
    x,
    root = FALSE,
    level = 0.95
) {
    x <- x[is.finite(x)]
    n_available <- length(x)

    if (n_available == 0L) {
        return(
            tibble::tibble(
                estimate = NA_real_,
                lower = NA_real_,
                upper = NA_real_,
                n_available = 0L
            )
        )
    }

    estimate_raw <- mean(x)

    if (n_available < 2L) {
        lower_raw <- NA_real_
        upper_raw <- NA_real_
    } else {
        se <- stats::sd(x) /
            sqrt(n_available)
        z <- stats::qnorm(
            1 - (1 - level) / 2
        )
        lower_raw <- estimate_raw - z * se
        upper_raw <- estimate_raw + z * se
    }

    if (isTRUE(root)) {
        estimate <- sqrt(
            pmax(estimate_raw, 0)
        )
        lower <- sqrt(
            pmax(lower_raw, 0)
        )
        upper <- sqrt(
            pmax(upper_raw, 0)
        )
    } else {
        estimate <- estimate_raw
        lower <- lower_raw
        upper <- upper_raw
    }

    tibble::tibble(
        estimate = estimate,
        lower = lower,
        upper = upper,
        n_available = n_available
    )
}


.paper_summarise_metric <- function(
    data,
    value_col,
    metric,
    root = FALSE,
    group_cols = c(
        "case",
        "subcase_id",
        "method"
    ),
    level = 0.95
) {
    missing_cols <- setdiff(
        c(group_cols, value_col),
        names(data)
    )

    if (length(missing_cols) > 0L) {
        stop(
            "Metric data are missing: ",
            paste(missing_cols, collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    data |>
        dplyr::group_by(
            dplyr::across(
                dplyr::all_of(group_cols)
            )
        ) |>
        dplyr::summarise(
            metric_summary = list(
                .paper_metric_summary(
                    .data[[value_col]],
                    root = root,
                    level = level
                )
            ),
            .groups = "drop"
        ) |>
        tidyr::unnest(
            metric_summary
        ) |>
        dplyr::mutate(
            metric = metric,
            .before = estimate
        )
}


.paper_first_row_per_run <- function(run_metrics) {
    run_keys <- intersect(
        c(
            "run_id",
            "case",
            "subcase_id",
            "method",
            "seed"
        ),
        names(run_metrics)
    )

    run_metrics |>
        dplyr::group_by(
            dplyr::across(
                dplyr::all_of(run_keys)
            )
        ) |>
        dplyr::slice(1L) |>
        dplyr::ungroup()
}


.paper_add_status_and_design <- function(
    metric_data,
    status_summary,
    subcases,
    ci_failure_threshold
) {
    metric_data |>
        dplyr::left_join(
            status_summary,
            by = c(
                "case",
                "subcase_id",
                "method"
            )
        ) |>
        add_subcase_info(
            subcases
        ) |>
        dplyr::mutate(
            prop_metric_missing =
                dplyr::if_else(
                    is.finite(n_expected) &
                        n_expected > 0,
                    1 -
                        n_available /
                        n_expected,
                    NA_real_
                ),

            show_metric = dplyr::case_when(
                metric %in%
                    c(
                        "coverage",
                        "ci_width"
                    ) ~
                    is.finite(n_expected) &
                    n_expected > 0 &
                    n_available >=
                        ceiling(
                            (1 -
                                ci_failure_threshold) *
                                n_expected
                        ),

                TRUE ~ n_available > 0L
            ),

            method_label =
                method_label(method)
        ) |>
        dplyr::arrange(
            case,
            n_obs_per_cluster,
            n_clusters,
            method,
            metric
        )
}


## Summarise one simulation case from scalar run metrics.
##
## This is the memory-light workhorse.  It expects the output of
## summarise_fit(..., emp_error = TRUE), which can be obtained either from a
## full in-memory simstudy object or by reading checkpoint files one at a time.
summarise_case_simstudy_from_run_metrics <- function(
    run_metrics,
    complete_status,
    subcases,
    case_name,
    ci_failure_threshold = 0.05,
    interval_level = 0.95
) {
    if (
        length(case_name) != 1L ||
        is.na(case_name) ||
        !nzchar(case_name)
    ) {
        stop(
            "case_name must be one non-empty string.",
            call. = FALSE
        )
    }

    if (
        length(ci_failure_threshold) != 1L ||
        !is.finite(ci_failure_threshold) ||
        ci_failure_threshold < 0 ||
        ci_failure_threshold >= 1
    ) {
        stop(
            "ci_failure_threshold must lie in [0, 1).",
            call. = FALSE
        )
    }

    metric_cases <- unique(run_metrics$case)

    if (
        length(metric_cases) != 1L ||
        !identical(
            as.character(metric_cases),
            case_name
        )
    ) {
        stop(
            "run_metrics does not contain only case ",
            case_name,
            ".",
            call. = FALSE
        )
    }

    individual_run_metrics <-
        .paper_first_row_per_run(
            run_metrics
        )

    population_run_metrics <-
        run_metrics |>
        filter_chosen_gp_method() |>
        dplyr::filter(
            !is.na(GP_method),
            GP_method != "Unavailable"
        )

    status_summary <-
        summarise_complete_run_status(
            complete_status
        ) |>
        add_subcase_info(
            subcases
        )

    RMISE_metrics <- individual_run_metrics |>
        dplyr::mutate(
            RMISE_squared = RMISE^2
        ) |>
        .paper_summarise_metric(
            value_col = "RMISE_squared",
            metric = "RMISE",
            root = TRUE,
            level = interval_level
        )

    coverage_metrics <-
        .paper_summarise_metric(
            individual_run_metrics,
            value_col = "coverage",
            metric = "coverage",
            root = FALSE,
            level = interval_level
        )

    ci_width_metrics <-
        .paper_summarise_metric(
            individual_run_metrics,
            value_col = "mean_CI_length",
            metric = "ci_width",
            root = FALSE,
            level = interval_level
        )

    time_metrics <-
        .paper_summarise_metric(
            individual_run_metrics,
            value_col = "time_fit",
            metric = "time",
            root = FALSE,
            level = interval_level
        )

    rmw2_metrics <-
        .paper_summarise_metric(
            population_run_metrics,
            value_col = "W2_bar",
            metric = "rmw2",
            root = TRUE,
            level = interval_level
        )

    cell_metrics <- dplyr::bind_rows(
        RMISE_metrics,
        coverage_metrics,
        ci_width_metrics,
        rmw2_metrics,
        time_metrics
    ) |>
        .paper_add_status_and_design(
            status_summary =
                status_summary,
            subcases = subcases,
            ci_failure_threshold =
                ci_failure_threshold
        )

    oracle_method <-
        paper_individual_oracle_method(
            case_name
        )

    individual_oracle_data <-
        cell_metrics |>
        dplyr::filter(
            method == oracle_method,
            metric %in%
                c(
                    "RMISE",
                    "coverage",
                    "ci_width",
                    "time"
                )
        ) |>
        dplyr::mutate(
            reference_type =
                "individual_oracle",
            reference_label =
                "Individual oracle bound"
        )

    population_oracle_run_metrics <-
        run_metrics |>
        dplyr::select(
            case,
            subcase_id,
            seed,
            W2_bar_emp
        ) |>
        dplyr::distinct()

    population_oracle_expected <-
        complete_status |>
        dplyr::distinct(
            case,
            subcase_id,
            seed
        ) |>
        dplyr::count(
            case,
            subcase_id,
            name = "n_expected"
        )

    population_oracle_data <-
        .paper_summarise_metric(
            population_oracle_run_metrics,
            value_col = "W2_bar_emp",
            metric = "rmw2",
            root = TRUE,
            group_cols = c(
                "case",
                "subcase_id"
            ),
            level = interval_level
        ) |>
        dplyr::left_join(
            population_oracle_expected,
            by = c(
                "case",
                "subcase_id"
            )
        ) |>
        add_subcase_info(
            subcases
        ) |>
        dplyr::mutate(
            method =
                "Population-oracle",
            method_label =
                "Population oracle bound",
            reference_type =
                "population_oracle",
            reference_label =
                "Population oracle bound",
            prop_metric_missing =
                dplyr::if_else(
                    n_expected > 0,
                    1 -
                        n_available /
                        n_expected,
                    NA_real_
                ),
            show_metric =
                n_available > 0L
        )

    method_plot_data <-
        cell_metrics |>
        dplyr::filter(
            method != oracle_method
        )

    oracle_reference_data <-
        dplyr::bind_rows(
            individual_oracle_data,
            population_oracle_data
        ) |>
        dplyr::arrange(
            metric,
            n_obs_per_cluster,
            n_clusters
        )

    metric_availability <-
        cell_metrics |>
        dplyr::select(
            case,
            subcase_id,
            n_clusters,
            n_obs_per_cluster,
            method,
            method_label,
            metric,
            n_expected,
            n_available,
            prop_metric_missing,
            show_metric,
            dplyr::any_of(
                c(
                    "n_oom",
                    "n_timeout",
                    "n_unknown_missing",
                    "n_sim_fail",
                    "n_fit_fail",
                    "n_pred_fail",
                    "n_gp_fail",
                    "n_ci_fail",
                    "prop_oom",
                    "prop_timeout",
                    "prop_unknown_missing",
                    "prop_sim_fail",
                    "prop_fit_fail",
                    "prop_pred_fail",
                    "prop_gp_fail",
                    "prop_ci_fail"
                )
            )
        )

    list(
        case = case_name,
        created_at = Sys.time(),
        settings = list(
            ci_failure_threshold =
                ci_failure_threshold,
            interval_level =
                interval_level,
            individual_oracle_method =
                oracle_method
        ),
        subcases = subcases,
        run_metrics = run_metrics,
        individual_run_metrics =
            individual_run_metrics,
        population_run_metrics =
            population_run_metrics,
        cell_metrics = cell_metrics,
        method_plot_data =
            method_plot_data,
        individual_oracle_data =
            individual_oracle_data,
        population_oracle_data =
            population_oracle_data,
        oracle_reference_data =
            oracle_reference_data,
        status_summary = status_summary,
        metric_availability =
            metric_availability
    )
}


## Summarise one simulation case from an in-memory simstudy object.
summarise_case_simstudy <- function(
    simstudy,
    complete_status,
    subcases,
    case_name,
    ci_failure_threshold = 0.05,
    interval_level = 0.95
) {
    simstudy_cases <- unique(
        simstudy$case
    )

    if (
        length(simstudy_cases) != 1L ||
        !identical(
            as.character(simstudy_cases),
            case_name
        )
    ) {
        stop(
            "simstudy does not contain only case ",
            case_name,
            ".",
            call. = FALSE
        )
    }

    run_metrics <- summarise_fit(
        simstudy,
        emp_error = TRUE
    )

    summarise_case_simstudy_from_run_metrics(
        run_metrics = run_metrics,
        complete_status = complete_status,
        subcases = subcases,
        case_name = case_name,
        ci_failure_threshold = ci_failure_threshold,
        interval_level = interval_level
    )
}


.paper_summarise_one_checkpoint <- function(
    file,
    expected_run_id = NULL
) {
    result <- readRDS(file)

    if (!is.data.frame(result) || nrow(result) != 1L) {
        stop(
            file,
            " does not contain exactly one result row.",
            call. = FALSE
        )
    }

    result <- tibble::as_tibble(result)

    if (!is.null(expected_run_id)) {
        if (!("run_id" %in% names(result))) {
            result$run_id <- expected_run_id
        }

        if (!identical(
            as.character(result$run_id[[1]]),
            as.character(expected_run_id)
        )) {
            stop(
                "run_id inside ",
                file,
                " does not match the checkpoint index.",
                call. = FALSE
            )
        }
    }

    summarise_fit(
        result,
        emp_error = TRUE
    )
}


## Read checkpoint files one at a time and retain only scalar run metrics.
##
## Use the checkpointed_result_index.rds written by the memory-light finalise
## stage.  This avoids constructing the full simstudy object in memory.
summarise_checkpointed_run_metrics <- function(
    checkpoint_index,
    output_file = NULL,
    verbose = TRUE
) {
    if (nrow(checkpoint_index) == 0L) {
        return(tibble::tibble())
    }

    if (isTRUE(verbose)) {
        message(
            "Summarising ",
            nrow(checkpoint_index),
            " checkpointed run files"
        )
    }

    run_metrics <- purrr::map2_dfr(
        checkpoint_index$file,
        checkpoint_index$run_id,
        function(file, run_id) {
            out <- tryCatch(
                .paper_summarise_one_checkpoint(
                    file = file,
                    expected_run_id = run_id
                ),
                error = function(e) {
                    warning(
                        "Could not summarise ",
                        file,
                        ": ",
                        conditionMessage(e),
                        call. = FALSE
                    )
                    NULL
                }
            )

            gc()
            out
        }
    )

    if (!is.null(output_file)) {
        saveRDS(
            run_metrics,
            output_file
        )
    }

    run_metrics
}


## Complete paper summary directly from checkpointed run files.
summarise_case_checkpointed_simstudy <- function(
    checkpoint_index,
    complete_status,
    subcases,
    case_name,
    output_dir = NULL,
    ci_failure_threshold = 0.05,
    interval_level = 0.95,
    verbose = TRUE
) {
    if (!is.null(output_dir)) {
        dir.create(
            output_dir,
            recursive = TRUE,
            showWarnings = FALSE
        )
    }

    run_metrics_file <- if (is.null(output_dir)) {
        NULL
    } else {
        file.path(
            output_dir,
            "run_metrics.rds"
        )
    }

    run_metrics <- summarise_checkpointed_run_metrics(
        checkpoint_index = checkpoint_index,
        output_file = run_metrics_file,
        verbose = verbose
    )

    analysis_results <- summarise_case_simstudy_from_run_metrics(
        run_metrics = run_metrics,
        complete_status = complete_status,
        subcases = subcases,
        case_name = case_name,
        ci_failure_threshold = ci_failure_threshold,
        interval_level = interval_level
    )

    if (!is.null(output_dir)) {
        write_case_simstudy_analysis(
            analysis_results,
            output_dir = output_dir
        )
    }

    analysis_results
}

## Write the analysis object and its scalar tables.
write_case_simstudy_analysis <- function(
    analysis_results,
    output_dir
) {
    dir.create(
        output_dir,
        recursive = TRUE,
        showWarnings = FALSE
    )

    saveRDS(
        analysis_results,
        file.path(
            output_dir,
            "analysis_results.rds"
        )
    )

    saveRDS(
        analysis_results$run_metrics,
        file.path(
            output_dir,
            "run_metrics.rds"
        )
    )

    readr::write_csv(
        analysis_results$cell_metrics,
        file.path(
            output_dir,
            "cell_metrics.csv"
        )
    )

    readr::write_csv(
        analysis_results$method_plot_data,
        file.path(
            output_dir,
            "paper_plot_data.csv"
        )
    )

    readr::write_csv(
        analysis_results$oracle_reference_data,
        file.path(
            output_dir,
            "oracle_reference_data.csv"
        )
    )

    readr::write_csv(
        analysis_results$status_summary,
        file.path(
            output_dir,
            "status_summary.csv"
        )
    )

    readr::write_csv(
        analysis_results$metric_availability,
        file.path(
            output_dir,
            "metric_availability.csv"
        )
    )

    invisible(analysis_results)
}
