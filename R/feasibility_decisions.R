## Helpers for the feasibility stage of the simulation studies.
##
## Inclusion is based on model-fitting time and reliability.  Total elapsed
## time is used only to estimate the computational cost of the full study and
## to help allocate Iridis jobs.

.feas_check_cols <- function(x, required, object_name = deparse(substitute(x))) {
    missing <- setdiff(required, names(x))

    if (length(missing) > 0L) {
        stop(
            object_name,
            " is missing required column",
            if (length(missing) == 1L) " " else "s ",
            paste(missing, collapse = ", "),
            call. = FALSE
        )
    }

    invisible(x)
}

.feas_safe_mean <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else mean(x)
}

.feas_safe_median <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else stats::median(x)
}

.feas_safe_sd <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 2L) NA_real_ else stats::sd(x)
}

.feas_safe_max <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else max(x)
}

.feas_max_to_median <- function(x) {
    x <- x[is.finite(x)]

    if (length(x) == 0L) {
        return(NA_real_)
    }

    med <- stats::median(x)

    if (!is.finite(med) || med <= 0) {
        return(NA_real_)
    }

    max(x) / med
}


## Summarise run-level feasibility results by case, subcase and method.
##
## complete_status should be the output from make_complete_run_status().
## A run is counted as complete only if all stages required for the full
## simulation study succeed.  k_status is reported but is not required for
## completion because K is not meaningful for every method.
summarise_feasibility_cells <- function(
    complete_status,
    require_gp = TRUE,
    require_ci = TRUE
) {
    required <- c(
        "case", "subcase_id", "method", "seed",
        "missing_from_simstudy", "oom_fail", "timeout_fail",
        "unknown_missing", "sim_fail", "fit_fail", "pred_fail",
        "gp_fail", "k_fail", "ci_fail", "has_CI",
        "sim_status", "fit_status", "pred_status", "gp_status",
        "time", "time_fit"
    )

    .feas_check_cols(complete_status, required)

    run_status <- complete_status |>
        dplyr::mutate(
            resource_fail =
                dplyr::coalesce(oom_fail, FALSE) |
                dplyr::coalesce(timeout_fail, FALSE),

            fit_success =
                !dplyr::coalesce(missing_from_simstudy, TRUE) &
                sim_status == "ok" &
                fit_status == "ok" &
                is.finite(time_fit),

            pipeline_complete =
                !dplyr::coalesce(missing_from_simstudy, TRUE) &
                sim_status == "ok" &
                fit_status == "ok" &
                pred_status == "ok" &
                (!require_gp | gp_status == "ok") &
                (!require_ci | dplyr::coalesce(has_CI, FALSE)),

            nonresource_incomplete =
                !pipeline_complete &
                !resource_fail &
                !dplyr::coalesce(unknown_missing, FALSE),

            successful_fit_time = dplyr::if_else(
                fit_success,
                time_fit,
                NA_real_
            ),

            observed_total_time = dplyr::if_else(
                !dplyr::coalesce(missing_from_simstudy, TRUE) &
                    is.finite(time),
                time,
                NA_real_
            ),

            complete_total_time = dplyr::if_else(
                pipeline_complete & is.finite(time),
                time,
                NA_real_
            )
        )

    run_status |>
        dplyr::group_by(case, subcase_id, method) |>
        dplyr::summarise(
            n_trials = dplyr::n(),
            n_read = sum(!missing_from_simstudy, na.rm = TRUE),
            n_complete = sum(pipeline_complete, na.rm = TRUE),

            n_oom = sum(oom_fail, na.rm = TRUE),
            n_timeout = sum(timeout_fail, na.rm = TRUE),
            n_resource_fail = sum(resource_fail, na.rm = TRUE),
            n_unknown_missing = sum(unknown_missing, na.rm = TRUE),

            n_sim_fail = sum(sim_fail, na.rm = TRUE),
            n_fit_fail = sum(fit_fail, na.rm = TRUE),
            n_pred_fail = sum(pred_fail, na.rm = TRUE),
            n_gp_fail = sum(gp_fail, na.rm = TRUE),
            n_k_fail = sum(k_fail, na.rm = TRUE),
            n_ci_fail = sum(ci_fail, na.rm = TRUE),
            n_nonresource_incomplete = sum(
                nonresource_incomplete,
                na.rm = TRUE
            ),

            n_successful_fits = sum(fit_success, na.rm = TRUE),
            prop_complete = n_complete / n_trials,
            prop_resource_fail = n_resource_fail / n_trials,

            mean_fit_time = .feas_safe_mean(successful_fit_time),
            median_fit_time = .feas_safe_median(successful_fit_time),
            sd_fit_time = .feas_safe_sd(successful_fit_time),
            max_fit_time = .feas_safe_max(successful_fit_time),
            max_to_median_fit_time = .feas_max_to_median(
                successful_fit_time
            ),

            mean_total_time_observed = .feas_safe_mean(
                observed_total_time
            ),
            max_total_time_observed = .feas_safe_max(
                observed_total_time
            ),
            mean_total_time_complete = .feas_safe_mean(
                complete_total_time
            ),
            max_total_time_complete = .feas_safe_max(
                complete_total_time
            ),

            .groups = "drop"
        )
}


## Apply the pre-specified feasibility rules.
##
## Initial stage (normally five seeds):
##   * clear include/exclude decisions are made immediately;
##   * ambiguous cells are labelled borderline for five additional seeds.
##
## Final stage (normally ten seeds for initially borderline cells):
##   * include if at least 80% complete, at most one resource failure,
##     and mean fitting time is no greater than the common cutoff;
##   * exclude for clear time, resource or reliability problems;
##   * otherwise request manual review.
classify_feasibility_cells <- function(
    feasibility_summary,
    fit_cutoff_minutes = 30,
    fit_time_exempt_methods = character(),
    stage = c("initial", "final"),
    borderline_fraction = 0.20,
    min_successful_fits_for_time = 3L
) {
    stage <- match.arg(stage)

    if (
        length(fit_cutoff_minutes) != 1L ||
        !is.finite(fit_cutoff_minutes) ||
        fit_cutoff_minutes <= 0
    ) {
        stop(
            "fit_cutoff_minutes must be one positive number.",
            call. = FALSE
        )
    }

    if (
        length(borderline_fraction) != 1L ||
        !is.finite(borderline_fraction) ||
        borderline_fraction < 0 ||
        borderline_fraction >= 1
    ) {
        stop(
            "borderline_fraction must lie in [0, 1).",
            call. = FALSE
        )
    }

    fit_time_exempt_methods <- unique(
        as.character(fit_time_exempt_methods)
    )

    if (
        anyNA(fit_time_exempt_methods) ||
        any(!nzchar(fit_time_exempt_methods))
    ) {
        stop(
            "fit_time_exempt_methods must contain non-empty method names.",
            call. = FALSE
        )
    }

    required <- c(
        "case", "subcase_id", "method", "n_trials", "n_complete",
        "n_resource_fail", "n_unknown_missing",
        "n_nonresource_incomplete", "n_successful_fits",
        "mean_fit_time", "max_fit_time"
    )

    .feas_check_cols(feasibility_summary, required)

    lower_fit_limit <- (1 - borderline_fraction) * fit_cutoff_minutes
    upper_fit_limit <- (1 + borderline_fraction) * fit_cutoff_minutes

    out <- feasibility_summary |>
        dplyr::mutate(
            fit_cutoff_minutes = fit_cutoff_minutes,
            fit_time_exempt =
                method %in% fit_time_exempt_methods,
            decision_stage = stage
        )

    if (stage == "initial") {
        out <- out |>
            dplyr::mutate(
                decision = dplyr::case_when(
                    n_resource_fail >= 2L ~ "exclude_resource",

                    n_complete <= floor(0.4 * n_trials) ~
                        "exclude_unreliable",

                    !fit_time_exempt &
                        n_successful_fits >= pmin(
                            min_successful_fits_for_time,
                            n_trials
                        ) &
                        is.finite(mean_fit_time) &
                        mean_fit_time > upper_fit_limit ~
                        "exclude_time",

                    n_complete == n_trials &
                        n_resource_fail == 0L &
                        n_unknown_missing == 0L &
                        n_nonresource_incomplete == 0L &
                        (
                            fit_time_exempt |
                                (
                                    is.finite(mean_fit_time) &
                                        mean_fit_time < lower_fit_limit &
                                        (
                                            !is.finite(max_fit_time) |
                                                max_fit_time <=
                                                    fit_cutoff_minutes
                                        )
                                )
                        ) ~
                        "include",

                    TRUE ~ "borderline"
                ),

                decision_reason = dplyr::case_when(
                    decision == "exclude_resource" ~
                        "At least two OOMs or four-hour timeouts",

                    decision == "exclude_unreliable" ~
                        "At most 40% of runs completed",

                    decision == "exclude_time" ~
                        "Mean fitting time clearly above the cutoff",

                    decision == "include" & fit_time_exempt ~
                        paste0(
                            "All runs completed; method exempt from the ",
                            "fitting-time cutoff"
                        ),

                    decision == "include" ~
                        paste0(
                            "All runs completed; mean fitting time below ",
                            format(lower_fit_limit, trim = TRUE),
                            " minutes and no fitting time exceeded the cutoff"
                        ),

                    n_resource_fail == 1L ~
                        "One OOM or four-hour timeout",

                    n_unknown_missing > 0L ~
                        "One or more unexplained missing runs",

                    n_nonresource_incomplete > 0L ~
                        "One or more non-resource pipeline failures",

                    !fit_time_exempt &
                        is.finite(mean_fit_time) &
                        mean_fit_time >= lower_fit_limit &
                        mean_fit_time <= upper_fit_limit ~
                        "Mean fitting time close to the cutoff",

                    !fit_time_exempt &
                        is.finite(max_fit_time) &
                        max_fit_time > fit_cutoff_minutes ~
                        "At least one fitting time exceeded the cutoff",

                    TRUE ~ "Additional runs required"
                )
            )
    } else {
        out <- out |>
            dplyr::mutate(
                min_complete_for_include = ceiling(0.8 * n_trials),
                min_complete_for_review = ceiling(0.7 * n_trials),
                max_resource_fail_for_include = floor(0.1 * n_trials),

                decision = dplyr::case_when(
                    n_resource_fail >= 2L ~ "exclude_resource",

                    !fit_time_exempt &
                        n_successful_fits >= pmin(
                            min_successful_fits_for_time,
                            n_trials
                        ) &
                        is.finite(mean_fit_time) &
                        mean_fit_time > fit_cutoff_minutes ~
                        "exclude_time",

                    n_complete < min_complete_for_review ~
                        "exclude_unreliable",

                    n_complete >= min_complete_for_include &
                        n_resource_fail <=
                            max_resource_fail_for_include &
                        n_unknown_missing == 0L &
                        (
                            fit_time_exempt |
                                (
                                    is.finite(mean_fit_time) &
                                        mean_fit_time <=
                                            fit_cutoff_minutes
                                )
                        ) ~
                        "include",

                    TRUE ~ "manual_review"
                ),

                decision_reason = dplyr::case_when(
                    decision == "exclude_resource" ~
                        "At least two OOMs or four-hour timeouts",

                    decision == "exclude_time" ~
                        "Mean fitting time above the common cutoff",

                    decision == "exclude_unreliable" ~
                        "Fewer than 70% of runs completed",

                    decision == "include" & fit_time_exempt ~
                        paste0(
                            "At least 80% of runs completed; method exempt ",
                            "from the fitting-time cutoff"
                        ),

                    decision == "include" ~
                        paste0(
                            "At least 80% of runs completed and mean ",
                            "fitting time did not exceed ",
                            format(fit_cutoff_minutes, trim = TRUE),
                            " minutes"
                        ),

                    n_unknown_missing > 0L ~
                        "One or more unexplained missing runs",

                    n_complete >= min_complete_for_review &
                        n_complete < min_complete_for_include ~
                        "Between 70% and 80% of runs completed",

                    n_resource_fail >
                        max_resource_fail_for_include ~
                        "Resource-failure rate requires review",

                    !fit_time_exempt &
                        !is.finite(mean_fit_time) ~
                        "Insufficient successful fits to estimate fitting time",

                    TRUE ~ "Manual review required"
                )
            ) |>
            dplyr::select(
                -min_complete_for_include,
                -min_complete_for_review,
                -max_resource_fail_for_include
            )
    }

    out
}

## Construct the extra-seed setup for cells classified as borderline.
make_borderline_setup <- function(
    decisions,
    extra_seeds = 10006:10010,
    decision_values = "borderline"
) {
    required <- c("case", "subcase_id", "method", "decision")
    .feas_check_cols(decisions, required)

    if (
        length(extra_seeds) == 0L ||
        anyNA(extra_seeds) ||
        anyDuplicated(extra_seeds)
    ) {
        stop(
            "extra_seeds must be a non-empty vector of distinct values.",
            call. = FALSE
        )
    }

    cells <- decisions |>
        dplyr::filter(decision %in% decision_values) |>
        dplyr::distinct(case, subcase_id, method)

    tidyr::crossing(
        cells,
        seed = extra_seeds
    )
}


## Replace the initially borderline decisions by the decisions based on all
## ten runs.  Initial non-borderline decisions are left unchanged.
finalise_feasibility_decisions <- function(
    initial_decisions,
    borderline_final_decisions
) {
    keys <- c("case", "subcase_id", "method")
    required <- c(keys, "decision", "decision_reason")

    .feas_check_cols(
        initial_decisions,
        required,
        "initial_decisions"
    )
    .feas_check_cols(
        borderline_final_decisions,
        required,
        "borderline_final_decisions"
    )

    initially_borderline <- initial_decisions |>
        dplyr::filter(decision == "borderline") |>
        dplyr::select(dplyr::all_of(keys))

    missing_final <- initially_borderline |>
        dplyr::anti_join(
            borderline_final_decisions |>
                dplyr::select(dplyr::all_of(keys)),
            by = keys
        )

    if (nrow(missing_final) > 0L) {
        stop(
            "Final decisions are missing for ",
            nrow(missing_final),
            " initially borderline cells.",
            call. = FALSE
        )
    }

    initial_decisions |>
        dplyr::filter(decision != "borderline") |>
        dplyr::bind_rows(borderline_final_decisions) |>
        dplyr::arrange(case, subcase_id, method)
}


## Estimate the total computational cost of the retained cells.
##
## The inclusion decision uses fitting time, but cost uses total observed
## elapsed time because simulation, prediction and GP calculations consume
## real Iridis resources.  The scheduling buffer is reported separately and
## does not alter the inclusion decision.
estimate_full_simulation_cost <- function(
    feasibility_summary,
    decisions,
    n_full_seeds,
    included_decisions = "include",
    scheduling_buffer = 1.25
) {
    keys <- c("case", "subcase_id", "method")

    .feas_check_cols(
        feasibility_summary,
        c(
            keys,
            "mean_total_time_observed",
            "max_total_time_observed"
        )
    )
    .feas_check_cols(decisions, c(keys, "decision"))

    if (
        length(n_full_seeds) != 1L ||
        !is.finite(n_full_seeds) ||
        n_full_seeds <= 0 ||
        n_full_seeds != as.integer(n_full_seeds)
    ) {
        stop(
            "n_full_seeds must be one positive integer.",
            call. = FALSE
        )
    }

    if (
        length(scheduling_buffer) != 1L ||
        !is.finite(scheduling_buffer) ||
        scheduling_buffer < 1
    ) {
        stop(
            "scheduling_buffer must be at least 1.",
            call. = FALSE
        )
    }

    by_cell <- feasibility_summary |>
        dplyr::inner_join(
            decisions |>
                dplyr::select(dplyr::all_of(keys), decision),
            by = keys
        ) |>
        dplyr::filter(decision %in% included_decisions) |>
        dplyr::mutate(
            n_full_seeds = as.integer(n_full_seeds),
            projected_total_hours =
                n_full_seeds * mean_total_time_observed / 60,
            projected_total_hours_buffered =
                scheduling_buffer * projected_total_hours,
            suggested_walltime_minutes =
                scheduling_buffer * max_total_time_observed
        )

    missing_time <- by_cell |>
        dplyr::filter(
            !is.finite(projected_total_hours) |
                !is.finite(suggested_walltime_minutes)
        )

    if (nrow(missing_time) > 0L) {
        warning(
            "Cost or walltime could not be estimated for ",
            nrow(missing_time),
            " retained cells.",
            call. = FALSE
        )
    }

    by_case <- by_cell |>
        dplyr::group_by(case) |>
        dplyr::summarise(
            n_cells = dplyr::n(),
            projected_total_hours = sum(
                projected_total_hours,
                na.rm = TRUE
            ),
            projected_total_hours_buffered = sum(
                projected_total_hours_buffered,
                na.rm = TRUE
            ),
            .groups = "drop"
        )

    overall <- by_cell |>
        dplyr::summarise(
            n_cells = dplyr::n(),
            projected_total_hours = sum(
                projected_total_hours,
                na.rm = TRUE
            ),
            projected_total_hours_buffered = sum(
                projected_total_hours_buffered,
                na.rm = TRUE
            )
        )

    list(
        by_cell = by_cell,
        by_case = by_case,
        overall = overall
    )
}


## Compare alternative common fitting-time cutoffs before fixing the final
## cutoff for all three simulation cases.
##
## At the initial stage, include and borderline cells are both counted in the
## provisional cost estimate.  At the final stage, only included cells are
## counted.
compare_feasibility_cutoffs <- function(
    feasibility_summary,
    cutoffs = c(20, 25, 30),
    n_full_seeds,
    fit_time_exempt_methods = character(),
    stage = c("initial", "final"),
    scheduling_buffer = 1.25,
    borderline_fraction = 0.20
) {
    stage <- match.arg(stage)

    if (
        length(cutoffs) == 0L ||
        anyNA(cutoffs) ||
        any(!is.finite(cutoffs)) ||
        any(cutoffs <= 0)
    ) {
        stop(
            "cutoffs must contain positive finite values.",
            call. = FALSE
        )
    }

    included_decisions <- if (stage == "initial") {
        c("include", "borderline")
    } else {
        "include"
    }

    results <- lapply(cutoffs, function(cutoff) {
        decisions <- classify_feasibility_cells(
            feasibility_summary = feasibility_summary,
            fit_cutoff_minutes = cutoff,
            fit_time_exempt_methods = fit_time_exempt_methods,
            stage = stage,
            borderline_fraction = borderline_fraction
        )

        cost <- estimate_full_simulation_cost(
            feasibility_summary = feasibility_summary,
            decisions = decisions,
            n_full_seeds = n_full_seeds,
            included_decisions = included_decisions,
            scheduling_buffer = scheduling_buffer
        )

        tibble::tibble(
            fit_cutoff_minutes = cutoff,
            n_cells = nrow(decisions),
            n_include = sum(decisions$decision == "include"),
            n_borderline = sum(decisions$decision == "borderline"),
            n_manual_review = sum(
                decisions$decision == "manual_review"
            ),
            n_exclude_time = sum(
                decisions$decision == "exclude_time"
            ),
            n_exclude_resource = sum(
                decisions$decision == "exclude_resource"
            ),
            n_exclude_unreliable = sum(
                decisions$decision == "exclude_unreliable"
            ),
            n_cells_in_cost_estimate = cost$overall$n_cells,
            projected_total_hours =
                cost$overall$projected_total_hours,
            projected_total_hours_buffered =
                cost$overall$projected_total_hours_buffered
        )
    })

    dplyr::bind_rows(results)
}
