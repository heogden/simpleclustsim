## Run-level simulation-study summaries.
##
## This file contains only the calculations which turn one fitted
## simulation run into scalar performance measures. Paper-specific
## aggregation and plotting are in simstudy_paper_summary.R and
## simstudy_paper_plots.R.


.analysis_first_non_null <- function(x) {
    keep <- !purrr::map_lgl(
        x,
        function(z) {
            is.null(z) || length(z) == 0L
        }
    )

    if (!any(keep)) {
        return(NULL)
    }

    x[[which(keep)[1L]]]
}


.analysis_has_prediction_estimates <- function(pred_data_each) {
    if (
        is.null(pred_data_each) ||
        !is.data.frame(pred_data_each) ||
        !("mu_c_hat" %in% names(pred_data_each))
    ) {
        return(FALSE)
    }

    out <- tryCatch(
        !all(is.na(pred_data_each$mu_c_hat$estimate)),
        error = function(e) FALSE
    )

    isTRUE(out)
}


.analysis_has_prediction_intervals <- function(pred_data_each) {
    if (!.analysis_has_prediction_estimates(pred_data_each)) {
        return(FALSE)
    }

    out <- tryCatch(
        !all(is.na(pred_data_each$mu_c_hat$lower)) &&
            !all(is.na(pred_data_each$mu_c_hat$upper)),
        error = function(e) FALSE
    )

    isTRUE(out)
}


add_subcase_info <- function(x, subcases) {
    if (
        all(
            c(
                "n_clusters",
                "n_obs_per_cluster"
            ) %in% names(x)
        )
    ) {
        out <- x
    } else {
        if (is.null(subcases)) {
            stop(
                "subcases must be supplied because x does not contain ",
                "n_clusters and n_obs_per_cluster.",
                call. = FALSE
            )
        }

        join_cols <- if (
            "case" %in% names(subcases) &&
            "case" %in% names(x)
        ) {
            c("case", "subcase_id")
        } else {
            "subcase_id"
        }

        out <- dplyr::left_join(
            x,
            subcases,
            by = join_cols
        )
    }

    out |>
        dplyr::mutate(
            subcase_label = paste0(
                "d=", n_clusters,
                ", n_i=", n_obs_per_cluster
            )
        )
}


find_RMISE_each <- function(pred_data_each) {
    if (!.analysis_has_prediction_estimates(pred_data_each)) {
        return(NA_real_)
    }

    x_grid <- sort(unique(pred_data_each$x))

    if (length(x_grid) < 2L) {
        return(NA_real_)
    }

    h <- x_grid[2L] - x_grid[1L]

    result <- pred_data_each |>
        dplyr::group_by(c) |>
        dplyr::mutate(
            SE = (mu_c_hat$estimate - mu_c)^2
        ) |>
        dplyr::summarise(
            ISE = sum(h * SE),
            .groups = "drop"
        ) |>
        dplyr::summarise(
            RMISE = sqrt(mean(ISE, na.rm = TRUE))
        )

    as.numeric(result$RMISE)
}


find_RMedISE_each <- function(pred_data_each) {
    if (!.analysis_has_prediction_estimates(pred_data_each)) {
        return(NA_real_)
    }

    x_grid <- sort(unique(pred_data_each$x))

    if (length(x_grid) < 2L) {
        return(NA_real_)
    }

    h <- x_grid[2L] - x_grid[1L]

    result <- pred_data_each |>
        dplyr::group_by(c) |>
        dplyr::mutate(
            SE = (mu_c_hat$estimate - mu_c)^2
        ) |>
        dplyr::summarise(
            ISE = sum(h * SE),
            .groups = "drop"
        ) |>
        dplyr::summarise(
            RMedISE = sqrt(
                stats::median(
                    ISE,
                    na.rm = TRUE
                )
            )
        )

    as.numeric(result$RMedISE)
}


find_coverage_each <- function(pred_data_each) {
    if (!.analysis_has_prediction_intervals(pred_data_each)) {
        return(NA_real_)
    }

    result <- pred_data_each |>
        dplyr::ungroup() |>
        dplyr::mutate(
            covers =
                mu_c_hat$lower < mu_c &
                mu_c_hat$upper > mu_c
        ) |>
        dplyr::summarise(
            coverage = mean(
                covers,
                na.rm = TRUE
            )
        )

    as.numeric(result$coverage)
}


find_mean_CI_length_each <- function(pred_data_each) {
    if (!.analysis_has_prediction_intervals(pred_data_each)) {
        return(NA_real_)
    }

    result <- pred_data_each |>
        dplyr::ungroup() |>
        dplyr::mutate(
            CI_length =
                mu_c_hat$upper -
                mu_c_hat$lower
        ) |>
        dplyr::summarise(
            mean_CI_length = mean(
                CI_length,
                na.rm = TRUE
            )
        )

    as.numeric(result$mean_CI_length)
}


find_emp_error <- function(simstudy) {
    required <- c(
        "case",
        "subcase_id",
        "seed",
        "GP0",
        "GP0_emp"
    )

    missing_cols <- setdiff(
        required,
        names(simstudy)
    )

    if (length(missing_cols) > 0L) {
        stop(
            "simstudy is missing: ",
            paste(missing_cols, collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    simstudy |>
        dplyr::group_by(
            case,
            subcase_id,
            seed
        ) |>
        dplyr::summarise(
            GP0 = list(
                .analysis_first_non_null(GP0)
            ),
            GP0_emp = list(
                .analysis_first_non_null(GP0_emp)
            ),
            .groups = "drop"
        ) |>
        dplyr::mutate(
            d_m_bar_emp = purrr::map2_dbl(
                GP0,
                GP0_emp,
                find_d_m_bar_GP
            ),
            d_C_bar_emp = purrr::map2_dbl(
                GP0,
                GP0_emp,
                find_d_C_bar_GP
            ),
            W2_bar_emp =
                d_m_bar_emp +
                d_C_bar_emp
        ) |>
        dplyr::select(
            case,
            subcase_id,
            seed,
            d_m_bar_emp,
            d_C_bar_emp,
            W2_bar_emp
        )
}


## Return one row per run and available GP reconstruction.
summarise_fit <- function(simstudy, emp_error = FALSE) {
    required <- c(
        "case",
        "subcase_id",
        "method",
        "seed",
        "pred_data",
        "GP0",
        "GP0_emp",
        "GP_hat",
        "k_hat",
        "time",
        "time_fit"
    )

    missing_cols <- setdiff(
        required,
        names(simstudy)
    )

    if (length(missing_cols) > 0L) {
        stop(
            "simstudy is missing: ",
            paste(missing_cols, collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    id_cols <- intersect(
        c(
            "run_id",
            "process_id",
            "case",
            "subcase_id",
            "method",
            "seed"
        ),
        names(simstudy)
    )

    simstudy_long <- simstudy |>
        tidyr::unnest_longer(
            GP_hat,
            values_to = "GP_hat_value",
            indices_to = "GP_method",
            keep_empty = TRUE
        )

    simstudy_summary <- simstudy_long |>
        dplyr::mutate(
            GP0_use = purrr::pmap(
                list(
                    case,
                    GP0,
                    GP0_emp
                ),
                function(
                    case_each,
                    GP0_each,
                    GP0_emp_each
                ) {
                    if (case_each %in% c("fat")) {
                        GP0_emp_each
                    } else {
                        GP0_each
                    }
                }
            ),
            RMISE = purrr::map_dbl(
                pred_data,
                find_RMISE_each
            ),
            RMedISE = purrr::map_dbl(
                pred_data,
                find_RMedISE_each
            ),
            coverage = purrr::map_dbl(
                pred_data,
                find_coverage_each
            ),
            mean_CI_length = purrr::map_dbl(
                pred_data,
                find_mean_CI_length_each
            ),
            d_m_bar = purrr::map2_dbl(
                GP0_use,
                GP_hat_value,
                find_d_m_bar_GP
            ),
            d_C_bar = purrr::map2_dbl(
                GP0_use,
                GP_hat_value,
                find_d_C_bar_GP
            ),
            W2_bar =
                d_m_bar +
                d_C_bar
        ) |>
        dplyr::select(
            dplyr::all_of(id_cols),
            GP_method,
            RMISE,
            RMedISE,
            coverage,
            mean_CI_length,
            d_m_bar,
            d_C_bar,
            W2_bar,
            k_hat,
            time,
            time_fit
        )

    if (!isTRUE(emp_error)) {
        return(simstudy_summary)
    }

    empirical_error <- find_emp_error(
        simstudy
    )

    dplyr::left_join(
        simstudy_summary,
        empirical_error,
        by = c(
            "case",
            "subcase_id",
            "seed"
        )
    )
}


mean_and_interval <- function(x, level = 0.95) {
    x <- x[is.finite(x)]

    if (length(x) == 0L) {
        return(
            list(
                list(
                    est = NA_real_,
                    lower = NA_real_,
                    upper = NA_real_
                )
            )
        )
    }

    estimate <- mean(x)

    if (length(x) < 2L) {
        lower <- NA_real_
        upper <- NA_real_
    } else {
        se <- stats::sd(x) / sqrt(length(x))
        z <- stats::qnorm(
            1 - (1 - level) / 2
        )
        lower <- estimate - z * se
        upper <- estimate + z * se
    }

    list(
        list(
            est = estimate,
            lower = lower,
            upper = upper
        )
    )
}


root_mean_and_interval <- function(x, level = 0.95) {
    out <- mean_and_interval(
        x,
        level = level
    )[[1L]]

    out$est <- sqrt(
        pmax(out$est, 0)
    )
    out$lower <- sqrt(
        pmax(out$lower, 0)
    )
    out$upper <- sqrt(
        pmax(out$upper, 0)
    )

    list(out)
}


## Legacy aggregate table retained for the pilot-analysis scripts.
find_summary_tab <- function(summary, subcases = NULL) {
    summary |>
        dplyr::group_by(
            case,
            method,
            GP_method,
            subcase_id
        ) |>
        dplyr::summarise(
            RMISE =
                root_mean_and_interval(
                    RMISE^2
                ),
            RMedISE =
                root_mean_and_interval(
                    RMedISE^2
                ),
            coverage =
                mean_and_interval(
                    coverage
                ),
            `ci-width` =
                mean_and_interval(
                    mean_CI_length
                ),
            rmdm =
                root_mean_and_interval(
                    d_m_bar
                ),
            rmdC =
                root_mean_and_interval(
                    d_C_bar
                ),
            rmw2 =
                root_mean_and_interval(
                    W2_bar
                ),
            rmw2_emp =
                root_mean_and_interval(
                    W2_bar_emp
                ),
            rmdm_emp =
                root_mean_and_interval(
                    d_m_bar_emp
                ),
            rmdC_emp =
                root_mean_and_interval(
                    d_C_bar_emp
                ),
            k_hat =
                mean_and_interval(
                    k_hat
                ),
            time =
                mean_and_interval(
                    time_fit
                ),
            .groups = "drop"
        ) |>
        add_subcase_info(
            subcases
        ) |>
        dplyr::arrange(
            case,
            n_clusters,
            n_obs_per_cluster,
            method,
            GP_method
        )
}


