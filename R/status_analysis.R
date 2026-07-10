safe_mean <- function(x) {
    if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

safe_min <- function(x) {
    if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
}

safe_max <- function(x) {
    if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
}

has_CI_predictions <- function(pred_data_each) {
    if (is.null(pred_data_each)) {
        return(FALSE)
    }

    if (!("mu_c_hat" %in% names(pred_data_each))) {
        return(FALSE)
    }

    out <- tryCatch({
        !all(is.na(pred_data_each$mu_c_hat$lower)) &&
            !all(is.na(pred_data_each$mu_c_hat$upper))
    }, error = function(e) FALSE)

    out
}

add_missing_cols <- function(x, cols, value = NA) {
    missing_cols <- setdiff(cols, names(x))

    for (cc in missing_cols) {
        x[[cc]] <- value
    }

    x
}

make_complete_run_status <- function(
    setup,
    simstudy,
    oom_ids = NULL,
    timeout_ids = NULL
) {
    run_keys <- c("case", "subcase_id", "method", "seed")

    setup2 <- setup

    if (!("process_id" %in% names(setup2))) {
        setup2 <- setup2 %>%
            dplyr::mutate(
                process_id = dplyr::row_number()
            )
    }

    needed_cols <- c(
        run_keys,
        "time", "time_sim", "time_fit", "time_pred", "time_gp",
        "sim_status", "fit_status", "pred_status", "gp_status", "k_status",
        "sim_error", "fit_error", "pred_error", "gp_error", "k_error",
        "pred_data"
    )

    simstudy2 <- simstudy %>%
        add_missing_cols(needed_cols)

    run_status <- simstudy2 %>%
        dplyr::mutate(
            has_CI = purrr::map_lgl(
                pred_data,
                has_CI_predictions
            )
        ) %>%
        dplyr::select(
            dplyr::all_of(run_keys),
            time, time_sim, time_fit, time_pred, time_gp,
            sim_status, fit_status, pred_status, gp_status, k_status,
            sim_error, fit_error, pred_error, gp_error, k_error,
            has_CI
        )

    setup2 %>%
        dplyr::left_join(run_status, by = run_keys) %>%
        dplyr::mutate(
            missing_from_simstudy = is.na(sim_status),

            oom_fail =
                missing_from_simstudy &
                !is.na(process_id) &
                process_id %in% oom_ids,

            timeout_fail =
                missing_from_simstudy &
                !oom_fail &
                !is.na(process_id) &
                process_id %in% timeout_ids,

            unknown_missing =
                missing_from_simstudy &
                !oom_fail &
                !timeout_fail,

            sim_fail =
                !missing_from_simstudy &
                sim_status != "ok",

            fit_fail =
                !missing_from_simstudy &
                sim_status == "ok" &
                fit_status != "ok",

            pred_fail =
                !missing_from_simstudy &
                sim_status == "ok" &
                fit_status == "ok" &
                pred_status != "ok",

            gp_fail =
                !missing_from_simstudy &
                sim_status == "ok" &
                fit_status == "ok" &
                pred_status == "ok" &
                gp_status != "ok",

            k_fail =
                !missing_from_simstudy &
                sim_status == "ok" &
                fit_status == "ok" &
                k_status != "ok",

            ci_fail =
                !missing_from_simstudy &
                sim_status == "ok" &
                fit_status == "ok" &
                pred_status == "ok" &
                !has_CI
        )
}

summarise_complete_run_status <- function(complete_status) {
    complete_status %>%
        dplyr::group_by(case, subcase_id, method) %>%
        dplyr::summarise(
            n_expected = dplyr::n(),
            n_read = sum(!missing_from_simstudy),
            n_missing = sum(missing_from_simstudy),

            n_oom = sum(oom_fail, na.rm = TRUE),
            n_timeout = sum(timeout_fail, na.rm = TRUE),
            n_unknown_missing = sum(
                unknown_missing,
                na.rm = TRUE
            ),

            n_sim_fail = sum(sim_fail, na.rm = TRUE),
            n_fit_fail = sum(fit_fail, na.rm = TRUE),
            n_pred_fail = sum(pred_fail, na.rm = TRUE),
            n_gp_fail = sum(gp_fail, na.rm = TRUE),
            n_k_fail = sum(k_fail, na.rm = TRUE),
            n_ci_fail = sum(ci_fail, na.rm = TRUE),

            prop_oom = n_oom / n_expected,
            prop_timeout = n_timeout / n_expected,
            prop_unknown_missing =
                n_unknown_missing / n_expected,

            prop_sim_fail = n_sim_fail / n_expected,
            prop_fit_fail = n_fit_fail / n_expected,
            prop_pred_fail = n_pred_fail / n_expected,
            prop_gp_fail = n_gp_fail / n_expected,
            prop_k_fail = n_k_fail / n_expected,
            prop_ci_fail = n_ci_fail / n_expected,

            mean_time_total = safe_mean(time),
            mean_time_fit = safe_mean(time_fit),
            max_time_total = safe_max(time),
            max_time_fit = safe_max(time_fit),

            .groups = "drop"
        )
}

summarise_method_failures <- function(complete_status) {
    complete_status %>%
        dplyr::group_by(case, method) %>%
        dplyr::summarise(
            n_expected = dplyr::n(),
            n_read = sum(!missing_from_simstudy),

            prop_oom = mean(oom_fail, na.rm = TRUE),
            prop_timeout = mean(
                timeout_fail,
                na.rm = TRUE
            ),
            prop_unknown_missing = mean(
                unknown_missing,
                na.rm = TRUE
            ),

            prop_fit_fail = mean(
                fit_fail,
                na.rm = TRUE
            ),
            prop_pred_fail = mean(
                pred_fail,
                na.rm = TRUE
            ),
            prop_ci_fail = mean(
                ci_fail,
                na.rm = TRUE
            ),

            mean_time_total = safe_mean(time),
            max_time_total = safe_max(time),

            .groups = "drop"
        )
}
