ensure_subcase_label <- function(x) {
    if ("subcase_label" %in% names(x)) {
        return(x)
    }

    if (all(c("n_clusters", "n_obs_per_cluster") %in% names(x))) {
        return(
            x %>%
                dplyr::mutate(
                    subcase_label = paste0(
                        "d=", n_clusters,
                        ", n_i=", n_obs_per_cluster
                    )
                )
        )
    }

    stop(
        "Need either subcase_label or both n_clusters and n_obs_per_cluster."
    )
}

make_pilot_performance_tab <- function(summary_tab_agg, complete_status) {
    summary_tab_agg <- summary_tab_agg %>%
        ensure_subcase_label()

    failure_tab <- summarise_complete_run_status(complete_status)

    subcase_info <- summary_tab_agg %>%
        dplyr::distinct(
            case,
            subcase_id,
            n_clusters,
            n_obs_per_cluster,
            subcase_label
        )

    summary_chosen <- summary_tab_agg %>%
        filter_chosen_gp_method() %>%
        dplyr::select(
            -dplyr::any_of(
                c("n_clusters", "n_obs_per_cluster", "subcase_label")
            )
        )

    failure_tab %>%
        dplyr::left_join(
            subcase_info,
            by = c("case", "subcase_id")
        ) %>%
        dplyr::left_join(
            summary_chosen,
            by = c("case", "subcase_id", "method")
        )
}
extract_est <- function(x) {
    purrr::map_dbl(x, function(z) {
        if (is.null(z) || length(z) == 0 || all(is.na(z))) {
            NA_real_
        } else if ("est" %in% names(z)) {
            z[["est"]]
        } else {
            NA_real_
        }
    })
}


make_pilot_metric_table <- function(pilot_performance_tab) {
    pilot_performance_tab <- pilot_performance_tab %>%
        ensure_subcase_label()

    pilot_performance_tab %>%
        dplyr::transmute(
            case,
            subcase_id,
            subcase = subcase_label,
            method,
            GP_method,

            RMISE = extract_est(RMISE),
            coverage = extract_est(coverage),
            CI_width = extract_est(`ci-width`),

            rmw2 = extract_est(rmw2),
            rmw2_emp = extract_est(rmw2_emp),

            n_expected,
            n_read,

            fit_fail = prop_fit_fail,
            pred_fail = prop_pred_fail,
            CI_fail = prop_ci_fail,
            OOM = prop_oom,
            Timeout = prop_timeout,
            unknown_missing = prop_unknown_missing,

            mean_time_total = mean_time_total,
            mean_time_fit = mean_time_fit,
            max_time_fit = max_time_fit
        )
}

add_relative_pilot_metrics <- function(pilot_metric_tab,
                                       individual_oracle_by_case =
                                           tibble::tribble(
                                               ~case, ~individual_oracle_method,
                                               "2re", "Oracle-GP",
                                               "3re", "Oracle-GP",
                                               "sitar", "Oracle-SITAR"
                                           )) {
    oracle_individual <- pilot_metric_tab %>%
        dplyr::inner_join(
            individual_oracle_by_case,
            by = "case"
        ) %>%
        dplyr::filter(method == individual_oracle_method) %>%
        dplyr::select(
            case,
            subcase_id,
            oracle_RMISE = RMISE
        )

    pilot_metric_tab %>%
        dplyr::left_join(
            oracle_individual,
            by = c("case", "subcase_id")
        ) %>%
        dplyr::mutate(
            rel_RMISE = dplyr::if_else(
                !is.na(RMISE) & !is.na(oracle_RMISE) & oracle_RMISE > 0,
                RMISE / oracle_RMISE,
                NA_real_
            ),

            rel_rmw2_emp_bound = dplyr::if_else(
                !is.na(rmw2) & !is.na(rmw2_emp) & rmw2_emp > 0,
                rmw2 / rmw2_emp,
                NA_real_
            )
        )
}

make_relative_method_screening_tab <- function(pilot_metric_tab) {
    rel_tab <- add_relative_pilot_metrics(pilot_metric_tab)

    rel_tab %>%
        dplyr::group_by(case, method) %>%
        dplyr::summarise(
            n_subcases = dplyr::n(),

            mean_rel_RMISE = safe_mean(rel_RMISE),
            min_rel_RMISE = safe_min(rel_RMISE),
            max_rel_RMISE = safe_max(rel_RMISE),

            mean_RMISE = safe_mean(RMISE),
            max_RMISE = safe_max(RMISE),

            mean_coverage = safe_mean(coverage),
            min_coverage = safe_min(coverage),
            mean_abs_coverage_error = safe_mean(abs(coverage - 0.95)),
            mean_CI_width = safe_mean(CI_width),

            mean_rel_rmw2_bound = safe_mean(rel_rmw2_emp_bound),
            min_rel_rmw2_bound = safe_min(rel_rmw2_emp_bound),
            max_rel_rmw2_bound = safe_max(rel_rmw2_emp_bound),
            mean_rmw2 = safe_mean(rmw2),

            fit_fail = safe_mean(fit_fail),
            pred_fail = safe_mean(pred_fail),
            CI_fail = safe_mean(CI_fail),
            OOM = safe_mean(OOM),
            Timeout = safe_mean(Timeout),
            unknown_missing = safe_mean(unknown_missing),
            mean_time = safe_mean(mean_time),

            .groups = "drop"
        ) %>%
        dplyr::arrange(case, mean_rel_RMISE)
}


summarise_bayesFPCA_L <- function(simstudy, true_L = NULL) {
    out <- simstudy %>%
        filter(grepl("^bayesFPCA", method)) %>%
        mutate(L_hat = k_hat)

    if (!is.null(true_L)) {
        out <- out %>%
            mutate(
                selection = case_when(
                    is.na(L_hat) ~ "missing",
                    L_hat < true_L ~ "under",
                    L_hat == true_L ~ "correct",
                    L_hat > true_L ~ "over"
                )
            )
    } else {
        out <- out %>%
            mutate(selection = NA_character_)
    }

    out %>%
        group_by(case, subcase_id, method) %>%
        summarise(
            n = n(),
            mean_L_hat = safe_mean(L_hat),
            min_L_hat = safe_min(L_hat),
            max_L_hat = safe_max(L_hat),
            prop_under = if (is.null(true_L)) NA_real_
                else mean(selection == "under", na.rm = TRUE),
            prop_correct = if (is.null(true_L)) NA_real_
                else mean(selection == "correct", na.rm = TRUE),
            prop_over = if (is.null(true_L)) NA_real_
                else mean(selection == "over", na.rm = TRUE),
            mean_time_fit = safe_mean(time_fit),
            .groups = "drop"
        )
}
