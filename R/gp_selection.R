choose_gp_method <- function(method) {
    dplyr::case_when(
        method %in% c("Oracle-GP", "Oracle-SITAR") ~ "Individual",

        method == "AdaStruMM" ~ "EmpiricalCorrected",

        grepl("^PACE", method) ~ "FPC",
        grepl("^FACE", method) ~ "FPC",

        method %in% c("HGAM-GS", "SITAR") ~ "Empirical",
        grepl("^bayesFPCA", method) ~ "FPC",
        grepl("^Di", method) ~ "Empirical",
        grepl("^Goldsmith", method) ~ "Empirical",
        grepl("^Local", method) ~ "Empirical",

        TRUE ~ NA_character_
    )
}


add_chosen_gp_method <- function(summary_tab) {
    summary_tab %>%
        dplyr::mutate(
            chosen_GP_method = choose_gp_method(method),
            chosen = !is.na(chosen_GP_method) &
                GP_method == chosen_GP_method
        )
}


filter_chosen_gp_method <- function(summary_tab) {
    summary_tab %>%
        add_chosen_gp_method() %>%
        dplyr::filter(
            is.na(chosen_GP_method) | chosen
        )
}


extract_gp_est <- function(x) {
    purrr::map_dbl(
        x,
        function(z) {
            if (is.null(z) ||
                length(z) == 0 ||
                all(is.na(z))) {
                return(NA_real_)
            }

            if ("est" %in% names(z)) {
                return(as.numeric(z[["est"]]))
            }

            NA_real_
        }
    )
}


gp_geometric_mean_positive <- function(x) {
    x <- x[is.finite(x) & x > 0]

    if (length(x) == 0) {
        return(NA_real_)
    }

    exp(mean(log(x)))
}


gp_safe_mean <- function(x) {
    x <- x[is.finite(x)]

    if (length(x) == 0) {
        return(NA_real_)
    }

    mean(x)
}


gp_safe_max <- function(x) {
    x <- x[is.finite(x)]

    if (length(x) == 0) {
        return(NA_real_)
    }

    max(x)
}


gp_safe_min <- function(x) {
    x <- x[is.finite(x)]

    if (length(x) == 0) {
        return(NA_real_)
    }

    min(x)
}


#' Construct cell-level GP-method comparison data
#'
#' Population error is first made relative to the empirical population
#' bound within each simulation cell. Each GP reconstruction is also
#' compared with the best available reconstruction for the same
#' case, subcase and fitted method.
make_gp_method_cell_tab <- function(summary_tab_agg) {
    summary_tab_agg %>%
        dplyr::filter(
            !is.na(GP_method),
            GP_method != "Unavailable"
        ) %>%
        dplyr::mutate(
            rmw2_est = extract_gp_est(rmw2),
            rmw2_emp_est = extract_gp_est(rmw2_emp),

            rel_rmw2_emp_bound = dplyr::if_else(
                is.finite(rmw2_est) &
                    is.finite(rmw2_emp_est) &
                    rmw2_emp_est > 0,
                rmw2_est / rmw2_emp_est,
                NA_real_
            )
        ) %>%
        dplyr::group_by(
            case,
            subcase_id,
            method
        ) %>%
        dplyr::mutate(
            best_rel_rmw2_emp_bound =
                gp_safe_min(rel_rmw2_emp_bound),

            rel_to_best_gp_method = dplyr::if_else(
                is.finite(rel_rmw2_emp_bound) &
                    is.finite(best_rel_rmw2_emp_bound) &
                    best_rel_rmw2_emp_bound > 0,
                rel_rmw2_emp_bound /
                    best_rel_rmw2_emp_bound,
                NA_real_
            )
        ) %>%
        dplyr::ungroup() %>%
        add_chosen_gp_method() %>%
        dplyr::select(
            case,
            subcase_id,
            dplyr::any_of(
                c(
                    "n_clusters",
                    "n_obs_per_cluster",
                    "subcase_label"
                )
            ),
            method,
            GP_method,
            chosen_GP_method,
            chosen,
            rmw2_est,
            rmw2_emp_est,
            rel_rmw2_emp_bound,
            best_rel_rmw2_emp_bound,
            rel_to_best_gp_method
        )
}


#' Summarise GP-method comparisons within each simulation case
make_gp_method_comparison_tab <- function(summary_tab_agg) {
    make_gp_method_cell_tab(summary_tab_agg) %>%
        dplyr::group_by(
            case,
            method,
            GP_method,
            chosen_GP_method,
            chosen
        ) %>%
        dplyr::summarise(
            n_subcases =
                sum(is.finite(rel_rmw2_emp_bound)),

            geometric_mean_rel_rmw2_emp_bound =
                gp_geometric_mean_positive(
                    rel_rmw2_emp_bound
                ),

            mean_rel_rmw2_emp_bound =
                gp_safe_mean(rel_rmw2_emp_bound),

            max_rel_rmw2_emp_bound =
                gp_safe_max(rel_rmw2_emp_bound),

            geometric_mean_rel_to_best =
                gp_geometric_mean_positive(
                    rel_to_best_gp_method
                ),

            mean_rel_to_best =
                gp_safe_mean(rel_to_best_gp_method),

            max_rel_to_best =
                gp_safe_max(rel_to_best_gp_method),

            n_best = sum(
                rel_to_best_gp_method <= 1 + 1e-8,
                na.rm = TRUE
            ),

            prop_within_05 = mean(
                rel_to_best_gp_method <= 1.05,
                na.rm = TRUE
            ),

            prop_within_10 = mean(
                rel_to_best_gp_method <= 1.10,
                na.rm = TRUE
            ),

            .groups = "drop"
        ) %>%
        dplyr::arrange(
            case,
            method,
            geometric_mean_rel_rmw2_emp_bound
        )
}


#' Summarise GP-method comparisons over all simulation cases
make_global_gp_method_comparison_tab <- function(
    all_gp_method_cell_tab
) {
    all_gp_method_cell_tab %>%
        dplyr::group_by(
            method,
            GP_method,
            chosen_GP_method,
            chosen
        ) %>%
        dplyr::summarise(
            n_cells =
                sum(is.finite(rel_rmw2_emp_bound)),

            n_cases = dplyr::n_distinct(
                case[is.finite(rel_rmw2_emp_bound)]
            ),

            geometric_mean_rel_rmw2_emp_bound =
                gp_geometric_mean_positive(
                    rel_rmw2_emp_bound
                ),

            mean_rel_rmw2_emp_bound =
                gp_safe_mean(rel_rmw2_emp_bound),

            max_rel_rmw2_emp_bound =
                gp_safe_max(rel_rmw2_emp_bound),

            geometric_mean_rel_to_best =
                gp_geometric_mean_positive(
                    rel_to_best_gp_method
                ),

            mean_rel_to_best =
                gp_safe_mean(rel_to_best_gp_method),

            max_rel_to_best =
                gp_safe_max(rel_to_best_gp_method),

            n_best = sum(
                rel_to_best_gp_method <= 1 + 1e-8,
                na.rm = TRUE
            ),

            prop_within_05 = mean(
                rel_to_best_gp_method <= 1.05,
                na.rm = TRUE
            ),

            prop_within_10 = mean(
                rel_to_best_gp_method <= 1.10,
                na.rm = TRUE
            ),

            .groups = "drop"
        ) %>%
        dplyr::arrange(
            method,
            geometric_mean_rel_rmw2_emp_bound
        )
}
