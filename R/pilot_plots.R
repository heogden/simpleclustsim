format_heatmap_value <- function(x, digits = 2, cap = NULL,
                                 label_capped = TRUE) {
    if (is.null(cap)) {
        return(
            ifelse(
                is.na(x),
                "",
                sprintf(paste0("%.", digits, "f"), x)
            )
        )
    }

    ifelse(
        is.na(x),
        "",
        ifelse(
            label_capped & x > cap,
            paste0(">", sprintf(paste0("%.", digits, "f"), cap)),
            sprintf(paste0("%.", digits, "f"), x)
        )
    )
}

plot_pilot_metric_heatmap <- function(pilot_metric_tab, metric,
                                      label = metric,
                                      digits = 2,
                                      cap = NULL,
                                      fill_min = NULL,
                                      drop_methods = NULL,
                                      label_capped = TRUE,
                                      method_order = NULL,
                                      subcase_order = NULL) {
    plot_dat <- pilot_metric_tab %>%
        dplyr::filter(
            !(method %in% drop_methods)
        ) %>%
        dplyr::mutate(
            value = .data[[metric]],
            value_for_fill = value
        )

    if (!is.null(cap)) {
        plot_dat <- plot_dat %>%
            dplyr::mutate(
                value_for_fill = pmin(value_for_fill, cap)
            )
    }

    if (!is.null(fill_min)) {
        plot_dat <- plot_dat %>%
            dplyr::mutate(
                value_for_fill = pmax(value_for_fill, fill_min)
            )
    }

    plot_dat <- plot_dat %>%
        dplyr::mutate(
            value_label = format_heatmap_value(
                value,
                digits = digits,
                cap = cap,
                label_capped = label_capped
            )
        )

    if (is.null(method_order)) {
        method_order <- plot_dat %>%
            dplyr::distinct(method) %>%
            dplyr::pull(method)
    } else {
        method_order <- setdiff(method_order, drop_methods)
    }

    if (is.null(subcase_order)) {
        subcase_order <- plot_dat %>%
            dplyr::arrange(subcase_id) %>%
            dplyr::distinct(subcase) %>%
            dplyr::pull(subcase)
    }

    p <- plot_dat %>%
        dplyr::mutate(
            method = factor(method, levels = rev(method_order)),
            subcase = factor(subcase, levels = subcase_order)
        ) %>%
        ggplot2::ggplot(
            ggplot2::aes(
                x = subcase,
                y = method,
                fill = value_for_fill
            )
        ) +
        ggplot2::geom_tile() +
        ggplot2::geom_text(
            ggplot2::aes(label = value_label),
            size = 2.5
        ) +
        ggplot2::facet_wrap(~ case) +
        ggplot2::labs(
            x = "Simulation subcase",
            y = "Method",
            fill = label
        )

    if (!is.null(cap)) {
        if (is.null(fill_min)) {
            fill_min <- suppressWarnings(
                min(plot_dat$value_for_fill, na.rm = TRUE)
            )
        }

        breaks <- pretty(c(fill_min, cap), n = 4)
        breaks <- breaks[breaks >= fill_min & breaks <= cap]
        breaks <- sort(unique(c(fill_min, breaks, cap)))

        break_labels <- function(x) {
            out <- sprintf(paste0("%.", digits, "f"), x)
            out[abs(x - cap) < 1e-8] <- paste0(
                "\u2265", sprintf(paste0("%.", digits, "f"), cap)
            )
            out
        }

        p <- p +
            ggplot2::scale_fill_viridis_c(
                limits = c(fill_min, cap),
                breaks = breaks,
                labels = break_labels,
                na.value = "grey80"
            )
    } else {
        p <- p +
            ggplot2::scale_fill_viridis_c(
                na.value = "grey80"
            )
    }

    p
}

plot_gp_method_heatmap <- function(gp_method_tab, digits = 3) {
    gp_method_tab %>%
        ggplot2::ggplot(
            ggplot2::aes(
                x = GP_method,
                y = method,
                fill = mean_rmw2
            )
        ) +
        ggplot2::geom_tile() +
        ggplot2::geom_text(
            ggplot2::aes(
                label = ifelse(
                    chosen,
                    paste0(sprintf(paste0("%.", digits, "f"), mean_rmw2), "*"),
                    sprintf(paste0("%.", digits, "f"), mean_rmw2)
                )
            ),
            size = 2.5
        ) +
        ggplot2::facet_wrap(~ case) +
        ggplot2::labs(
            x = "Population-summary method",
            y = "Fitting method",
            fill = "Mean RMW2",
            caption = "* chosen for main comparison"
        )
}

plot_relative_pilot_heatmaps <- function(relative_metric_tab) {
    list(
        rel_RMISE = plot_pilot_metric_heatmap(
            relative_metric_tab,
            metric = "rel_RMISE",
            label = "RMISE / Oracle-GP",
            digits = 1
        ),

        rel_rmw2 = plot_pilot_metric_heatmap(
            relative_metric_tab,
            metric = "rel_rmw2_emp_bound",
            label = "RMW2 / empirical bound",
            digits = 1
        ),

        coverage = plot_pilot_metric_heatmap(
            relative_metric_tab,
            metric = "coverage",
            label = "Coverage",
            digits = 2
        ),

        CI_width = plot_pilot_metric_heatmap(
            relative_metric_tab,
            metric = "CI_width",
            label = "CI width",
            digits = 2
        )
    )
}


plot_pilot_failure_heatmap <- function(
    pilot_metric_tab,
    digits = 2
) {
    pilot_metric_tab %>%
        dplyr::select(
            case,
            subcase_id,
            subcase,
            method,
            fit_fail,
            pred_fail,
            CI_fail,
            OOM,
            Timeout,
            unknown_missing
        ) %>%
        tidyr::pivot_longer(
            cols = c(
                fit_fail,
                pred_fail,
                CI_fail,
                OOM,
                Timeout,
                unknown_missing
            ),
            names_to = "metric",
            values_to = "value"
        ) %>%
        dplyr::mutate(
            metric = factor(
                metric,
                levels = c(
                    "fit_fail",
                    "pred_fail",
                    "CI_fail",
                    "OOM",
                    "Timeout",
                    "unknown_missing"
                ),
                labels = c(
                    "Fitting failure",
                    "Prediction failure",
                    "CI failure",
                    "Out of memory",
                    "Timeout",
                    "Unknown missing"
                )
            ),

            value_label = format_heatmap_value(
                value,
                digits = digits
            ),

            subcase = factor(
                subcase,
                levels = unique(
                    subcase[order(subcase_id)]
                )
            )
        ) %>%
        ggplot2::ggplot(
            ggplot2::aes(
                x = subcase,
                y = method,
                fill = value
            )
        ) +
        ggplot2::geom_tile() +
        ggplot2::geom_text(
            ggplot2::aes(label = value_label),
            size = 2.3
        ) +
        ggplot2::facet_grid(
            metric ~ case,
            scales = "free_x"
        ) +
        ggplot2::scale_fill_viridis_c(
            limits = c(0, 1),
            na.value = "grey80"
        ) +
        ggplot2::labs(
            x = "Simulation subcase",
            y = "Method",
            fill = "Proportion"
        )
}

plot_pilot_time_heatmap <- function(
    pilot_metric_tab,
    digits = 1,
    cap = NULL
) {
    plot_pilot_metric_heatmap(
        pilot_metric_tab,
        metric = "mean_time",
        label = "Mean completed-run time (minutes)",
        digits = digits,
        cap = cap,
        fill_min = 0
    )
}
