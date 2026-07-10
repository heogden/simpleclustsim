## Paper plots for the full simulation studies.
##
## The named colour and shape vectors ensure that each method keeps the
## same appearance across cases.

paper_method_order <- function() {
    c(
        "AdaStruMM",
        "Local-GAM",
        "HGAM-GS",
        "PACE",
        "FACE",
        "bayesFPCA",
        "SITAR"
    )
}


paper_method_colours <- function() {
    c(
        "AdaStruMM" = "#D55E00",
        "Local-GAM" = "#E69F00",
        "HGAM-GS" = "#009E73",
        "PACE" = "#0072B2",
        "FACE" = "#CC79A7",
        "bayesFPCA" = "#56B4E9",
        "SITAR" = "#666666"
    )
}

paper_method_shapes <- function() {
    c(
        "AdaStruMM" = 16,
        "Local-GAM" = 17,
        "HGAM-GS" = 15,
        "PACE" = 18,
        "FACE" = 3,
        "bayesFPCA" = 7,
        "SITAR" = 4
    )
}


paper_method_key <- function() {
    tibble::tibble(
        method = paper_method_order(),
        method_label = method_label(method),
        colour = unname(
            paper_method_colours()[method]
        ),
        shape = unname(
            paper_method_shapes()[method]
        )
    )
}


paper_metric_y_label <- function(
    metric,
    tex = TRUE
) {
    if (isTRUE(tex)) {
        switch(
            metric,
            "RMISE" = "RMISE",
            "coverage" = "Coverage",
            "ci_width" =
                "Average CI width",
            "rmw2" =
                "$\\mathrm{RM}\\overline{W}_2$",
            "time" =
                "Fitting time (minutes)",
            metric
        )
    } else {
        switch(
            metric,
            "RMISE" = "RMISE",
            "coverage" = "Coverage",
            "ci_width" =
                "Average CI width",
            "rmw2" =
                "Root mean W2 error",
            "time" =
                "Fitting time (minutes)",
            metric
        )
    }
}


paper_method_display_labels <- function() {
    labels <- method_labels()

    unname(
        labels[
            paper_method_order()
        ]
    ) |>
        stats::setNames(
            paper_method_order()
        )
}


theme_simstudy_paper <- function(
    base_size = 9.5
) {
    ggplot2::theme_bw(
        base_size = base_size
    ) +
        ggplot2::theme(
            panel.grid.minor =
                ggplot2::element_blank(),
            strip.background =
                ggplot2::element_blank(),
            strip.placement = "outside",
            legend.title =
                ggplot2::element_blank(),
            legend.key.width =
                grid::unit(1.1, "lines"),
            plot.margin =
                ggplot2::margin(
                    2,
                    3,
                    2,
                    2,
                    unit = "pt"
                )
        )
}


.paper_facet_labeller <- function(
    tex = TRUE
) {
    if (isTRUE(tex)) {
        ggplot2::as_labeller(
            function(x) {
                paste0(
                    "$n_i = ",
                    x,
                    "$"
                )
            }
        )
    } else {
        ggplot2::as_labeller(
            function(x) {
                paste0(
                    "n_i = ",
                    x
                )
            }
        )
    }
}


.paper_blank_facet_labeller <- function() {
    ggplot2::as_labeller(
        function(x) {
            rep(
                "",
                length(x)
            )
        }
    )
}

.paper_apply_plot_display_rules <- function(
    method_data,
    metric,
    local_gam_min_ni_for_error = 5,
    local_gam_show_ci = FALSE
) {
    if (nrow(method_data) == 0L) {
        return(
            list(
                method_data = method_data,
                exclusions = tibble::tibble()
            )
        )
    }

    method_data <- method_data |>
        dplyr::mutate(
                   plot_suppression_reason =
                       dplyr::case_when(
                                  method == "Local-GAM" &
                                  metric %in% c("RMISE", "rmw2") &
                                  n_obs_per_cluster < local_gam_min_ni_for_error &
                                  dplyr::coalesce(show_metric, TRUE) ~
                                      paste0(
                                          "Local-GAM suppressed for ",
                                          "n_i < ",
                                          local_gam_min_ni_for_error
                                      ),

                                  method == "Local-GAM" &
                                  metric %in% c("coverage", "ci_width") &
                                  !isTRUE(local_gam_show_ci) &
                                  dplyr::coalesce(show_metric, TRUE) ~
                                      "Local-GAM CI results suppressed",

                                  method == "SITAR" &
                                  metric %in% c("RMISE", "rmw2", "time") &
                                  n_obs_per_cluster < 5 &
                                  dplyr::coalesce(show_metric, TRUE) ~
                                      "SITAR suppressed for n_i < 5 because fitting was unreliable",

                                  TRUE ~ NA_character_
                              ),
            show_metric = dplyr::if_else(
                !is.na(plot_suppression_reason),
                FALSE,
                dplyr::coalesce(show_metric, TRUE)
            )
        )

    exclusions <- method_data |>
        dplyr::filter(
            !is.na(plot_suppression_reason)
        ) |>
        dplyr::mutate(
            metric = .env$metric,
            reason = plot_suppression_reason
        )

    if (nrow(exclusions) > 0L) {
        if (!("case" %in% names(exclusions))) {
            exclusions$case <- NA_character_
        }

        if (!("subcase_id" %in% names(exclusions))) {
            exclusions$subcase_id <- NA_integer_
        }

        exclusions <- exclusions |>
            dplyr::select(
                dplyr::any_of(
                    c(
                        "case",
                        "subcase_id",
                        "metric",
                        "method",
                        "n_obs_per_cluster",
                        "n_clusters",
                        "estimate",
                        "lower",
                        "upper",
                        "reason"
                    )
                )
            ) |>
            dplyr::arrange(
                metric,
                method,
                n_obs_per_cluster,
                n_clusters
            )
    } else {
        exclusions <- tibble::tibble()
    }

    method_data <- method_data |>
        dplyr::select(
            -plot_suppression_reason
        )

    list(
        method_data = method_data,
        exclusions = exclusions
    )
}


.paper_warn_plot_exclusions <- function(
    exclusions,
    metric
) {
    if (
        is.null(exclusions) ||
        nrow(exclusions) == 0L
    ) {
        return(invisible(exclusions))
    }

    exclusion_summary <- exclusions |>
        dplyr::count(
            method,
            reason,
            name = "n_points"
        ) |>
        dplyr::arrange(
            method,
            reason
        ) |>
        dplyr::mutate(
            text = paste0(
                method,
                " (",
                reason,
                ": ",
                n_points,
                ")"
            )
        )

    warning(
        "For metric ",
        metric,
        ", suppressed points using paper plotting rules: ",
        paste(
            exclusion_summary$text,
            collapse = "; "
        ),
        call. = FALSE
    )

    invisible(exclusions)
}


.paper_write_plot_exclusions <- function(
    exclusions,
    plot_exclusions_file = NULL
) {
    if (
        is.null(plot_exclusions_file) ||
        is.null(exclusions) ||
        nrow(exclusions) == 0L
    ) {
        return(invisible(exclusions))
    }

    dir.create(
        dirname(plot_exclusions_file),
        recursive = TRUE,
        showWarnings = FALSE
    )

    file_exists <- file.exists(plot_exclusions_file)

    readr::write_csv(
        exclusions,
        file = plot_exclusions_file,
        append = file_exists,
        col_names = !file_exists
    )

    invisible(exclusions)
}



.paper_prepare_plot_data <- function(
    analysis_results,
    metric,
    apply_plot_rules = TRUE,
    local_gam_min_ni_for_error = 5,
    local_gam_show_ci = FALSE
) {
    method_data <-
        analysis_results$method_plot_data |>
        dplyr::filter(
            .data$metric == .env$metric
        )

    exclusions <- tibble::tibble()

    if (isTRUE(apply_plot_rules)) {
        display_rule_results <-
            .paper_apply_plot_display_rules(
                method_data = method_data,
                metric = metric,
                local_gam_min_ni_for_error =
                    local_gam_min_ni_for_error,
                local_gam_show_ci =
                    local_gam_show_ci
            )

        method_data <-
            display_rule_results$method_data
        exclusions <-
            display_rule_results$exclusions
    }

    method_data <-
        method_data |>
        dplyr::mutate(
            method = factor(
                method,
                levels =
                    paper_method_order()
            ),
            estimate = dplyr::if_else(
                dplyr::coalesce(show_metric, FALSE),
                estimate,
                NA_real_
            ),
            lower = dplyr::if_else(
                dplyr::coalesce(show_metric, FALSE),
                lower,
                NA_real_
            ),
            upper = dplyr::if_else(
                dplyr::coalesce(show_metric, FALSE),
                upper,
                NA_real_
            )
        )

    oracle_data <-
        analysis_results$oracle_reference_data |>
        dplyr::filter(
            .data$metric == .env$metric
        ) |>
        dplyr::mutate(
            estimate = dplyr::if_else(
                dplyr::coalesce(show_metric, FALSE),
                estimate,
                NA_real_
            )
        )

    list(
        method_data = method_data,
        oracle_data = oracle_data,
        exclusions = exclusions
    )
}


.paper_metric_uses_oracle_limit <- function(metric) {
    metric %in% c(
        "RMISE",
        "ci_width",
        "rmw2"
    )
}


.paper_add_plot_values <- function(data) {
    data |>
        dplyr::mutate(
            estimate_original = estimate,
            lower_original = lower,
            upper_original = upper,
            estimate_plot = estimate,
            lower_plot = lower,
            upper_plot = upper,
            estimate_omitted_for_plot = FALSE,
            interval_clipped_for_plot = FALSE
        )
}


.paper_oracle_y_upper <- function(
    oracle_data,
    metric,
    oracle_y_multiplier
) {
    if (
        !.paper_metric_uses_oracle_limit(metric) ||
        !is.finite(oracle_y_multiplier) ||
        oracle_y_multiplier <= 0 ||
        nrow(oracle_data) == 0L
    ) {
        return(NA_real_)
    }

    oracle_estimates <- oracle_data$estimate[
        is.finite(oracle_data$estimate) &
            dplyr::coalesce(
                oracle_data$show_metric,
                TRUE
            )
    ]

    if (length(oracle_estimates) == 0L) {
        return(NA_real_)
    }

    y_upper <- oracle_y_multiplier *
        max(
            oracle_estimates,
            na.rm = TRUE
        )

    if (!is.finite(y_upper) || y_upper <= 0) {
        return(NA_real_)
    }

    y_upper
}


.paper_needs_oracle_y_limit <- function(
    method_data,
    oracle_data,
    y_upper
) {
    if (!is.finite(y_upper)) {
        return(FALSE)
    }

    y_values <- c(
        method_data$estimate_original,
        method_data$lower_original,
        method_data$upper_original,
        oracle_data$estimate_original,
        oracle_data$lower_original,
        oracle_data$upper_original
    )

    y_values <- y_values[
        is.finite(y_values)
    ]

    if (length(y_values) == 0L) {
        return(FALSE)
    }

    max(
        y_values,
        na.rm = TRUE
    ) > y_upper
}


.paper_apply_oracle_y_limit <- function(
    method_data,
    y_upper
) {
    if (!is.finite(y_upper)) {
        return(method_data)
    }

    ## Do not filter or clip the data here.  The plot uses
    ## coord_cartesian(), so values outside the y-range are retained in
    ## the ggplot object and only clipped visually.  This means that
    ## line segments can still visibly run upwards out of the panel.
    method_data |>
        dplyr::mutate(
            estimate_omitted_for_plot =
                is.finite(estimate_original) &
                estimate_original > y_upper,

            interval_clipped_for_plot =
                is.finite(upper_original) &
                upper_original > y_upper,

            estimate_plot = estimate_original,
            lower_plot = lower_original,
            upper_plot = upper_original
        )
}


.paper_warn_y_omissions <- function(
    method_data,
    metric,
    y_upper
) {
    if (!is.finite(y_upper)) {
        return(invisible(method_data))
    }

    affected <- method_data |>
        dplyr::filter(
            estimate_omitted_for_plot |
                interval_clipped_for_plot
        )

    if (nrow(affected) == 0L) {
        return(invisible(method_data))
    }

    affected_summary <- affected |>
        dplyr::mutate(
            action = dplyr::case_when(
                estimate_omitted_for_plot ~
                    "estimate above plotting range",
                interval_clipped_for_plot ~
                    "interval extends above plotting range",
                TRUE ~ "affected"
            ),
            method_label = method_label(
                as.character(method)
            )
        ) |>
        dplyr::count(
            method_label,
            action,
            name = "n_points"
        ) |>
        dplyr::arrange(
            method_label,
            action
        ) |>
        dplyr::mutate(
            text = paste0(
                method_label,
                " (",
                action,
                ": ",
                n_points,
                ")"
            )
        )

    warning(
        "For metric ",
        metric,
        ", the y-axis upper limit is ",
        signif(y_upper, 4),
        ". Some values extend beyond the plotted y-range for: ",
        paste(
            affected_summary$text,
            collapse = "; "
        ),
        call. = FALSE
    )

    invisible(method_data)
}


.paper_coverage_y_limits <- function(
    method_data,
    oracle_data,
    coverage_ylim = c(0.7, 1),
    coverage_lower_padding = 0.02
) {
    if (!is.null(coverage_ylim)) {
        if (
            !is.numeric(coverage_ylim) ||
            length(coverage_ylim) != 2L ||
            any(!is.finite(coverage_ylim)) ||
            coverage_ylim[1] >= coverage_ylim[2]
        ) {
            stop(
                "coverage_ylim must be NULL or a numeric vector ",
                "of length 2 with increasing finite values.",
                call. = FALSE
            )
        }

        return(coverage_ylim)
    }

    coverage_values <- c(
        method_data$lower_plot,
        method_data$estimate_plot,
        oracle_data$estimate_plot
    )

    coverage_values <- coverage_values[
        is.finite(coverage_values)
    ]

    if (length(coverage_values) == 0L) {
        return(NULL)
    }

    lower <- min(
        coverage_values,
        na.rm = TRUE
    ) - coverage_lower_padding

    lower <- max(
        0,
        lower
    )

    c(
        lower,
        1
    )
}


.paper_warn_coverage_y_limits <- function(
    method_data,
    oracle_data,
    coverage_y_limits
) {
    if (
        is.null(coverage_y_limits) ||
        length(coverage_y_limits) != 2L
    ) {
        return(invisible(NULL))
    }

    lower <- coverage_y_limits[[1]]
    upper <- coverage_y_limits[[2]]

    coverage_values <- c(
        method_data$lower_plot,
        method_data$estimate_plot,
        method_data$upper_plot,
        oracle_data$estimate_plot
    )

    coverage_values <- coverage_values[
        is.finite(coverage_values)
    ]

    if (length(coverage_values) == 0L) {
        return(invisible(NULL))
    }

    if (
        any(coverage_values < lower) ||
        any(coverage_values > upper)
    ) {
        warning(
            "Some displayed coverage values extend outside ",
            "coverage_ylim = c(",
            paste(
                signif(coverage_y_limits, 4),
                collapse = ", "
            ),
            "). They are retained in the plot object and clipped ",
            "visually by coord_cartesian().",
            call. = FALSE
        )
    }

    invisible(NULL)
}


## Make one paper panel.
plot_simstudy_metric <- function(
    analysis_results,
    metric = c(
        "RMISE",
        "coverage",
        "ci_width",
        "rmw2",
        "time"
    ),
    show_interval = TRUE,
    show_oracle = TRUE,
    show_nominal_coverage = TRUE,
    legend = FALSE,
    tex = TRUE,
    log_time = FALSE,
    oracle_y_multiplier = 5,
    clip_to_oracle_multiple = TRUE,
    coverage_ylim = c(0.7, 1),
    coverage_lower_padding = 0.02,
    base_size = 9.5,
    show_facet_labels = NULL,
    show_axis_labels = FALSE,
    apply_plot_rules = TRUE,
    local_gam_min_ni_for_error = 5,
    local_gam_show_ci = FALSE,
    plot_exclusions_file = NULL
) {
    metric <- match.arg(metric)

    if (is.null(show_facet_labels)) {
        show_facet_labels <- identical(
            metric,
            "RMISE"
        )
    }

    plot_data <-
        .paper_prepare_plot_data(
            analysis_results = analysis_results,
            metric = metric,
            apply_plot_rules =
                apply_plot_rules,
            local_gam_min_ni_for_error =
                local_gam_min_ni_for_error,
            local_gam_show_ci =
                local_gam_show_ci
        )

    .paper_warn_plot_exclusions(
        exclusions = plot_data$exclusions,
        metric = metric
    )

    .paper_write_plot_exclusions(
        exclusions = plot_data$exclusions,
        plot_exclusions_file =
            plot_exclusions_file
    )

    method_data <- plot_data$method_data |>
        .paper_add_plot_values()

    oracle_data <- plot_data$oracle_data |>
        .paper_add_plot_values()

    y_upper <- NA_real_

    if (
        isTRUE(clip_to_oracle_multiple) &&
        .paper_metric_uses_oracle_limit(metric)
    ) {
        candidate_y_upper <- .paper_oracle_y_upper(
            oracle_data = oracle_data,
            metric = metric,
            oracle_y_multiplier =
                oracle_y_multiplier
        )

        if (
            .paper_needs_oracle_y_limit(
                method_data = method_data,
                oracle_data = oracle_data,
                y_upper = candidate_y_upper
            )
        ) {
            y_upper <- candidate_y_upper

            method_data <- .paper_apply_oracle_y_limit(
                method_data = method_data,
                y_upper = y_upper
            )

            .paper_warn_y_omissions(
                method_data = method_data,
                metric = metric,
                y_upper = y_upper
            )
        }
    }

    coverage_y_limits <- NULL

    if (identical(metric, "coverage")) {
        coverage_y_limits <- .paper_coverage_y_limits(
            method_data = method_data,
            oracle_data = oracle_data,
            coverage_ylim =
                coverage_ylim,
            coverage_lower_padding =
                coverage_lower_padding
        )

        .paper_warn_coverage_y_limits(
            method_data = method_data,
            oracle_data = oracle_data,
            coverage_y_limits =
                coverage_y_limits
        )
    }

    if (nrow(method_data) == 0L) {
        stop(
            "There are no ordinary-method data for metric ",
            metric,
            ".",
            call. = FALSE
        )
    }

    visible_method_data <- method_data |>
        dplyr::filter(
            is.finite(estimate_plot) |
                is.finite(lower_plot) |
                is.finite(upper_plot)
        )

    present_methods <-
        paper_method_order()[
            paper_method_order() %in%
                as.character(
                    unique(visible_method_data$method)
                )
        ]

    p <- ggplot2::ggplot(
        method_data,
        ggplot2::aes(
            x = n_clusters,
            y = estimate_plot,
            group = method,
            colour = method,
            shape = method
        )
    )

    if (
        isTRUE(show_nominal_coverage) &&
        identical(metric, "coverage")
    ) {
        p <- p +
            ggplot2::geom_hline(
                yintercept = 0.95,
                linetype = "dotted",
                linewidth = 0.35,
                colour = "grey30"
            )
    }

    if (isTRUE(show_interval)) {
        p <- p +
            ggplot2::geom_ribbon(
                ggplot2::aes(
                    ymin = lower_plot,
                    ymax = upper_plot,
                    fill = method,
                    colour = NULL
                ),
                alpha = 0.12,
                na.rm = TRUE,
                show.legend = FALSE
            )
    }

    p <- p +
        ggplot2::geom_line(
            linewidth = 0.45,
            na.rm = TRUE
        ) +
        ggplot2::geom_point(
            size = 1.55,
            stroke = 0.55,
            na.rm = TRUE
        )

    if (
        isTRUE(show_oracle) &&
        !identical(metric, "time") &&
        nrow(oracle_data) > 0L
    ) {
        p <- p +
            ggplot2::geom_line(
                data = oracle_data,
                mapping = ggplot2::aes(
                    x = n_clusters,
                    y = estimate_plot,
                    group =
                        interaction(
                            reference_type,
                            n_obs_per_cluster
                        )
                ),
                inherit.aes = FALSE,
                linetype = "dashed",
                linewidth = 0.55,
                colour = "black",
                na.rm = TRUE
            )
    }

    p <- p +
        ggplot2::facet_wrap(
            ggplot2::vars(
                n_obs_per_cluster
            ),
            nrow = 1,
            strip.position = "top",
            labeller = if (isTRUE(show_facet_labels)) {
                .paper_facet_labeller(
                    tex = tex
                )
            } else {
                .paper_blank_facet_labeller()
            }
        ) +
        ggplot2::scale_x_continuous(
            breaks = c(
                100,
                200,
                300,
                400,
                500
            )
        ) +
        ggplot2::scale_colour_manual(
            values =
                paper_method_colours(),
            breaks = present_methods,
            labels =
                paper_method_display_labels()[
                    present_methods
                ],
            drop = FALSE
        ) +
        ggplot2::scale_fill_manual(
            values =
                paper_method_colours(),
            breaks = present_methods,
            labels =
                paper_method_display_labels()[
                    present_methods
                ],
            drop = FALSE
        ) +
        ggplot2::scale_shape_manual(
            values =
                paper_method_shapes(),
            breaks = present_methods,
            labels =
                paper_method_display_labels()[
                    present_methods
                ],
            drop = FALSE
        ) +
        ggplot2::labs(
            x = if (isTRUE(show_axis_labels)) {
                if (isTRUE(tex)) {
                    "$d$"
                } else {
                    "Number of clusters, d"
                }
            } else {
                NULL
            },
            y = if (isTRUE(show_axis_labels)) {
                paper_metric_y_label(
                    metric,
                    tex = tex
                )
            } else {
                NULL
            }
        ) +
        theme_simstudy_paper(
            base_size = base_size
        )

    if (!isTRUE(show_facet_labels)) {
        p <- p +
            ggplot2::theme(
                strip.text =
                    ggplot2::element_blank(),
                strip.background =
                    ggplot2::element_blank()
            )
    }

    if (
        .paper_metric_uses_oracle_limit(metric) &&
        is.finite(y_upper)
    ) {
        p <- p +
            ggplot2::coord_cartesian(
                ylim = c(
                    0,
                    y_upper
                )
            )
    } else if (
        metric %in% c(
            "RMISE",
            "rmw2"
        )
    ) {
        ## For error plots, use zero as the lower limit even when the
        ## upper limit is left data-driven.  This keeps the RMISE and
        ## population-error panels comparable and avoids a visually
        ## exaggerated y-axis when all methods are close together.
        p <- p +
            ggplot2::coord_cartesian(
                ylim = c(
                    0,
                    NA_real_
                )
            )
    }

    if (
        identical(metric, "coverage") &&
        !is.null(coverage_y_limits)
    ) {
        p <- p +
            ggplot2::coord_cartesian(
                ylim = coverage_y_limits
            )
    }

    if (identical(metric, "time")) {
        time_values <- c(
            method_data$estimate_plot,
            method_data$lower_plot,
            method_data$upper_plot
        )

        time_values <- time_values[
            is.finite(time_values)
        ]

        if (
            length(time_values) > 0L &&
            max(
                time_values,
                na.rm = TRUE
            ) > 30
        ) {
            p <- p +
                ggplot2::coord_cartesian(
                    ylim = if (isTRUE(log_time)) {
                        c(
                            NA_real_,
                            30
                        )
                    } else {
                        c(
                            0,
                            30
                        )
                    }
                )
        } else if (!isTRUE(log_time)) {
            p <- p +
                ggplot2::coord_cartesian(
                    ylim = c(
                        0,
                        NA_real_
                    )
                )
        }
    }

    if (
        identical(metric, "time") &&
        isTRUE(log_time)
    ) {
        p <- p +
            ggplot2::scale_y_log10()
    }

    if (!isTRUE(legend)) {
        p <- p +
            ggplot2::theme(
                legend.position = "none"
            )
    } else {
        p <- p +
            ggplot2::guides(
                colour =
                    ggplot2::guide_legend(
                        order = 1
                    ),
                shape =
                    ggplot2::guide_legend(
                        order = 1
                    ),
                fill = "none"
            )
    }

    p
}


make_case_simstudy_plots <- function(
    analysis_results,
    legend = FALSE,
    tex = TRUE,
    show_interval = TRUE,
    log_time = FALSE,
    oracle_y_multiplier = 5,
    clip_to_oracle_multiple = TRUE,
    coverage_ylim = c(0.7, 1),
    coverage_lower_padding = 0.02,
    base_size = 9.5,
    show_axis_labels = FALSE,
    show_facet_labels = c(
        RMISE = TRUE,
        coverage = FALSE,
        ci_width = FALSE,
        rmw2 = FALSE,
        time = FALSE
    ),
    apply_plot_rules = TRUE,
    local_gam_min_ni_for_error = 5,
    local_gam_show_ci = FALSE,
    plot_exclusions_file = NULL
) {
    if (
        !is.null(plot_exclusions_file) &&
        file.exists(plot_exclusions_file)
    ) {
        file.remove(plot_exclusions_file)
    }

    metrics <- c(
        "RMISE",
        "coverage",
        "ci_width",
        "rmw2",
        "time"
    )

    plots <- lapply(
        metrics,
        function(metric) {
            plot_simstudy_metric(
                analysis_results =
                    analysis_results,
                metric = metric,
                show_interval =
                    show_interval,
                show_oracle = TRUE,
                show_nominal_coverage =
                    TRUE,
                legend = legend,
                tex = tex,
                log_time =
                    log_time,
                oracle_y_multiplier =
                    oracle_y_multiplier,
                clip_to_oracle_multiple =
                    clip_to_oracle_multiple,
                coverage_ylim =
                    coverage_ylim,
                coverage_lower_padding =
                    coverage_lower_padding,
                base_size =
                    base_size,
                show_axis_labels =
                    show_axis_labels,
                show_facet_labels =
                    isTRUE(
                        show_facet_labels[[metric]]
                    ),
                apply_plot_rules =
                    apply_plot_rules,
                local_gam_min_ni_for_error =
                    local_gam_min_ni_for_error,
                local_gam_show_ci =
                    local_gam_show_ci,
                plot_exclusions_file =
                    plot_exclusions_file
            )
        }
    )

    names(plots) <- metrics
    plots
}


## Save a ggplot as TikZ/LaTeX code, PDF, or both.
save_simstudy_plot <- function(
    plot,
    filename_stem,
    output_dir,
    formats = c(
        "tikz",
        "pdf"
    ),
    width = 6.00,
    height = 1.3,
    stand_alone = FALSE,
    tikz_packages = NULL
) {
    formats <- match.arg(
        formats,
        choices = c(
            "tikz",
            "pdf"
        ),
        several.ok = TRUE
    )

    dir.create(
        output_dir,
        recursive = TRUE,
        showWarnings = FALSE
    )

    written <- character()

    if ("tikz" %in% formats) {
        if (
            !requireNamespace(
                "tikzDevice",
                quietly = TRUE
            )
        ) {
            stop(
                "Package tikzDevice is required for TikZ output.",
                call. = FALSE
            )
        }

        tikz_file <- file.path(
            output_dir,
            paste0(
                filename_stem,
                ".tex"
            )
        )

        tikz_args <- list(
            file = tikz_file,
            width = width,
            height = height,
            standAlone = stand_alone
        )

        if (!is.null(tikz_packages)) {
            tikz_args$packages <-
                tikz_packages
        }

        do.call(
            tikzDevice::tikz,
            tikz_args
        )

        tryCatch(
            print(plot),
            finally = {
                grDevices::dev.off()
            }
        )

        written <- c(
            written,
            tikz_file
        )
    }

    if ("pdf" %in% formats) {
        pdf_file <- file.path(
            output_dir,
            paste0(
                filename_stem,
                ".pdf"
            )
        )

        ggplot2::ggsave(
            filename = pdf_file,
            plot = plot,
            width = width,
            height = height,
            units = "in"
        )

        written <- c(
            written,
            pdf_file
        )
    }

    invisible(written)
}


## Save all five case-level panels.
save_case_simstudy_plots <- function(
    plots,
    case_name,
    output_dir,
    formats = c(
        "tikz",
        "pdf"
    ),
    width = 6.00,
    heights = c(
        RMISE = 1.5,
        coverage = 1.3,
        ci_width = 1.3,
        rmw2 = 1.3,
        time = 1.3
    ),
    stand_alone = FALSE,
    tikz_packages = NULL
) {
    required_plots <- c(
        "RMISE",
        "coverage",
        "ci_width",
        "rmw2",
        "time"
    )

    missing_plots <- setdiff(
        required_plots,
        names(plots)
    )

    if (length(missing_plots) > 0L) {
        stop(
            "plots is missing: ",
            paste(missing_plots, collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    filename_metric <- c(
        RMISE = "RMISE",
        coverage = "coverage",
        ci_width = "ci-width",
        rmw2 = "rmw2",
        time = "time"
    )

    written <- purrr::map(
        required_plots,
        function(metric) {
            save_simstudy_plot(
                plot = plots[[metric]],
                filename_stem =
                    paste0(
                        filename_metric[[metric]],
                        "-",
                        case_name
                    ),
                output_dir = output_dir,
                formats = formats,
                width = width,
                height =
                    unname(
                        heights[[metric]]
                    ),
                stand_alone =
                    stand_alone,
                tikz_packages =
                    tikz_packages
            )
        }
    )

    names(written) <- required_plots
    invisible(written)
}
