## Extract numerical summaries used to describe simulation results in the paper.
##
## Intended location:
##   instances/<case>/analysis/describing_results.R
##
## Change case_name below by hand, then source this file.
##
## The script reads:
##   instances/<case>/analysis/output/analysis_results.rds
##
## and writes:
##   instances/<case>/analysis/output/describing_results_summary.csv
##   instances/<case>/analysis/output/describing_results_details.csv
##   instances/<case>/analysis/output/describing_results_text.txt

project_root <- here::here()
devtools::load_all(project_root)

## Change this by hand for the case you want to summarise.
case_name <- "sitar"


analysis_dir <- here::here(
    "instances",
    case_name,
    "analysis",
    "output"
)

analysis_results_file <- file.path(
    analysis_dir,
    "analysis_results.rds"
)

if (!file.exists(analysis_results_file)) {
    stop(
        "Cannot find ",
        analysis_results_file,
        ". Run the case analysis script first.",
        call. = FALSE
    )
}

analysis_results <- readRDS(
    analysis_results_file
)

method_plot_data <- analysis_results$method_plot_data
oracle_reference_data <- analysis_results$oracle_reference_data

if (is.null(method_plot_data)) {
    stop(
        "analysis_results$method_plot_data is missing.",
        call. = FALSE
    )
}

if (
    is.null(oracle_reference_data) ||
    nrow(oracle_reference_data) == 0L
) {
    ## Fallback for older outputs: try to find oracle rows in method_plot_data.
    oracle_reference_data <- method_plot_data |>
        dplyr::filter(
            grepl(
                "^Oracle",
                method
            )
        )
}

required_method_cols <- c(
    "metric",
    "method",
    "n_obs_per_cluster",
    "n_clusters",
    "estimate"
)

missing_method_cols <- setdiff(
    required_method_cols,
    names(method_plot_data)
)

if (length(missing_method_cols) > 0L) {
    stop(
        "method_plot_data is missing required columns: ",
        paste(missing_method_cols, collapse = ", "),
        ".",
        call. = FALSE
    )
}

required_oracle_cols <- c(
    "metric",
    "n_obs_per_cluster",
    "n_clusters",
    "estimate"
)

missing_oracle_cols <- setdiff(
    required_oracle_cols,
    names(oracle_reference_data)
)

if (length(missing_oracle_cols) > 0L) {
    stop(
        "oracle_reference_data is missing required columns: ",
        paste(missing_oracle_cols, collapse = ", "),
        ".",
        call. = FALSE
    )
}


if (!("show_metric" %in% names(method_plot_data))) {
    method_plot_data$show_metric <- TRUE
}

if (!("show_metric" %in% names(oracle_reference_data))) {
    oracle_reference_data$show_metric <- TRUE
}


## These are the quantities used in the narrative text.
quantity_specs <- tibble::tibble(
    quantity = c(
        "Next-best general-purpose RMISE / AdaStruMM RMISE",
        "AdaStruMM RMISE / oracle RMISE",
        "AdaStruMM CI coverage",
        "Narrowest general-purpose CI width / AdaStruMM CI width",
        "AdaStruMM CI width / oracle CI width",
        "Next-best general-purpose RM2W / AdaStruMM RM2W",
        "AdaStruMM RM2W / oracle RM2W"
    ),
    metric = c(
        "RMISE",
        "RMISE",
        "coverage",
        "ci_width",
        "ci_width",
        "rmw2",
        "rmw2"
    ),
    order = seq_len(7L)
)


general_purpose_methods <- function(method) {
    !is.na(method) &
        method != "AdaStruMM" &
        method != "SITAR" &
        !grepl(
            "^Oracle",
            method
        )
}


clean_method_data <- method_plot_data |>
    dplyr::mutate(
        show_metric =
            dplyr::coalesce(
                show_metric,
                TRUE
            )
    ) |>
    dplyr::filter(
        show_metric,
        is.finite(estimate)
    )


clean_oracle_data <- oracle_reference_data |>
    dplyr::mutate(
        show_metric =
            dplyr::coalesce(
                show_metric,
                TRUE
            )
    ) |>
    dplyr::filter(
        show_metric,
        is.finite(estimate)
    ) |>
    dplyr::group_by(
        metric,
        n_obs_per_cluster,
        n_clusters
    ) |>
    dplyr::summarise(
        oracle_estimate = min(
            estimate,
            na.rm = TRUE
        ),
        .groups = "drop"
    )


best_general_method <- function(metric_name) {
    clean_method_data |>
        dplyr::filter(
            .data$metric == .env$metric_name,
            general_purpose_methods(method)
        ) |>
        dplyr::group_by(
            n_obs_per_cluster,
            n_clusters
        ) |>
        dplyr::slice_min(
            order_by = estimate,
            n = 1,
            with_ties = FALSE
        ) |>
        dplyr::ungroup() |>
        dplyr::transmute(
            metric = .env$metric_name,
            n_obs_per_cluster,
            n_clusters,
            comparison_method = method,
            comparison_estimate = estimate
        )
}


adastrumm_metric <- function(metric_name) {
    clean_method_data |>
        dplyr::filter(
            .data$metric == .env$metric_name,
            method == "AdaStruMM"
        ) |>
        dplyr::transmute(
            metric,
            n_obs_per_cluster,
            n_clusters,
            adastrumm_estimate = estimate
        )
}


oracle_metric <- function(metric_name) {
    clean_oracle_data |>
        dplyr::filter(
            .data$metric == .env$metric_name
        ) |>
        dplyr::transmute(
            metric,
            n_obs_per_cluster,
            n_clusters,
            oracle_estimate
        )
}


make_best_general_ratio <- function(
    metric_name,
    quantity,
    numerator_label,
    denominator_label
) {
    best_general_method(metric_name) |>
        dplyr::inner_join(
            adastrumm_metric(metric_name),
            by = c(
                "metric",
                "n_obs_per_cluster",
                "n_clusters"
            )
        ) |>
        dplyr::filter(
            is.finite(comparison_estimate),
            is.finite(adastrumm_estimate),
            adastrumm_estimate > 0
        ) |>
        dplyr::mutate(
            case = case_name,
            quantity = quantity,
            numerator = numerator_label,
            denominator = denominator_label,
            numerator_estimate = comparison_estimate,
            denominator_estimate = adastrumm_estimate,
            value = comparison_estimate / adastrumm_estimate
        ) |>
        dplyr::relocate(
            case,
            quantity,
            metric,
            n_obs_per_cluster,
            n_clusters
        )
}


make_adastrumm_oracle_ratio <- function(
    metric_name,
    quantity,
    numerator_label,
    denominator_label
) {
    adastrumm_metric(metric_name) |>
        dplyr::inner_join(
            oracle_metric(metric_name),
            by = c(
                "metric",
                "n_obs_per_cluster",
                "n_clusters"
            )
        ) |>
        dplyr::filter(
            is.finite(adastrumm_estimate),
            is.finite(oracle_estimate),
            oracle_estimate > 0
        ) |>
        dplyr::mutate(
            case = case_name,
            quantity = quantity,
            numerator = numerator_label,
            denominator = denominator_label,
            numerator_estimate = adastrumm_estimate,
            denominator_estimate = oracle_estimate,
            comparison_method = "Oracle",
            value = adastrumm_estimate / oracle_estimate
        ) |>
        dplyr::relocate(
            case,
            quantity,
            metric,
            n_obs_per_cluster,
            n_clusters
        )
}


make_adastrumm_coverage <- function() {
    clean_method_data |>
        dplyr::filter(
            metric == "coverage",
            method == "AdaStruMM"
        ) |>
        dplyr::transmute(
            case = case_name,
            quantity = "AdaStruMM CI coverage",
            metric,
            n_obs_per_cluster,
            n_clusters,
            comparison_method = NA_character_,
            comparison_estimate = NA_real_,
            adastrumm_estimate = estimate,
            oracle_estimate = NA_real_,
            numerator = "AdaStruMM coverage",
            denominator = NA_character_,
            numerator_estimate = estimate,
            denominator_estimate = NA_real_,
            value = estimate
        )
}


details <- dplyr::bind_rows(
    make_best_general_ratio(
        metric_name = "RMISE",
        quantity = "Next-best general-purpose RMISE / AdaStruMM RMISE",
        numerator_label =
            "minimum RMISE among non-AdaStruMM general-purpose methods",
        denominator_label =
            "AdaStruMM RMISE"
    ),
    make_adastrumm_oracle_ratio(
        metric_name = "RMISE",
        quantity = "AdaStruMM RMISE / oracle RMISE",
        numerator_label =
            "AdaStruMM RMISE",
        denominator_label =
            "oracle RMISE"
    ),
    make_adastrumm_coverage(),
    make_best_general_ratio(
        metric_name = "ci_width",
        quantity =
            "Narrowest general-purpose CI width / AdaStruMM CI width",
        numerator_label =
            "minimum CI width among non-AdaStruMM general-purpose methods",
        denominator_label =
            "AdaStruMM CI width"
    ),
    make_adastrumm_oracle_ratio(
        metric_name = "ci_width",
        quantity = "AdaStruMM CI width / oracle CI width",
        numerator_label =
            "AdaStruMM CI width",
        denominator_label =
            "oracle CI width"
    ),
    make_best_general_ratio(
        metric_name = "rmw2",
        quantity = "Next-best general-purpose RM2W / AdaStruMM RM2W",
        numerator_label =
            "minimum RM2W among non-AdaStruMM general-purpose methods",
        denominator_label =
            "AdaStruMM RM2W"
    ),
    make_adastrumm_oracle_ratio(
        metric_name = "rmw2",
        quantity = "AdaStruMM RM2W / oracle RM2W",
        numerator_label =
            "AdaStruMM RM2W",
        denominator_label =
            "oracle RM2W"
    )
) |>
    dplyr::arrange(
        match(
            quantity,
            quantity_specs$quantity
        ),
        n_obs_per_cluster,
        n_clusters
    )


summary_raw <- details |>
    dplyr::group_by(
        case,
        quantity,
        metric
    ) |>
    dplyr::summarise(
        n_cells = sum(
            is.finite(value)
        ),
        min = if (sum(is.finite(value)) > 0L) {
            min(
                value[is.finite(value)],
                na.rm = TRUE
            )
        } else {
            NA_real_
        },
        max = if (sum(is.finite(value)) > 0L) {
            max(
                value[is.finite(value)],
                na.rm = TRUE
            )
        } else {
            NA_real_
        },
        .groups = "drop"
    )


summary_table <- quantity_specs |>
    dplyr::mutate(
        case = case_name,
        .before = quantity
    ) |>
    dplyr::left_join(
        summary_raw,
        by = c(
            "case",
            "quantity",
            "metric"
        )
    ) |>
    dplyr::mutate(
        n_cells = dplyr::coalesce(
            n_cells,
            0L
        ),
        min_1dp = dplyr::if_else(
            is.finite(min),
            formatC(
                min,
                format = "f",
                digits = 1
            ),
            NA_character_
        ),
        max_1dp = dplyr::if_else(
            is.finite(max),
            formatC(
                max,
                format = "f",
                digits = 1
            ),
            NA_character_
        ),
        min_2dp = dplyr::if_else(
            is.finite(min),
            formatC(
                min,
                format = "f",
                digits = 2
            ),
            NA_character_
        ),
        max_2dp = dplyr::if_else(
            is.finite(max),
            formatC(
                max,
                format = "f",
                digits = 2
            ),
            NA_character_
        )
    ) |>
    dplyr::arrange(
        order
    ) |>
    dplyr::select(
        -order
    )


readr::write_csv(
    summary_table,
    file.path(
        analysis_dir,
        "describing_results_summary.csv"
    )
)

readr::write_csv(
    details,
    file.path(
        analysis_dir,
        "describing_results_details.csv"
    )
)


get_summary_row <- function(quantity_name) {
    row <- summary_table |>
        dplyr::filter(
            .data$quantity == .env$quantity_name
        )

    if (nrow(row) == 0L) {
        return(NULL)
    }

    row[1L, , drop = FALSE]
}


format_range <- function(row, digits = 1L) {
    if (
        is.null(row) ||
        nrow(row) == 0L ||
        !is.finite(row$min) ||
        !is.finite(row$max)
    ) {
        return("[missing]")
    }

    if (identical(digits, 1L)) {
        paste0(
            row$min_1dp,
            "--",
            row$max_1dp
        )
    } else {
        paste0(
            row$min_2dp,
            "--",
            row$max_2dp
        )
    }
}


text_lines <- c(
    paste0(
        "Case: ",
        case_name
    ),
    "",
    paste0(
        "RMISE improvement over next-best general-purpose method: ",
        format_range(
            get_summary_row(
                "Next-best general-purpose RMISE / AdaStruMM RMISE"
            ),
            digits = 1L
        ),
        " times."
    ),
    paste0(
        "AdaStruMM RMISE / individual oracle RMISE: ",
        format_range(
            get_summary_row(
                "AdaStruMM RMISE / oracle RMISE"
            ),
            digits = 1L
        ),
        "."
    ),
    paste0(
        "AdaStruMM CI coverage range: ",
        format_range(
            get_summary_row(
                "AdaStruMM CI coverage"
            ),
            digits = 2L
        ),
        "."
    ),
    paste0(
        "CI width improvement over next-narrowest general-purpose method: ",
        format_range(
            get_summary_row(
                "Narrowest general-purpose CI width / AdaStruMM CI width"
            ),
            digits = 1L
        ),
        " times."
    ),
    paste0(
        "AdaStruMM CI width / oracle CI width: ",
        format_range(
            get_summary_row(
                "AdaStruMM CI width / oracle CI width"
            ),
            digits = 1L
        ),
        "."
    ),
    paste0(
        "RM2W improvement over next-best general-purpose method: ",
        format_range(
            get_summary_row(
                "Next-best general-purpose RM2W / AdaStruMM RM2W"
            ),
            digits = 1L
        ),
        " times."
    ),
    paste0(
        "AdaStruMM RM2W / population oracle RM2W: ",
        format_range(
            get_summary_row(
                "AdaStruMM RM2W / oracle RM2W"
            ),
            digits = 1L
        ),
        "."
    )
)

writeLines(
    text_lines,
    con = file.path(
        analysis_dir,
        "describing_results_text.txt"
    )
)

print(summary_table)

message(
    "\nWrote:\n  ",
    file.path(
        analysis_dir,
        "describing_results_summary.csv"
    ),
    "\n  ",
    file.path(
        analysis_dir,
        "describing_results_details.csv"
    ),
    "\n  ",
    file.path(
        analysis_dir,
        "describing_results_text.txt"
    )
)
