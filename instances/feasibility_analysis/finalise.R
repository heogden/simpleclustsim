library(tidyverse)
devtools::load_all()

analysis_output_dir <- here::here(
    "instances",
    "feasibility_analysis",
    "output"
)

decisions_initial_file <- file.path(
    analysis_output_dir,
    "feasibility_decisions_initial.csv"
)

cutoff_file <- file.path(
    analysis_output_dir,
    "chosen_fit_cutoff.csv"
)

required_files <- c(
    decisions_initial_file,
    cutoff_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
    stop(
        "Missing required files:\n",
        paste(missing_files, collapse = "\n"),
        call. = FALSE
    )
}

decisions_initial <- readr::read_csv(
    decisions_initial_file,
    show_col_types = FALSE
)

feasibility_settings <- readr::read_csv(
    cutoff_file,
    show_col_types = FALSE
)

fit_cutoffs <- unique(
    feasibility_settings$fit_cutoff_minutes
)

n_full_seeds_values <- unique(
    feasibility_settings$n_full_seeds
)

scheduling_buffers <- unique(
    feasibility_settings$scheduling_buffer
)

if (!"fit_time_exempt_methods" %in% names(feasibility_settings)) {
    stop(
        "chosen_fit_cutoff.csv is missing fit_time_exempt_methods. ",
        "Rerun assess_initial.R with the oracle exemptions recorded.",
        call. = FALSE
    )
}

fit_time_exempt_values <- unique(
    feasibility_settings$fit_time_exempt_methods
)

if (length(fit_cutoffs) != 1L) {
    stop(
        "chosen_fit_cutoff.csv must contain one fitting-time cutoff.",
        call. = FALSE
    )
}

if (length(n_full_seeds_values) != 1L) {
    stop(
        "chosen_fit_cutoff.csv must contain one value of n_full_seeds.",
        call. = FALSE
    )
}

if (length(scheduling_buffers) != 1L) {
    stop(
        "chosen_fit_cutoff.csv must contain one scheduling buffer.",
        call. = FALSE
    )
}

if (length(fit_time_exempt_values) != 1L) {
    stop(
        "chosen_fit_cutoff.csv must contain one value of ",
        "fit_time_exempt_methods.",
        call. = FALSE
    )
}

chosen_fit_cutoff <- fit_cutoffs[[1]]
n_full_seeds <- as.integer(n_full_seeds_values[[1]])
scheduling_buffer <- scheduling_buffers[[1]]

fit_time_exempt_value <- fit_time_exempt_values[[1]]

if (
    is.na(fit_time_exempt_value) ||
    !nzchar(trimws(fit_time_exempt_value))
) {
    fit_time_exempt_methods <- character()
} else {
    fit_time_exempt_methods <- trimws(
        strsplit(
            fit_time_exempt_value,
            split = ";",
            fixed = TRUE
        )[[1]]
    )
}


## ==================================================================
## Read the initial feasibility evidence
## ==================================================================

initial_status_files <- c(
    "2re" = here::here(
        "instances", "2re", "feasibility", "output",
        "complete_status_initial.rds"
    ),
    "3re" = here::here(
        "instances", "3re", "feasibility", "output",
        "complete_status_initial.rds"
    ),
    "sitar" = here::here(
        "instances", "sitar", "feasibility", "output",
        "complete_status_initial.rds"
    )
)

missing_initial_files <- initial_status_files[
    !file.exists(initial_status_files)
]

if (length(missing_initial_files) > 0L) {
    stop(
        "Missing initial feasibility files:\n",
        paste(missing_initial_files, collapse = "\n"),
        call. = FALSE
    )
}

complete_status_initial <- purrr::map_dfr(
    initial_status_files,
    readRDS
)


## ==================================================================
## Read additional evidence for initially borderline cells
## ==================================================================

borderline_cells <- decisions_initial %>%
    dplyr::filter(
        decision == "borderline"
    ) %>%
    dplyr::distinct(
        case,
        subcase_id,
        method
    )

borderline_status_files <- c(
    "2re" = here::here(
        "instances", "2re", "feasibility", "output",
        "complete_status_borderline.rds"
    ),
    "3re" = here::here(
        "instances", "3re", "feasibility", "output",
        "complete_status_borderline.rds"
    ),
    "sitar" = here::here(
        "instances", "sitar", "feasibility", "output",
        "complete_status_borderline.rds"
    )
)

borderline_cases <- unique(
    borderline_cells$case
)

if (length(borderline_cases) > 0L) {
    required_borderline_files <- borderline_status_files[
        borderline_cases
    ]

    missing_borderline_files <- required_borderline_files[
        !file.exists(required_borderline_files)
    ]

    if (length(missing_borderline_files) > 0L) {
        stop(
            "Borderline feasibility runs are required, but these files ",
            "are missing:\n",
            paste(missing_borderline_files, collapse = "\n"),
            call. = FALSE
        )
    }

    complete_status_borderline <- purrr::map_dfr(
        required_borderline_files,
        readRDS
    )
} else {
    complete_status_borderline <- tibble::tibble()
}

complete_status_all <- dplyr::bind_rows(
    complete_status_initial,
    complete_status_borderline
)

duplicated_runs <- complete_status_all %>%
    dplyr::count(
        case,
        subcase_id,
        method,
        seed,
        name = "n_records"
    ) %>%
    dplyr::filter(
        n_records > 1L
    )

if (nrow(duplicated_runs) > 0L) {
    print(duplicated_runs)

    stop(
        "The combined feasibility evidence contains duplicated runs.",
        call. = FALSE
    )
}

saveRDS(
    complete_status_all,
    file.path(
        analysis_output_dir,
        "complete_status_all_cases.rds"
    )
)


## ==================================================================
## Re-summarise using five or ten runs as appropriate
## ==================================================================

cell_summary_final <- summarise_feasibility_cells(
    complete_status_all
) %>%
    dplyr::arrange(
        case,
        subcase_id,
        method
    )

readr::write_csv(
    cell_summary_final,
    file.path(
        analysis_output_dir,
        "feasibility_cell_summary_final.csv"
    )
)


## ==================================================================
## Make final decisions for the initially borderline cells
## ==================================================================

if (nrow(borderline_cells) > 0L) {
    borderline_summary_final <- cell_summary_final %>%
        dplyr::semi_join(
            borderline_cells,
            by = c(
                "case",
                "subcase_id",
                "method"
            )
        )

    unexpected_trial_counts <- borderline_summary_final %>%
        dplyr::filter(
            n_trials != 10L
        )

    if (nrow(unexpected_trial_counts) > 0L) {
        print(unexpected_trial_counts)

        stop(
            "Each initially borderline cell should have ten feasibility ",
            "runs before final decisions are made.",
            call. = FALSE
        )
    }

    borderline_decisions_final <- classify_feasibility_cells(
        feasibility_summary = borderline_summary_final,
        fit_cutoff_minutes = chosen_fit_cutoff,
        fit_time_exempt_methods = fit_time_exempt_methods,
        stage = "final"
    )

    decisions_final <- finalise_feasibility_decisions(
        initial_decisions = decisions_initial,
        borderline_final_decisions = borderline_decisions_final
    )
} else {
    borderline_decisions_final <- tibble::tibble()
    decisions_final <- decisions_initial
}

decisions_final <- decisions_final %>%
    dplyr::arrange(
        case,
        subcase_id,
        method
    )

readr::write_csv(
    borderline_decisions_final,
    file.path(
        analysis_output_dir,
        "feasibility_decisions_borderline_final.csv"
    )
)

readr::write_csv(
    decisions_final,
    file.path(
        analysis_output_dir,
        "feasibility_decisions.csv"
    )
)


## ==================================================================
## Final projected cost
## ==================================================================

final_cost <- estimate_full_simulation_cost(
    feasibility_summary = cell_summary_final,
    decisions = decisions_final,
    n_full_seeds = n_full_seeds,
    included_decisions = "include",
    scheduling_buffer = scheduling_buffer
)

readr::write_csv(
    final_cost$by_cell,
    file.path(
        analysis_output_dir,
        "projected_full_cost_by_cell.csv"
    )
)

readr::write_csv(
    final_cost$by_case,
    file.path(
        analysis_output_dir,
        "projected_full_cost_by_case.csv"
    )
)

readr::write_csv(
    final_cost$overall,
    file.path(
        analysis_output_dir,
        "projected_full_cost_overall.csv"
    )
)


## ==================================================================
## Report unresolved cells
## ==================================================================

decision_counts <- decisions_final %>%
    dplyr::count(
        case,
        decision,
        name = "n_cells"
    ) %>%
    dplyr::arrange(
        case,
        decision
    )

readr::write_csv(
    decision_counts,
    file.path(
        analysis_output_dir,
        "feasibility_decision_counts_final.csv"
    )
)

print(decision_counts)
print(final_cost$by_case)
print(final_cost$overall)

manual_review_cells <- decisions_final %>%
    dplyr::filter(
        decision == "manual_review"
    )

if (nrow(manual_review_cells) > 0L) {
    readr::write_csv(
        manual_review_cells,
        file.path(
            analysis_output_dir,
            "feasibility_manual_review.csv"
        )
    )

    warning(
        nrow(manual_review_cells),
        " cells require manual review before starting the full runs.",
        call. = FALSE
    )
} else {
    message(
        "All feasibility decisions are resolved. ",
        "The full simulation runs can now be prepared."
    )
}
