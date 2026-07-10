library(tidyverse)
devtools::load_all()


## ==================================================================
## User choices
## ==================================================================

## Number of seeds which will ultimately be used for each retained
## method--subcase cell in the full simulation studies.
n_full_seeds <- 100L

## Candidate common fitting-time cutoffs, in minutes.
candidate_fit_cutoffs <- c(20, 25, 30)

## Oracle benchmark methods are retained irrespective of fitting time.
## They are still subject to the memory, timeout and reliability rules,
## and their elapsed time remains in the projected computational cost.
fit_time_exempt_methods <- c(
    "Oracle-GP",
    "Oracle-SITAR"
)

## First run:
##   leave this as NA_real_ and inspect cutoff_cost_comparison.csv.
##
## Second run:
##   replace NA_real_ by the chosen common cutoff, for example 30.
chosen_fit_cutoff <- 30

## Extra feasibility seeds for cells classified as borderline.
borderline_seeds <- 10006:10010

## Buffer applied when estimating the resources required for the full
## simulation study.
scheduling_buffer <- 1.25


## ==================================================================
## Paths
## ==================================================================

analysis_dir <- here::here(
    "instances",
    "feasibility_analysis"
)

analysis_output_dir <- file.path(
    analysis_dir,
    "output"
)

dir.create(
    analysis_output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

initial_status_files <- c(
    "2re" = here::here(
        "instances",
        "2re",
        "feasibility",
        "output",
        "complete_status_initial.rds"
    ),
    "3re" = here::here(
        "instances",
        "3re",
        "feasibility",
        "output",
        "complete_status_initial.rds"
    ),
    "sitar" = here::here(
        "instances",
        "sitar",
        "feasibility",
        "output",
        "complete_status_initial.rds"
    )
)


## ==================================================================
## Read the initial feasibility results
## ==================================================================

missing_status_files <- initial_status_files[
    !file.exists(initial_status_files)
]

if (length(missing_status_files) > 0L) {
    stop(
        paste0(
            "The following initial feasibility files are missing:\n",
            paste(
                paste0(
                    names(missing_status_files),
                    ": ",
                    missing_status_files
                ),
                collapse = "\n"
            )
        ),
        call. = FALSE
    )
}

complete_status_initial <- purrr::map_dfr(
    initial_status_files,
    readRDS
)

required_status_columns <- c(
    "case",
    "subcase_id",
    "method",
    "seed"
)

missing_status_columns <- setdiff(
    required_status_columns,
    names(complete_status_initial)
)

if (length(missing_status_columns) > 0L) {
    stop(
        "The combined complete-status table is missing: ",
        paste(missing_status_columns, collapse = ", "),
        call. = FALSE
    )
}

duplicated_runs <- complete_status_initial %>%
    dplyr::count(
        case,
        subcase_id,
        method,
        seed,
        name = "n_records"
    ) %>%
    dplyr::filter(n_records > 1L)

if (nrow(duplicated_runs) > 0L) {
    print(duplicated_runs)

    stop(
        "The initial feasibility results contain duplicated runs.",
        call. = FALSE
    )
}

saveRDS(
    complete_status_initial,
    file.path(
        analysis_output_dir,
        "complete_status_initial_all_cases.rds"
    )
)


## ==================================================================
## Construct one summary row per method--subcase cell
## ==================================================================

cell_summary_initial <- summarise_feasibility_cells(
    complete_status_initial
) %>%
    dplyr::arrange(
        case,
        subcase_id,
        method
    )

readr::write_csv(
    cell_summary_initial,
    file.path(
        analysis_output_dir,
        "feasibility_cell_summary_initial.csv"
    )
)


## ==================================================================
## Compare candidate common fitting-time cutoffs
## ==================================================================

cutoff_comparison <- compare_feasibility_cutoffs(
    feasibility_summary = cell_summary_initial,
    cutoffs = candidate_fit_cutoffs,
    n_full_seeds = n_full_seeds,
    fit_time_exempt_methods = fit_time_exempt_methods,
    stage = "initial",
    scheduling_buffer = scheduling_buffer
)

readr::write_csv(
    cutoff_comparison,
    file.path(
        analysis_output_dir,
        "cutoff_cost_comparison.csv"
    )
)

print(cutoff_comparison)


## ==================================================================
## Apply the chosen cutoff and create borderline setups
## ==================================================================

if (is.na(chosen_fit_cutoff)) {
    message(
        "\nCreated cutoff_cost_comparison.csv.\n",
        "Inspect that table, set chosen_fit_cutoff near the top of ",
        "this script, and run the script again to create the initial ",
        "decisions and borderline setups."
    )
} else {
    if (
        length(chosen_fit_cutoff) != 1L ||
        !is.finite(chosen_fit_cutoff) ||
        chosen_fit_cutoff <= 0
    ) {
        stop(
            "chosen_fit_cutoff must be one positive finite number.",
            call. = FALSE
        )
    }

    if (!chosen_fit_cutoff %in% candidate_fit_cutoffs) {
        stop(
            "chosen_fit_cutoff must currently be one of: ",
            paste(candidate_fit_cutoffs, collapse = ", "),
            ". Add it to candidate_fit_cutoffs first if another ",
            "cutoff is required.",
            call. = FALSE
        )
    }

    decisions_initial <- classify_feasibility_cells(
        feasibility_summary = cell_summary_initial,
        fit_cutoff_minutes = chosen_fit_cutoff,
        fit_time_exempt_methods = fit_time_exempt_methods,
        stage = "initial"
    ) %>%
        dplyr::arrange(
            case,
            subcase_id,
            method
        )

    readr::write_csv(
        decisions_initial,
        file.path(
            analysis_output_dir,
            "feasibility_decisions_initial.csv"
        )
    )

    chosen_cutoff_record <- tibble::tibble(
        fit_cutoff_minutes = chosen_fit_cutoff,
        n_full_seeds = n_full_seeds,
        scheduling_buffer = scheduling_buffer,
        fit_time_exempt_methods = paste(
            fit_time_exempt_methods,
            collapse = ";"
        ),
        n_cells = nrow(decisions_initial),
        n_include = sum(
            decisions_initial$decision == "include"
        ),
        n_borderline = sum(
            decisions_initial$decision == "borderline"
        ),
        n_excluded = sum(
            grepl(
                "^exclude_",
                decisions_initial$decision
            )
        ),
        created_at = Sys.time()
    )

    readr::write_csv(
        chosen_cutoff_record,
        file.path(
            analysis_output_dir,
            "chosen_fit_cutoff.csv"
        )
    )

    borderline_setup_all <- make_borderline_setup(
        decisions = decisions_initial,
        extra_seeds = borderline_seeds
    ) %>%
        dplyr::arrange(
            case,
            subcase_id,
            method,
            seed
        )

    readr::write_csv(
        borderline_setup_all,
        file.path(
            analysis_output_dir,
            "borderline_setup_all_cases.csv"
        )
    )


    ## --------------------------------------------------------------
    ## Write the case-specific files read by run_borderline.R
    ## --------------------------------------------------------------

    case_borderline_files <- c(
        "2re" = here::here(
            "instances",
            "2re",
            "feasibility",
            "output",
            "borderline_setup.csv"
        ),
        "3re" = here::here(
            "instances",
            "3re",
            "feasibility",
            "output",
            "borderline_setup.csv"
        ),
        "sitar" = here::here(
            "instances",
            "sitar",
            "feasibility",
            "output",
            "borderline_setup.csv"
        )
    )

    for (case_name in names(case_borderline_files)) {
        output_file <- case_borderline_files[[case_name]]

        setup_case <- borderline_setup_all %>%
            dplyr::filter(case == case_name) %>%
            dplyr::select(
                case,
                subcase_id,
                method,
                seed
            ) %>%
            dplyr::arrange(
                subcase_id,
                method,
                seed
            )

        dir.create(
            dirname(output_file),
            recursive = TRUE,
            showWarnings = FALSE
        )

        readr::write_csv(
            setup_case,
            output_file
        )

        message(
            "Wrote ",
            nrow(setup_case),
            " additional feasibility runs for ",
            case_name,
            " to ",
            output_file
        )
    }


    ## --------------------------------------------------------------
    ## Brief decision summary
    ## --------------------------------------------------------------

    decision_counts <- decisions_initial %>%
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
            "feasibility_decision_counts_initial.csv"
        )
    )

    print(decision_counts)

    if (nrow(borderline_setup_all) == 0L) {
        message(
            "\nNo cells were classified as borderline. ",
            "No additional feasibility runs are required."
        )
    } else {
        message(
            "\nCreated ",
            nrow(borderline_setup_all),
            " additional runs across ",
            dplyr::n_distinct(
                paste(
                    borderline_setup_all$case,
                    borderline_setup_all$subcase_id,
                    borderline_setup_all$method,
                    sep = "::"
                )
            ),
            " borderline method--subcase cells."
        )
    }
}
