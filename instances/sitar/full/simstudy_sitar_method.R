## Additional full simulation run for the fitted SITAR method.
##
## This is needed because fitted SITAR was accidentally excluded from the
## main SITAR full simulation after the feasibility filter treated missing
## CIs as making the method unreliable.
##
## This script writes a separate checkpointed Iridis directory:
##
##   instances/sitar/full/iridis_sitar_method
##
## Run one stage at a time:
##
##   "prepare"
##   "after_4h"
##   "after_8h"
##   "after_recovery"
##
## After all runs are either checkpointed or definitively classified, run
## instances/sitar/full/incorporate_sitar_method.R.

project_root <- here::here()
devtools::load_all(project_root)

case_to_run <- "sitar"
stage <- "after_4h"

initial_walltime_minutes <- 240
extended_walltime_minutes <- 480
mem_gb <- 8
binpacking_back <- 0
scheduling_buffer_if_missing <- 1.25

iridis_dir <- here::here(
    "instances",
    "sitar",
    "full",
    "iridis_sitar_method"
)

recovery_iridis_dir <- here::here(
    "instances",
    "sitar",
    "full",
    "iridis_sitar_method_recovery"
)

diagnostic_dir <- here::here(
    "instances",
    "sitar",
    "full",
    "output",
    "sitar_method_extra"
)

dir.create(
    diagnostic_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

source(
    here::here(
        "instances",
        case_to_run,
        "design.R"
    )
)

if (!identical(case_name, case_to_run)) {
    stop(
        "The sourced design has case_name = ",
        case_name,
        ", but case_to_run = ",
        case_to_run,
        ".",
        call. = FALSE
    )
}


load_rda_object <- function(file, object_name) {
    env <- new.env(parent = emptyenv())
    loaded_names <- load(file, envir = env)

    if (!(object_name %in% loaded_names)) {
        stop(
            file,
            " does not contain an object called ",
            object_name,
            ".",
            call. = FALSE
        )
    }

    env[[object_name]]
}


get_unique_setting <- function(x, column, default = NULL) {
    if (!(column %in% names(x))) {
        if (!is.null(default)) {
            return(default)
        }

        stop(
            "Settings file is missing column ",
            column,
            ".",
            call. = FALSE
        )
    }

    values <- unique(x[[column]])

    if (
        length(values) != 1L ||
        is.na(values[[1]])
    ) {
        if (!is.null(default)) {
            return(default)
        }

        stop(
            "Expected exactly one non-missing value of ",
            column,
            ".",
            call. = FALSE
        )
    }

    values[[1]]
}


read_setup <- function(dir) {
    load_rda_object(
        file.path(
            dir,
            "setup.Rda"
        ),
        "setup"
    )
}


write_diagnostics <- function(
    prefix,
    setup,
    dir,
    time_limit_minutes
) {
    missing_runs <- find_missing_runs_iridis_checkpointed(
        dir,
        setup = setup
    )

    missing_processes <-
        find_processes_with_missing_runs_iridis_checkpointed(
            dir,
            setup = setup
        )

    oom_processes <- find_oom_ids_iridis(
        dir,
        min_mem_gb = mem_gb,
        mem_limit_if_missing_gb = mem_gb
    )

    timeout_processes <- find_timeout_ids_iridis(
        dir,
        min_time_minutes = time_limit_minutes,
        time_limit_if_missing_minutes = time_limit_minutes
    )

    active_runs <- find_active_runs_at_failure_iridis_checkpointed(
        dir,
        process_ids = unique(
            missing_processes$process_id
        )
    )

    readr::write_csv(
        missing_runs,
        file.path(
            diagnostic_dir,
            paste0(prefix, "_missing_runs.csv")
        )
    )

    readr::write_csv(
        missing_processes,
        file.path(
            diagnostic_dir,
            paste0(prefix, "_processes_with_missing_runs.csv")
        )
    )

    readr::write_csv(
        oom_processes,
        file.path(
            diagnostic_dir,
            paste0(prefix, "_oom_processes.csv")
        )
    )

    readr::write_csv(
        timeout_processes,
        file.path(
            diagnostic_dir,
            paste0(prefix, "_timeout_processes.csv")
        )
    )

    readr::write_csv(
        active_runs,
        file.path(
            diagnostic_dir,
            paste0(prefix, "_active_runs_at_failure.csv")
        )
    )

    list(
        missing_runs = missing_runs,
        missing_processes = missing_processes,
        oom_processes = oom_processes,
        timeout_processes = timeout_processes,
        active_runs = active_runs
    )
}


make_sitar_timing_status <- function(
    complete_status_feasibility,
    setup_unallocated,
    diagnostic_dir
) {
    required_columns <- c(
        "case",
        "subcase_id",
        "method",
        "sim_status",
        "fit_status",
        "pred_status",
        "gp_status",
        "time"
    )

    missing_columns <- setdiff(
        required_columns,
        names(complete_status_feasibility)
    )

    if (length(missing_columns) > 0L) {
        stop(
            "complete_status_feasibility is missing: ",
            paste(missing_columns, collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    ## The fitted SITAR feasibility runs can have missing GP summaries or
    ## confidence intervals even when the model fitting/prediction time is
    ## usable.  The Iridis allocation helper treats gp_status == "ok" as
    ## part of "successful timing", so for scheduling only we set gp_status
    ## to "ok" whenever simulation, fitting and prediction succeeded and
    ## a finite elapsed time is available.
    timing_status <- complete_status_feasibility |>
        dplyr::mutate(
            timing_ok =
                !dplyr::coalesce(missing_from_simstudy, TRUE) &
                sim_status == "ok" &
                fit_status == "ok" &
                pred_status == "ok" &
                is.finite(time),

            gp_status = dplyr::if_else(
                timing_ok,
                "ok",
                gp_status
            )
        )

    timing_check <- timing_status |>
        dplyr::group_by(
            case,
            subcase_id,
            method
        ) |>
        dplyr::summarise(
            n_timing_runs = sum(
                timing_ok,
                na.rm = TRUE
            ),
            min_time = suppressWarnings(
                min(
                    time[timing_ok],
                    na.rm = TRUE
                )
            ),
            median_time = suppressWarnings(
                stats::median(
                    time[timing_ok],
                    na.rm = TRUE
                )
            ),
            max_time = suppressWarnings(
                max(
                    time[timing_ok],
                    na.rm = TRUE
                )
            ),
            .groups = "drop"
        ) |>
        dplyr::mutate(
            dplyr::across(
                c(min_time, median_time, max_time),
                ~ dplyr::if_else(
                    is.infinite(.x),
                    NA_real_,
                    .x
                )
            )
        )

    readr::write_csv(
        timing_check,
        file.path(
            diagnostic_dir,
            "sitar_method_feasibility_timing_check.csv"
        )
    )

    print(timing_check)

    missing_timing_cells <- setup_unallocated |>
        dplyr::distinct(
            case,
            subcase_id,
            method
        ) |>
        dplyr::anti_join(
            timing_check |>
                dplyr::filter(
                    n_timing_runs > 0L,
                    is.finite(max_time)
                ) |>
                dplyr::select(
                    case,
                    subcase_id,
                    method
                ),
            by = c(
                "case",
                "subcase_id",
                "method"
            )
        )

    if (nrow(missing_timing_cells) > 0L) {
        fallback_time <- max(
            timing_check$max_time,
            na.rm = TRUE
        )

        if (!is.finite(fallback_time)) {
            stop(
                "No fitted-SITAR feasibility timings are available at all.",
                call. = FALSE
            )
        }

        warning(
            "No fitted-SITAR feasibility timing is available for ",
            nrow(missing_timing_cells),
            " method--subcase cells. Using the largest observed fitted-SITAR ",
            "feasibility time, ",
            signif(fallback_time, 4),
            " minutes, for scheduling those cells.",
            call. = FALSE
        )

        timing_fallback_rows <- missing_timing_cells |>
            dplyr::mutate(
                seed = NA_integer_,
                missing_from_simstudy = FALSE,
                sim_status = "ok",
                fit_status = "ok",
                pred_status = "ok",
                gp_status = "ok",
                time = fallback_time,
                time_fit = NA_real_
            )

        timing_status <- dplyr::bind_rows(
            timing_status,
            timing_fallback_rows
        )

        readr::write_csv(
            missing_timing_cells,
            file.path(
                diagnostic_dir,
                "sitar_method_missing_timing_cells.csv"
            )
        )
    }

    timing_status
}


if (stage == "prepare") {
    settings_file <- here::here(
        "instances",
        "feasibility_analysis",
        "output",
        "chosen_fit_cutoff.csv"
    )

    feasibility_status_file <- here::here(
        "instances",
        "feasibility_analysis",
        "output",
        "complete_status_all_cases.rds"
    )

    if (!file.exists(settings_file)) {
        stop(
            "Missing ",
            settings_file,
            ". Run the feasibility finalisation first.",
            call. = FALSE
        )
    }

    if (!file.exists(feasibility_status_file)) {
        stop(
            "Missing ",
            feasibility_status_file,
            ". Run the feasibility finalisation first.",
            call. = FALSE
        )
    }

    feasibility_settings <- readr::read_csv(
        settings_file,
        show_col_types = FALSE
    )

    n_full_seeds <- as.integer(
        get_unique_setting(
            feasibility_settings,
            "n_full_seeds"
        )
    )

    scheduling_buffer <- as.numeric(
        get_unique_setting(
            feasibility_settings,
            "scheduling_buffer",
            default = scheduling_buffer_if_missing
        )
    )

    if (
        !is.finite(n_full_seeds) ||
        n_full_seeds <= 0L
    ) {
        stop(
            "n_full_seeds must be a positive integer.",
            call. = FALSE
        )
    }

    sitar_extra_cells <- applicable_method_cells |>
        dplyr::filter(
            case == case_to_run,
            method == "SITAR"
        ) |>
        dplyr::distinct()

    if (nrow(sitar_extra_cells) == 0L) {
        stop(
            "No applicable fitted-SITAR cells were found.",
            call. = FALSE
        )
    }

    setup_unallocated <- tidyr::crossing(
        sitar_extra_cells,
        seed = seq_len(n_full_seeds)
    ) |>
        dplyr::arrange(
            subcase_id,
            seed
        )

    complete_status_feasibility <- readRDS(
        feasibility_status_file
    ) |>
        dplyr::filter(
            case == case_to_run,
            method == "SITAR"
        )

    if (nrow(complete_status_feasibility) == 0L) {
        stop(
            "No fitted-SITAR feasibility status rows were found.",
            call. = FALSE
        )
    }

    complete_status_feasibility_for_timing <- make_sitar_timing_status(
        complete_status_feasibility =
            complete_status_feasibility,
        setup_unallocated = setup_unallocated,
        diagnostic_dir = diagnostic_dir
    )

    setup_extra <- find_iridis_setup_from_feasibility_checkpointed(
        setup = setup_unallocated,
        simstudy_feasibility =
            complete_status_feasibility_for_timing,
        time_each = initial_walltime_minutes,
        buffer = scheduling_buffer,
        back = binpacking_back,
        pack_by = "method"
    )

    if (dir.exists(iridis_dir)) {
        stop(
            "Directory ",
            iridis_dir,
            " already exists. Delete or rename it before preparing ",
            "a new fitted-SITAR run.",
            call. = FALSE
        )
    }

    readr::write_csv(
        setup_extra,
        file.path(
            diagnostic_dir,
            "setup_sitar_method.csv"
        )
    )

    saveRDS(
        setup_extra,
        file.path(
            diagnostic_dir,
            "setup_sitar_method.rds"
        )
    )

    setup_written <- write_simstudy_iridis_checkpointed(
        setup = setup_extra,
        subcases = subcases,
        simulation_fun = run_simstudy_each,
        simulation_fun_name = "run_simstudy_each",
        time_each = initial_walltime_minutes,
        mem = mem_gb,
        iridis_dir = iridis_dir
    )

    readr::write_csv(
        setup_written,
        file.path(
            diagnostic_dir,
            "setup_sitar_method_written.csv"
        )
    )

    message(
        "\nCreated the fitted-SITAR run directory:\n  ",
        iridis_dir,
        "\n\nSubmit with:\n\n  cd ",
        iridis_dir,
        "\n  sbatch simstudy_run.slurm\n\n",
        "After it has finished, rerun this script with stage <- \"after_4h\"."
    )
}


if (stage == "after_4h") {
    setup <- read_setup(
        iridis_dir
    )

    diagnostics <- write_diagnostics(
        prefix = "after_4h",
        setup = setup,
        dir = iridis_dir,
        time_limit_minutes = initial_walltime_minutes
    )

    missing_process_ids <- unique(
        diagnostics$missing_processes$process_id
    )

    oom_process_ids <- intersect(
        missing_process_ids,
        unique(diagnostics$oom_processes$process_id)
    )

    timeout_process_ids <- intersect(
        missing_process_ids,
        unique(diagnostics$timeout_processes$process_id)
    )

    update_process_ids <- setdiff(
        missing_process_ids,
        oom_process_ids
    )

    update_plan <- tibble::tibble(
        process_id = missing_process_ids
    ) |>
        dplyr::mutate(
            initial_oom = process_id %in% oom_process_ids,
            initial_timeout =
                process_id %in% timeout_process_ids,
            action = dplyr::case_when(
                initial_oom ~
                    "defer_to_isolated_recovery",
                TRUE ~
                    "rerun_packed_for_8_hours"
            )
        )

    readr::write_csv(
        update_plan,
        file.path(
            diagnostic_dir,
            "after_4h_update_plan.csv"
        )
    )

    if (length(missing_process_ids) == 0L) {
        message(
            "All fitted-SITAR runs have checkpoints. ",
            "Proceed to incorporate_sitar_method.R."
        )
    } else if (length(update_process_ids) == 0L) {
        message(
            "All missing fitted-SITAR processes were OOM. ",
            "Rerun this script with stage <- \"after_8h\" to create ",
            "isolated recovery jobs."
        )
    } else {
        write_update_script_iridis_checkpointed(
            process_ids = update_process_ids,
            time_each = extended_walltime_minutes,
            mem = mem_gb,
            iridis_dir = iridis_dir,
            script_name = "simstudy_update_extended.slurm"
        )

        message(
            "\nCreated an 8-hour update script. Submit with:\n\n  cd ",
            iridis_dir,
            "\n  sbatch simstudy_update_extended.slurm\n\n",
            "After it has finished, rerun this script with stage <- \"after_8h\"."
        )
    }
}


if (stage == "after_8h") {
    setup <- read_setup(
        iridis_dir
    )

    diagnostics <- write_diagnostics(
        prefix = "after_8h",
        setup = setup,
        dir = iridis_dir,
        time_limit_minutes = extended_walltime_minutes
    )

    if (nrow(diagnostics$missing_runs) == 0L) {
        message(
            "All fitted-SITAR runs now have checkpoints. ",
            "Proceed to incorporate_sitar_method.R."
        )
    } else {
        if (dir.exists(recovery_iridis_dir)) {
            stop(
                "Directory ",
                recovery_iridis_dir,
                " already exists. Delete or rename it before creating ",
                "a new fitted-SITAR recovery run.",
                call. = FALSE
            )
        }

        completed_run_ids <-
            find_completed_run_ids_iridis_checkpointed(
                iridis_dir,
                validate = TRUE
            )

        recovery_setup <- make_isolated_recovery_setup(
            setup = setup,
            completed_run_ids = completed_run_ids,
            isolate = TRUE
        )

        write_simstudy_iridis_checkpointed(
            setup = recovery_setup,
            subcases = subcases,
            simulation_fun = run_simstudy_each,
            simulation_fun_name = "run_simstudy_each",
            time_each = extended_walltime_minutes,
            mem = mem_gb,
            iridis_dir = recovery_iridis_dir
        )

        readr::write_csv(
            recovery_setup,
            file.path(
                diagnostic_dir,
                "isolated_recovery_setup.csv"
            )
        )

        saveRDS(
            recovery_setup,
            file.path(
                diagnostic_dir,
                "isolated_recovery_setup.rds"
            )
        )

        message(
            "\nCreated one-run-per-process fitted-SITAR recovery jobs. ",
            "Submit with:\n\n  cd ",
            recovery_iridis_dir,
            "\n  sbatch simstudy_run.slurm\n\n",
            "After it has finished, rerun this script with ",
            "stage <- \"after_recovery\"."
        )
    }
}


if (stage == "after_recovery") {
    if (!dir.exists(recovery_iridis_dir)) {
        message(
            "No fitted-SITAR recovery directory exists. ",
            "Proceed to incorporate_sitar_method.R."
        )
    } else {
        recovery_setup <- read_setup(
            recovery_iridis_dir
        )

        diagnostics <- write_diagnostics(
            prefix = "after_recovery",
            setup = recovery_setup,
            dir = recovery_iridis_dir,
            time_limit_minutes = extended_walltime_minutes
        )

        resource_failure_process_ids <- union(
            unique(diagnostics$oom_processes$process_id),
            unique(diagnostics$timeout_processes$process_id)
        )

        unknown_missing_process_ids <- setdiff(
            unique(
                diagnostics$missing_processes$process_id
            ),
            resource_failure_process_ids
        )

        retry_plan <- tibble::tibble(
            process_id = unique(
                diagnostics$missing_processes$process_id
            )
        ) |>
            dplyr::mutate(
                definitive_oom =
                    process_id %in%
                    unique(diagnostics$oom_processes$process_id),
                definitive_timeout =
                    process_id %in%
                    unique(diagnostics$timeout_processes$process_id),
                action = dplyr::case_when(
                    definitive_oom ~
                        "record_run_level_oom",
                    definitive_timeout ~
                        "record_run_level_timeout",
                    TRUE ~
                        "retry_isolated_for_8_hours"
                )
            )

        readr::write_csv(
            retry_plan,
            file.path(
                diagnostic_dir,
                "after_recovery_retry_plan.csv"
            )
        )

        if (nrow(diagnostics$missing_runs) == 0L) {
            message(
                "All fitted-SITAR recovery runs completed. ",
                "Proceed to incorporate_sitar_method.R."
            )
        } else if (length(unknown_missing_process_ids) == 0L) {
            message(
                "All remaining fitted-SITAR missing recovery runs have ",
                "definitive OOM or timeout classifications. Proceed to ",
                "incorporate_sitar_method.R."
            )
        } else {
            write_update_script_iridis_checkpointed(
                process_ids = unknown_missing_process_ids,
                time_each = extended_walltime_minutes,
                mem = mem_gb,
                iridis_dir = recovery_iridis_dir,
                script_name = "simstudy_retry_unknown_extended.slurm"
            )

            message(
                "\nCreated an 8-hour retry for unexplained missing ",
                "fitted-SITAR isolated runs. Submit with:\n\n  cd ",
                recovery_iridis_dir,
                "\n  sbatch simstudy_retry_unknown_extended.slurm\n\n",
                "After it has finished, rerun stage <- \"after_recovery\"."
            )
        }
    }
}
