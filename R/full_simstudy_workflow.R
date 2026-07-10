## Workflow helpers for the checkpointed full simulation studies.
##
## Low-level checkpointing and recovery are implemented in
## R/iridis_checkpointed.R.  This file coordinates the common staged
## workflow used by the 2re, 3re and SITAR full simulation studies.


.full_simstudy_require_files <- function(files, description) {
    missing_files <- files[!file.exists(files)]

    if (length(missing_files) > 0L) {
        stop(
            description,
            " cannot proceed because these files are missing:\n",
            paste(missing_files, collapse = "\n"),
            call. = FALSE
        )
    }

    invisible(files)
}


.full_simstudy_load_rda_object <- function(file, object_name) {
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


.full_simstudy_get_unique_setting <- function(
    x,
    column,
    description
) {
    if (!(column %in% names(x))) {
        stop(
            description,
            " is missing column ",
            column,
            ".",
            call. = FALSE
        )
    }

    values <- unique(x[[column]])

    if (
        length(values) != 1L ||
        length(values[[1]]) != 1L ||
        is.na(values[[1]])
    ) {
        stop(
            description,
            " must contain exactly one non-missing value of ",
            column,
            ".",
            call. = FALSE
        )
    }

    values[[1]]
}


.full_simstudy_context <- function(
    case_name,
    subcases,
    applicable_method_cells,
    simulation_fun,
    simulation_fun_name,
    initial_walltime_minutes,
    extended_walltime_minutes,
    mem_gb,
    binpacking_back,
    pack_by,
    project_root
) {
    if (
        length(case_name) != 1L ||
        is.na(case_name) ||
        !nzchar(case_name)
    ) {
        stop(
            "case_name must be one non-empty string.",
            call. = FALSE
        )
    }

    if (
        length(simulation_fun_name) != 1L ||
        is.na(simulation_fun_name) ||
        !grepl(
            "^[A-Za-z.][A-Za-z0-9._]*$",
            simulation_fun_name
        )
    ) {
        stop(
            "simulation_fun_name must be the name of a package function.",
            call. = FALSE
        )
    }

    if (!is.function(simulation_fun)) {
        stop(
            "simulation_fun must be a function.",
            call. = FALSE
        )
    }

    if (
        length(initial_walltime_minutes) != 1L ||
        !is.finite(initial_walltime_minutes) ||
        initial_walltime_minutes <= 0
    ) {
        stop(
            "initial_walltime_minutes must be positive.",
            call. = FALSE
        )
    }

    if (
        length(extended_walltime_minutes) != 1L ||
        !is.finite(extended_walltime_minutes) ||
        extended_walltime_minutes <= initial_walltime_minutes
    ) {
        stop(
            "extended_walltime_minutes must exceed ",
            "initial_walltime_minutes.",
            call. = FALSE
        )
    }

    if (
        length(mem_gb) != 1L ||
        !is.finite(mem_gb) ||
        mem_gb <= 0
    ) {
        stop(
            "mem_gb must be positive.",
            call. = FALSE
        )
    }

    required_cells <- c(
        "case",
        "subcase_id",
        "method"
    )

    missing_cell_columns <- setdiff(
        required_cells,
        names(applicable_method_cells)
    )

    if (length(missing_cell_columns) > 0L) {
        stop(
            "applicable_method_cells is missing: ",
            paste(missing_cell_columns, collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    cases_in_design <- unique(
        applicable_method_cells$case
    )

    if (
        length(cases_in_design) != 1L ||
        !identical(as.character(cases_in_design), as.character(case_name))
    ) {
        stop(
            "applicable_method_cells does not describe only case ",
            case_name,
            ".",
            call. = FALSE
        )
    }

    full_dir <- file.path(
        project_root,
        "instances",
        case_name,
        "full"
    )

    output_dir <- file.path(
        full_dir,
        "output"
    )

    list(
        case_name = case_name,
        subcases = subcases,
        applicable_method_cells =
            tibble::as_tibble(applicable_method_cells),
        simulation_fun = simulation_fun,
        simulation_fun_name = simulation_fun_name,

        initial_walltime_minutes =
            initial_walltime_minutes,
        extended_walltime_minutes =
            extended_walltime_minutes,
        mem_gb = mem_gb,
        binpacking_back = binpacking_back,
        pack_by = pack_by,

        project_root = project_root,
        full_dir = full_dir,
        output_dir = output_dir,

        primary_iridis_dir = file.path(
            full_dir,
            "iridis"
        ),
        recovery_iridis_dir = file.path(
            full_dir,
            "iridis_recovery"
        ),

        decisions_file = file.path(
            project_root,
            "instances",
            "feasibility_analysis",
            "output",
            "feasibility_decisions.csv"
        ),
        feasibility_status_file = file.path(
            project_root,
            "instances",
            "feasibility_analysis",
            "output",
            "complete_status_all_cases.rds"
        ),
        feasibility_settings_file = file.path(
            project_root,
            "instances",
            "feasibility_analysis",
            "output",
            "chosen_fit_cutoff.csv"
        )
    )
}


.full_simstudy_submission_message <- function(
    iridis_dir,
    script_name
) {
    paste0(
        "\nRun:\n\n",
        "  cd ",
        normalizePath(
            iridis_dir,
            mustWork = FALSE
        ),
        "\n",
        "  sbatch ",
        script_name,
        "\n"
    )
}


.full_simstudy_write_process_diagnostics <- function(
    context,
    diagnostics,
    prefix
) {
    readr::write_csv(
        diagnostics$oom,
        file.path(
            context$output_dir,
            paste0(prefix, "_oom_processes.csv")
        )
    )

    readr::write_csv(
        diagnostics$timeout,
        file.path(
            context$output_dir,
            paste0(prefix, "_timeout_processes.csv")
        )
    )

    readr::write_csv(
        diagnostics$missing_processes,
        file.path(
            context$output_dir,
            paste0(
                prefix,
                "_processes_with_missing_runs.csv"
            )
        )
    )

    readr::write_csv(
        diagnostics$missing_runs,
        file.path(
            context$output_dir,
            paste0(prefix, "_missing_runs.csv")
        )
    )

    readr::write_csv(
        diagnostics$active_runs,
        file.path(
            context$output_dir,
            paste0(prefix, "_active_runs_at_failure.csv")
        )
    )

    invisible(diagnostics)
}


.full_simstudy_read_setup <- function(iridis_dir) {
    setup_file <- file.path(
        iridis_dir,
        "setup.Rda"
    )

    .full_simstudy_require_files(
        setup_file,
        "Reading the Iridis setup"
    )

    .full_simstudy_load_rda_object(
        setup_file,
        "setup"
    )
}


.full_simstudy_process_diagnostics <- function(
    iridis_dir,
    setup,
    mem_gb,
    time_limit_minutes
) {
    missing_runs <- find_missing_runs_iridis_checkpointed(
        iridis_dir,
        setup = setup
    )

    missing_processes <-
        find_processes_with_missing_runs_iridis_checkpointed(
            iridis_dir,
            setup = setup
        )

    oom <- find_oom_ids_iridis(
        iridis_dir,
        min_mem_gb = mem_gb,
        mem_limit_if_missing_gb = mem_gb
    )

    timeout <- find_timeout_ids_iridis(
        iridis_dir,
        min_time_minutes = time_limit_minutes,
        time_limit_if_missing_minutes =
            time_limit_minutes
    )

    active_runs <-
        find_active_runs_at_failure_iridis_checkpointed(
            iridis_dir,
            process_ids = unique(
                missing_processes$process_id
            )
        )

    list(
        missing_runs = missing_runs,
        missing_processes = missing_processes,
        oom = oom,
        timeout = timeout,
        active_runs = active_runs
    )
}


.full_simstudy_write_checkpointed_run <- function(
    context,
    setup,
    iridis_dir,
    time_each
) {
    write_simstudy_iridis_checkpointed(
        setup = setup,
        subcases = context$subcases,
        simulation_fun = context$simulation_fun,
        simulation_fun_name =
            context$simulation_fun_name,
        time_each = time_each,
        mem = context$mem_gb,
        iridis_dir = iridis_dir
    )
}


.full_simstudy_prepare <- function(context) {
    .full_simstudy_require_files(
        c(
            context$decisions_file,
            context$feasibility_status_file,
            context$feasibility_settings_file
        ),
        "The prepare stage"
    )

    decisions <- readr::read_csv(
        context$decisions_file,
        show_col_types = FALSE
    )

    feasibility_settings <- readr::read_csv(
        context$feasibility_settings_file,
        show_col_types = FALSE
    )

    n_full_seeds <- as.integer(
        .full_simstudy_get_unique_setting(
            feasibility_settings,
            "n_full_seeds",
            "chosen_fit_cutoff.csv"
        )
    )

    scheduling_buffer <- as.numeric(
        .full_simstudy_get_unique_setting(
            feasibility_settings,
            "scheduling_buffer",
            "chosen_fit_cutoff.csv"
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

    if (
        !is.finite(scheduling_buffer) ||
        scheduling_buffer <= 0
    ) {
        stop(
            "scheduling_buffer must be positive.",
            call. = FALSE
        )
    }

    decisions_case <- decisions |>
        dplyr::filter(
            case == context$case_name
        )

    if (nrow(decisions_case) == 0L) {
        stop(
            "No feasibility decisions were found for case ",
            context$case_name,
            ".",
            call. = FALSE
        )
    }

    unresolved_cells <- decisions_case |>
        dplyr::filter(
            decision %in% c(
                "borderline",
                "manual_review"
            )
        )

    if (nrow(unresolved_cells) > 0L) {
        print(unresolved_cells)

        stop(
            "Some feasibility decisions remain unresolved.",
            call. = FALSE
        )
    }

    cells_full <- decisions_case |>
        dplyr::filter(
            decision == "include"
        ) |>
        dplyr::select(
            case,
            subcase_id,
            method
        ) |>
        dplyr::distinct() |>
        dplyr::arrange(
            method,
            subcase_id
        )

    if (nrow(cells_full) == 0L) {
        stop(
            "No method--subcase cells are included for case ",
            context$case_name,
            ".",
            call. = FALSE
        )
    }

    invalid_cells <- cells_full |>
        dplyr::anti_join(
            context$applicable_method_cells,
            by = c(
                "case",
                "subcase_id",
                "method"
            )
        )

    if (nrow(invalid_cells) > 0L) {
        print(invalid_cells)

        stop(
            "The final decisions contain cells which are not ",
            "listed in applicable_method_cells.",
            call. = FALSE
        )
    }

    setup_unallocated <- tidyr::crossing(
        cells_full,
        seed = seq_len(n_full_seeds)
    )

    complete_status_feasibility <- readRDS(
        context$feasibility_status_file
    ) |>
        dplyr::filter(
            case == context$case_name
        ) |>
        dplyr::semi_join(
            cells_full,
            by = c(
                "case",
                "subcase_id",
                "method"
            )
        )

    setup_full <-
        find_iridis_setup_from_feasibility_checkpointed(
            setup = setup_unallocated,
            simstudy_feasibility =
                complete_status_feasibility,
            time_each =
                context$initial_walltime_minutes,
            buffer = scheduling_buffer,
            back = context$binpacking_back,
            pack_by = context$pack_by
        )

    if (dir.exists(context$primary_iridis_dir)) {
        stop(
            "Directory ",
            context$primary_iridis_dir,
            " already exists. Delete or rename it before ",
            "preparing a new full simulation run.",
            call. = FALSE
        )
    }

    dir.create(
        context$output_dir,
        recursive = TRUE,
        showWarnings = FALSE
    )

    readr::write_csv(
        cells_full,
        file.path(
            context$output_dir,
            "full_run_cells.csv"
        )
    )

    readr::write_csv(
        setup_full,
        file.path(
            context$output_dir,
            "setup_full.csv"
        )
    )

    saveRDS(
        setup_full,
        file.path(
            context$output_dir,
            "setup_full.rds"
        )
    )

    run_settings <- tibble::tibble(
        case = context$case_name,
        n_full_seeds = n_full_seeds,
        scheduling_buffer = scheduling_buffer,
        initial_walltime_minutes =
            context$initial_walltime_minutes,
        extended_walltime_minutes =
            context$extended_walltime_minutes,
        mem_gb = context$mem_gb,
        pack_by = paste(
            context$pack_by,
            collapse = ","
        )
    )

    readr::write_csv(
        run_settings,
        file.path(
            context$output_dir,
            "full_run_settings.csv"
        )
    )

    .full_simstudy_write_checkpointed_run(
        context = context,
        setup = setup_full,
        iridis_dir = context$primary_iridis_dir,
        time_each =
            context$initial_walltime_minutes
    )

    message(
        "\nCreated the initial checkpointed full run for ",
        context$case_name,
        ".",
        .full_simstudy_submission_message(
            context$primary_iridis_dir,
            "simstudy_run.slurm"
        ),
        "\nAfter it has finished, use stage = \"after_4h\"."
    )

    invisible(setup_full)
}


.full_simstudy_after_4h <- function(context) {
    setup <- .full_simstudy_read_setup(
        context$primary_iridis_dir
    )

    diagnostics <- .full_simstudy_process_diagnostics(
        iridis_dir = context$primary_iridis_dir,
        setup = setup,
        mem_gb = context$mem_gb,
        time_limit_minutes =
            context$initial_walltime_minutes
    )

    .full_simstudy_write_process_diagnostics(
        context,
        diagnostics,
        prefix = "after_4h"
    )

    missing_process_ids <- unique(
        diagnostics$missing_processes$process_id
    )

    oom_process_ids <- intersect(
        missing_process_ids,
        unique(diagnostics$oom$process_id)
    )

    timeout_process_ids <- intersect(
        missing_process_ids,
        unique(diagnostics$timeout$process_id)
    )

    ## Packed processes with missing runs are retried for eight hours,
    ## except those whose most recent attempt ended in OOM.  The latter
    ## are deferred to one-run-per-process recovery.
    update_process_ids <- setdiff(
        missing_process_ids,
        oom_process_ids
    )

    update_plan <- tibble::tibble(
        process_id = missing_process_ids
    ) |>
        dplyr::mutate(
            initial_oom =
                process_id %in% oom_process_ids,
            initial_timeout =
                process_id %in% timeout_process_ids,
            action = dplyr::case_when(
                initial_oom ~
                    "defer_to_isolated_recovery",
                TRUE ~
                    "rerun_packed_for_extended_time"
            )
        )

    readr::write_csv(
        update_plan,
        file.path(
            context$output_dir,
            "after_4h_update_plan.csv"
        )
    )

    if (length(missing_process_ids) == 0L) {
        message(
            "All runs have valid checkpoints. ",
            "No extended packed update is required; ",
            "use stage = \"finalise\"."
        )

        return(invisible(update_plan))
    }

    if (length(update_process_ids) == 0L) {
        message(
            "All processes with missing runs were OOM. ",
            "No packed update was created; ",
            "use stage = \"after_8h\" to create isolated recovery."
        )

        return(invisible(update_plan))
    }

    write_update_script_iridis_checkpointed(
        process_ids = update_process_ids,
        time_each =
            context$extended_walltime_minutes,
        mem = context$mem_gb,
        iridis_dir = context$primary_iridis_dir,
        script_name =
            "simstudy_update_extended.slurm"
    )

    message(
        "\nCreated the extended packed update for ",
        length(update_process_ids),
        " process IDs.",
        .full_simstudy_submission_message(
            context$primary_iridis_dir,
            "simstudy_update_extended.slurm"
        ),
        "\nCompleted run checkpoints will be skipped. ",
        "After it has finished, use stage = \"after_8h\"."
    )

    invisible(update_plan)
}


.full_simstudy_after_8h <- function(context) {
    setup <- .full_simstudy_read_setup(
        context$primary_iridis_dir
    )

    diagnostics <- .full_simstudy_process_diagnostics(
        iridis_dir = context$primary_iridis_dir,
        setup = setup,
        mem_gb = context$mem_gb,
        time_limit_minutes =
            context$extended_walltime_minutes
    )

    .full_simstudy_write_process_diagnostics(
        context,
        diagnostics,
        prefix = "after_8h"
    )

    if (nrow(diagnostics$missing_runs) == 0L) {
        message(
            "All runs now have valid checkpoints. ",
            "No isolated recovery is required; ",
            "use stage = \"finalise\"."
        )

        return(invisible(diagnostics$missing_runs))
    }

    if (dir.exists(context$recovery_iridis_dir)) {
        stop(
            "Directory ",
            context$recovery_iridis_dir,
            " already exists. Delete or rename it before creating ",
            "a new isolated recovery run.",
            call. = FALSE
        )
    }

    completed_run_ids <-
        find_completed_run_ids_iridis_checkpointed(
            context$primary_iridis_dir,
            validate = TRUE
        )

    recovery_setup <- make_isolated_recovery_setup(
        setup = setup,
        completed_run_ids = completed_run_ids,
        isolate = TRUE
    )

    if (nrow(recovery_setup) == 0L) {
        stop(
            "Missing runs were reported, but no recovery setup ",
            "could be constructed.",
            call. = FALSE
        )
    }

    .full_simstudy_write_checkpointed_run(
        context = context,
        setup = recovery_setup,
        iridis_dir = context$recovery_iridis_dir,
        time_each =
            context$extended_walltime_minutes
    )

    readr::write_csv(
        recovery_setup,
        file.path(
            context$output_dir,
            "isolated_recovery_setup.csv"
        )
    )

    saveRDS(
        recovery_setup,
        file.path(
            context$output_dir,
            "isolated_recovery_setup.rds"
        )
    )

    message(
        "\nCreated one-run-per-process recovery jobs for ",
        nrow(recovery_setup),
        " incomplete runs.",
        .full_simstudy_submission_message(
            context$recovery_iridis_dir,
            "simstudy_run.slurm"
        ),
        "\nAfter recovery has finished, ",
        "use stage = \"after_recovery\"."
    )

    invisible(recovery_setup)
}


.full_simstudy_after_recovery <- function(context) {
    if (!dir.exists(context$recovery_iridis_dir)) {
        message(
            "No isolated recovery directory exists. ",
            "There is no recovery run to inspect; ",
            "use stage = \"finalise\"."
        )

        return(invisible(NULL))
    }

    recovery_setup <- .full_simstudy_read_setup(
        context$recovery_iridis_dir
    )

    diagnostics <- .full_simstudy_process_diagnostics(
        iridis_dir = context$recovery_iridis_dir,
        setup = recovery_setup,
        mem_gb = context$mem_gb,
        time_limit_minutes =
            context$extended_walltime_minutes
    )

    .full_simstudy_write_process_diagnostics(
        context,
        diagnostics,
        prefix = "after_recovery"
    )

    resource_failure_process_ids <- union(
        unique(diagnostics$oom$process_id),
        unique(diagnostics$timeout$process_id)
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
                unique(diagnostics$oom$process_id),
            definitive_timeout =
                process_id %in%
                unique(diagnostics$timeout$process_id),
            action = dplyr::case_when(
                definitive_oom ~
                    "record_run_level_oom",
                definitive_timeout ~
                    "record_run_level_timeout",
                TRUE ~
                    "retry_isolated_for_extended_time"
            )
        )

    readr::write_csv(
        retry_plan,
        file.path(
            context$output_dir,
            "after_recovery_retry_plan.csv"
        )
    )

    if (nrow(diagnostics$missing_runs) == 0L) {
        message(
            "All isolated recovery runs completed; ",
            "use stage = \"finalise\"."
        )

        return(invisible(retry_plan))
    }

    if (length(unknown_missing_process_ids) == 0L) {
        message(
            "All remaining missing isolated runs have definitive ",
            "OOM or timeout classifications; ",
            "use stage = \"finalise\"."
        )

        return(invisible(retry_plan))
    }

    write_update_script_iridis_checkpointed(
        process_ids = unknown_missing_process_ids,
        time_each =
            context$extended_walltime_minutes,
        mem = context$mem_gb,
        iridis_dir = context$recovery_iridis_dir,
        script_name =
            "simstudy_retry_unknown_extended.slurm"
    )

    message(
        "\nCreated an extended retry for ",
        length(unknown_missing_process_ids),
        " unexplained missing isolated runs.",
        .full_simstudy_submission_message(
            context$recovery_iridis_dir,
            "simstudy_retry_unknown_extended.slurm"
        ),
        "\nAfter it has finished, rerun ",
        "stage = \"after_recovery\"."
    )

    invisible(retry_plan)
}


.full_simstudy_finalise <- function(context) {
    setup <- .full_simstudy_read_setup(
        context$primary_iridis_dir
    )

    checkpoint_dirs <- context$primary_iridis_dir

    if (dir.exists(context$recovery_iridis_dir)) {
        checkpoint_dirs <- c(
            checkpoint_dirs,
            context$recovery_iridis_dir
        )
    }

    ## Do not read the full checkpointed simulation outputs here.  Those
    ## contain large pred_data/GP list columns and can exceed memory when all
    ## runs are combined.  Finalisation only needs scalar status fields plus an
    ## index of checkpoint files for the later analysis stage.
    checkpoint_index <- checkpointed_result_file_index(
        checkpoint_dirs,
        prefer = "last"
    )

    run_status_results <- read_checkpointed_run_status(
        checkpoint_index = checkpoint_index,
        verbose = TRUE
    )

    if (dir.exists(context$recovery_iridis_dir)) {
        oom_run_ids <-
            find_oom_run_ids_iridis_checkpointed(
                context$recovery_iridis_dir,
                min_mem_gb = context$mem_gb,
                mem_limit_if_missing_gb =
                    context$mem_gb
            )

        timeout_run_ids <-
            find_timeout_run_ids_iridis_checkpointed(
                context$recovery_iridis_dir,
                min_time_minutes =
                    context$extended_walltime_minutes,
                time_limit_if_missing_minutes =
                    context$extended_walltime_minutes
            )

        ## Give OOM precedence if logs happened to classify a run
        ## both ways.
        timeout_run_ids <- setdiff(
            timeout_run_ids,
            oom_run_ids
        )
    } else {
        oom_run_ids <- character()
        timeout_run_ids <- character()
    }

    completed_run_ids <- unique(
        as.character(run_status_results$run_id)
    )

    resolved_run_ids <- union(
        completed_run_ids,
        union(
            oom_run_ids,
            timeout_run_ids
        )
    )

    unresolved_runs <- setup |>
        dplyr::filter(
            !as.character(run_id) %in%
                resolved_run_ids
        )

    readr::write_csv(
        unresolved_runs,
        file.path(
            context$output_dir,
            "unresolved_runs.csv"
        )
    )

    if (nrow(unresolved_runs) > 0L) {
        print(
            unresolved_runs |>
                dplyr::select(
                    run_id,
                    process_id,
                    case,
                    subcase_id,
                    method,
                    seed
                )
        )

        stop(
            nrow(unresolved_runs),
            " runs are neither checkpointed nor definitively ",
            "classified as OOM or timeout. Inspect recovery logs ",
            "and use stage = \"after_recovery\" before finalising.",
            call. = FALSE
        )
    }

    complete_status_full <-
        make_complete_run_status_checkpointed(
            setup = setup,
            simstudy = run_status_results,
            oom_run_ids = oom_run_ids,
            timeout_run_ids = timeout_run_ids
        )

    full_status_summary <-
        summarise_complete_run_status(
            complete_status_full
        ) |>
        dplyr::left_join(
            context$subcases,
            by = "subcase_id"
        ) |>
        dplyr::relocate(
            case,
            subcase_id,
            n_clusters,
            n_obs_per_cluster,
            method
        )

    oom_runs <- setup |>
        dplyr::filter(
            as.character(run_id) %in%
                oom_run_ids
        ) |>
        dplyr::arrange(
            method,
            subcase_id,
            seed
        )

    timeout_runs <- setup |>
        dplyr::filter(
            as.character(run_id) %in%
                timeout_run_ids
        ) |>
        dplyr::arrange(
            method,
            subcase_id,
            seed
        )

    worker_error_runs <- run_status_results |>
        dplyr::filter(
            dplyr::coalesce(
                checkpoint_worker_status == "error",
                FALSE
            )
        ) |>
        dplyr::select(
            dplyr::any_of(
                c(
                    "run_id",
                    "process_id",
                    "case",
                    "subcase_id",
                    "method",
                    "seed",
                    "checkpoint_worker_error"
                )
            )
        )

    ## Save the small objects needed by downstream scripts.  The full
    ## checkpointed simulation results remain on disk as one file per run; the
    ## analysis pipeline should read/summarise them one at a time using
    ## checkpointed_result_index.rds.
    saveRDS(
        checkpoint_index,
        file.path(
            context$output_dir,
            "checkpointed_result_index.rds"
        )
    )

    readr::write_csv(
        checkpoint_index,
        file.path(
            context$output_dir,
            "checkpointed_result_index.csv"
        )
    )

    saveRDS(
        run_status_results,
        file.path(
            context$output_dir,
            "run_status_from_checkpoints.rds"
        )
    )

    saveRDS(
        complete_status_full,
        file.path(
            context$output_dir,
            "complete_status_full.rds"
        )
    )

    readr::write_csv(
        full_status_summary,
        file.path(
            context$output_dir,
            "full_status_summary.csv"
        )
    )

    readr::write_csv(
        oom_runs,
        file.path(
            context$output_dir,
            "oom_runs.csv"
        )
    )

    readr::write_csv(
        timeout_runs,
        file.path(
            context$output_dir,
            "timeout_runs.csv"
        )
    )

    readr::write_csv(
        worker_error_runs,
        file.path(
            context$output_dir,
            "worker_error_runs.csv"
        )
    )

    final_counts <- tibble::tibble(
        case = context$case_name,
        n_expected = nrow(setup),
        n_checkpointed =
            length(completed_run_ids),
        n_oom = length(oom_run_ids),
        n_timeout = length(timeout_run_ids),
        n_worker_error =
            nrow(worker_error_runs),
        n_unresolved = nrow(unresolved_runs)
    )

    readr::write_csv(
        final_counts,
        file.path(
            context$output_dir,
            "full_run_final_counts.csv"
        )
    )

    print(final_counts)

    if (length(oom_run_ids) > 0L) {
        warning(
            length(oom_run_ids),
            " isolated runs exceeded the ",
            context$mem_gb,
            " GB memory limit.",
            call. = FALSE
        )
    }

    if (length(timeout_run_ids) > 0L) {
        warning(
            length(timeout_run_ids),
            " isolated runs exceeded the ",
            context$extended_walltime_minutes / 60,
            "-hour time limit.",
            call. = FALSE
        )
    }

    message(
        "Final checkpoint index and run-level statuses ",
        "were written to ",
        normalizePath(
            context$output_dir,
            mustWork = FALSE
        ),
        ". The full run outputs remain as one checkpoint file per run."
    )

    invisible(
        list(
            checkpoint_index = checkpoint_index,
            run_status = run_status_results,
            complete_status = complete_status_full,
            status_summary = full_status_summary,
            counts = final_counts
        )
    )
}


## Run one stage of a checkpointed full simulation study.
##
## The case-specific script supplies the design objects and changes only
## `stage` as the workflow progresses.
run_full_simstudy_stage <- function(
    case_name,
    subcases,
    applicable_method_cells,
    simulation_fun,
    simulation_fun_name = NULL,
    stage = c(
        "prepare",
        "after_4h",
        "after_8h",
        "after_recovery",
        "finalise"
    ),
    initial_walltime_minutes = 240,
    extended_walltime_minutes = 480,
    mem_gb = 8,
    binpacking_back = 0,
    pack_by = "method",
    project_root = here::here()
) {
    stage <- match.arg(stage)

    if (is.null(simulation_fun_name)) {
        simulation_fun_name <- deparse(
            substitute(simulation_fun)
        )
    }

    context <- .full_simstudy_context(
        case_name = case_name,
        subcases = subcases,
        applicable_method_cells =
            applicable_method_cells,
        simulation_fun = simulation_fun,
        simulation_fun_name =
            simulation_fun_name,
        initial_walltime_minutes =
            initial_walltime_minutes,
        extended_walltime_minutes =
            extended_walltime_minutes,
        mem_gb = mem_gb,
        binpacking_back = binpacking_back,
        pack_by = pack_by,
        project_root = project_root
    )

    dir.create(
        context$full_dir,
        recursive = TRUE,
        showWarnings = FALSE
    )

    dir.create(
        context$output_dir,
        recursive = TRUE,
        showWarnings = FALSE
    )

    switch(
        stage,
        prepare =
            .full_simstudy_prepare(context),
        after_4h =
            .full_simstudy_after_4h(context),
        after_8h =
            .full_simstudy_after_8h(context),
        after_recovery =
            .full_simstudy_after_recovery(context),
        finalise =
            .full_simstudy_finalise(context)
    )
}
