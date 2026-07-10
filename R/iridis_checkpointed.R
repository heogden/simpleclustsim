## Checkpointed Iridis helpers for the full simulation studies.
##
## The permanent unit of work is one method--subcase--seed run, identified by
## run_id.  process_id is only a temporary scheduling bin: several runs may be
## assigned to one Slurm array task, but each completed run is saved
## immediately to its own RDS file.
##
## This file deliberately does not replace the older helpers used by the
## feasibility jobs which may already be running.


.checkpoint_check_cols <- function(
    x,
    required,
    object_name = deparse(substitute(x))
) {
    missing <- setdiff(required, names(x))

    if (length(missing) > 0L) {
        stop(
            object_name,
            " is missing required column",
            if (length(missing) == 1L) " " else "s ",
            paste(missing, collapse = ", "),
            call. = FALSE
        )
    }

    invisible(x)
}


.checkpoint_safe_component <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- "NA"
    gsub("[^A-Za-z0-9._-]+", "_", x)
}


.checkpoint_run_filename <- function(run_id) {
    paste0("run-", .checkpoint_safe_component(run_id), ".rds")
}


.checkpoint_run_path <- function(iridis_dir, run_id) {
    file.path(
        iridis_dir,
        "storage",
        "runs",
        .checkpoint_run_filename(run_id)
    )
}


.checkpoint_progress_path <- function(iridis_dir, process_id) {
    file.path(
        iridis_dir,
        "storage",
        "progress",
        paste0("progress-", process_id, ".tsv")
    )
}


.checkpoint_scalar_chr <- function(x) {
    if (length(x) == 0L || is.null(x) || is.na(x[[1]])) {
        ""
    } else {
        as.character(x[[1]])
    }
}


.checkpoint_append_event <- function(
    iridis_dir,
    event,
    process_id,
    setup_row = NULL,
    message = ""
) {
    progress_file <- .checkpoint_progress_path(
        iridis_dir,
        process_id
    )

    dir.create(
        dirname(progress_file),
        recursive = TRUE,
        showWarnings = FALSE
    )

    get_value <- function(name) {
        if (is.null(setup_row) || !(name %in% names(setup_row))) {
            ""
        } else {
            .checkpoint_scalar_chr(setup_row[[name]])
        }
    }

    event_row <- data.frame(
        timestamp = format(
            Sys.time(),
            "%Y-%m-%d %H:%M:%S %z"
        ),
        event = as.character(event),
        process_id = as.character(process_id),
        run_id = get_value("run_id"),
        case = get_value("case"),
        subcase_id = get_value("subcase_id"),
        method = get_value("method"),
        seed = get_value("seed"),
        slurm_job_id = Sys.getenv("SLURM_JOB_ID", unset = ""),
        slurm_array_job_id = Sys.getenv(
            "SLURM_ARRAY_JOB_ID",
            unset = ""
        ),
        message = as.character(message),
        stringsAsFactors = FALSE
    )

    write.table(
        event_row,
        file = progress_file,
        append = file.exists(progress_file),
        sep = "\t",
        row.names = FALSE,
        col.names = !file.exists(progress_file),
        quote = TRUE,
        na = ""
    )

    invisible(event_row)
}


.checkpoint_output_is_valid <- function(file, expected_run_id = NULL) {
    if (!file.exists(file)) {
        return(FALSE)
    }

    out <- tryCatch(
        readRDS(file),
        error = function(e) NULL
    )

    if (!is.data.frame(out) || nrow(out) != 1L) {
        return(FALSE)
    }

    if (!is.null(expected_run_id)) {
        if (!("run_id" %in% names(out))) {
            return(FALSE)
        }

        if (!identical(
            as.character(out$run_id[[1]]),
            as.character(expected_run_id)
        )) {
            return(FALSE)
        }
    }

    TRUE
}


.checkpoint_worker_error_result <- function(setup_row, error_message) {
    key_names <- intersect(
        c("case", "subcase_id", "method", "seed"),
        names(setup_row)
    )

    out <- tibble::as_tibble(
        setup_row[, key_names, drop = FALSE]
    )

    out$data <- list(NULL)
    out$pred_data <- list(NULL)
    out$GP0 <- list(NULL)
    out$GP0_emp <- list(NULL)
    out$GP_hat <- list(list(Unavailable = NA_real_))
    out$k_hat <- NA_real_

    out$time <- NA_real_
    out$time_sim <- NA_real_
    out$time_fit <- NA_real_
    out$time_pred <- NA_real_
    out$time_gp <- NA_real_

    out$sim_status <- "worker_error"
    out$fit_status <- "skipped"
    out$pred_status <- "skipped"
    out$gp_status <- "skipped"
    out$k_status <- "skipped"

    out$sim_error <- as.character(error_message)
    out$fit_error <- NA_character_
    out$pred_error <- NA_character_
    out$gp_error <- NA_character_
    out$k_error <- NA_character_

    out$sim_warnings <- list(character())
    out$fit_warnings <- list(character())
    out$pred_warnings <- list(character())
    out$gp_warnings <- list(character())
    out$k_warnings <- list(character())

    out$n_sim_warnings <- 0L
    out$n_fit_warnings <- 0L
    out$n_pred_warnings <- 0L
    out$n_gp_warnings <- 0L
    out$n_k_warnings <- 0L

    out
}


## Add permanent run identifiers and within-process execution order.
##
## If run_id is absent, identifiers are constructed deterministically from
## case, subcase_id, method and seed.  Existing run_id values are retained,
## which is essential when constructing recovery setups.
add_checkpoint_run_ids <- function(
    setup,
    run_keys = c("case", "subcase_id", "method", "seed")
) {
    .checkpoint_check_cols(setup, run_keys)

    duplicated_key_rows <- duplicated(setup[, run_keys, drop = FALSE]) |
        duplicated(setup[, run_keys, drop = FALSE], fromLast = TRUE)

    if (any(duplicated_key_rows)) {
        duplicated_keys <- setup[
            duplicated_key_rows,
            run_keys,
            drop = FALSE
        ]
        print(duplicated_keys)
        stop(
            "setup contains duplicated run keys.",
            call. = FALSE
        )
    }

    out <- tibble::as_tibble(setup)

    if (!("run_id" %in% names(out))) {
        out <- out |>
            dplyr::mutate(
                run_id = paste(
                    .checkpoint_safe_component(case),
                    paste0(
                        "subcase-",
                        sprintf("%03d", as.integer(subcase_id))
                    ),
                    .checkpoint_safe_component(method),
                    paste0("seed-", .checkpoint_safe_component(seed)),
                    sep = "__"
                )
            )
    }

    if (anyNA(out$run_id) || any(out$run_id == "")) {
        stop(
            "run_id values must be non-missing and non-empty.",
            call. = FALSE
        )
    }

    if (any(as.character(out$run_id) != .checkpoint_safe_component(out$run_id))) {
        stop(
            "run_id values may contain only letters, numbers, periods, ",
            "underscores and hyphens.",
            call. = FALSE
        )
    }

    if (anyDuplicated(out$run_id)) {
        duplicates <- out |>
            dplyr::count(run_id, name = "n_records") |>
            dplyr::filter(n_records > 1L)

        print(duplicates)
        stop(
            "run_id values are not unique. Supply explicit unique run_id ",
            "values if sanitised method names collide.",
            call. = FALSE
        )
    }

    if (!("process_id" %in% names(out))) {
        out <- out |>
            dplyr::mutate(process_id = dplyr::row_number())
    }

    if (anyNA(out$process_id)) {
        stop(
            "process_id values must not be missing.",
            call. = FALSE
        )
    }

    if (!("run_order" %in% names(out))) {
        out <- out |>
            dplyr::group_by(process_id) |>
            dplyr::mutate(run_order = dplyr::row_number()) |>
            dplyr::ungroup()
    }

    out
}


## Allocate full-study runs to checkpointed Slurm processes.
##
## This is the checkpointed counterpart of
## find_iridis_setup_from_feasibility().  By default, runs from different
## method families are never packed into the same process.  Within a process,
## the run with the largest predicted elapsed time is attempted first.
find_iridis_setup_from_feasibility_checkpointed <- function(
    setup,
    simstudy_feasibility,
    time_each = 240,
    buffer = 1.5,
    back = 0,
    pack_by = "method",
    run_keys = c("case", "subcase_id", "method", "seed")
) {
    if (
        length(time_each) != 1L ||
        !is.finite(time_each) ||
        time_each <= 0
    ) {
        stop("time_each must be one positive number.", call. = FALSE)
    }

    if (
        length(buffer) != 1L ||
        !is.finite(buffer) ||
        buffer <= 0
    ) {
        stop("buffer must be one positive number.", call. = FALSE)
    }

    .checkpoint_check_cols(
        simstudy_feasibility,
        c(
            "case", "subcase_id", "method", "time",
            "sim_status", "fit_status", "pred_status", "gp_status"
        ),
        "simstudy_feasibility"
    )

    set.seed(1)

    timing_tab <- simstudy_feasibility |>
        dplyr::filter(
            sim_status == "ok",
            fit_status == "ok",
            pred_status == "ok",
            gp_status == "ok",
            is.finite(time)
        ) |>
        dplyr::group_by(case, subcase_id, method) |>
        dplyr::summarise(
            n_successful_trials = dplyr::n(),
            median_time = stats::median(time),
            max_time = max(time),
            .groups = "drop"
        ) |>
        dplyr::mutate(allowed_time = buffer * max_time)

    setup_with_times <- setup |>
        dplyr::left_join(
            timing_tab,
            by = c("case", "subcase_id", "method")
        )

    missing_times <- setup_with_times |>
        dplyr::filter(!is.finite(allowed_time))

    if (nrow(missing_times) > 0L) {
        print(
            missing_times |>
                dplyr::distinct(case, subcase_id, method)
        )
        stop(
            "No successful feasibility timing is available for ",
            nrow(missing_times),
            " full-simulation rows.",
            call. = FALSE
        )
    }

    if (max(setup_with_times$allowed_time) > time_each) {
        stop(
            "time_each must be at least ",
            ceiling(max(setup_with_times$allowed_time)),
            " minutes.",
            call. = FALSE
        )
    }

    if (is.null(pack_by) || length(pack_by) == 0L) {
        setup_with_times$.checkpoint_pack_group <- "all"
        pack_by_internal <- ".checkpoint_pack_group"
    } else {
        .checkpoint_check_cols(setup_with_times, pack_by, "setup")
        pack_by_internal <- pack_by
    }

    groups <- setup_with_times |>
        dplyr::group_by(
            dplyr::across(dplyr::all_of(pack_by_internal))
        ) |>
        dplyr::group_split(.keep = TRUE)

    process_offset <- 0L
    allocated <- vector("list", length(groups))

    for (ii in seq_along(groups)) {
        group_i <- groups[[ii]]

        local_process_id <- find_process_ids(
            allowed_time = group_i$allowed_time,
            time_each = time_each,
            back = back
        )

        group_i$process_id <- local_process_id + process_offset
        process_offset <- max(group_i$process_id)
        allocated[[ii]] <- group_i
    }

    out <- dplyr::bind_rows(allocated)

    if (".checkpoint_pack_group" %in% names(out)) {
        out <- dplyr::select(out, -.checkpoint_pack_group)
    }

    out <- out |>
        dplyr::group_by(process_id) |>
        dplyr::arrange(
            dplyr::desc(allowed_time),
            .by_group = TRUE
        ) |>
        dplyr::mutate(run_order = dplyr::row_number()) |>
        dplyr::ungroup()

    add_checkpoint_run_ids(out, run_keys = run_keys)
}


.checkpoint_format_slurm_time <- function(minutes) {
    if (
        length(minutes) != 1L ||
        !is.finite(minutes) ||
        minutes <= 0
    ) {
        stop("minutes must be one positive number.", call. = FALSE)
    }

    total_seconds <- ceiling(minutes * 60)
    days <- total_seconds %/% (24 * 60 * 60)
    remainder <- total_seconds %% (24 * 60 * 60)
    hours <- remainder %/% (60 * 60)
    remainder <- remainder %% (60 * 60)
    mins <- remainder %/% 60
    secs <- remainder %% 60

    if (days > 0L) {
        sprintf("%d-%02d:%02d:%02d", days, hours, mins, secs)
    } else {
        sprintf("%02d:%02d:%02d", hours, mins, secs)
    }
}


write_shell_script_iridis_checkpointed <- function(
    runs,
    time_each,
    mem,
    iridis_dir,
    script_name = "simstudy_run.slurm",
    r_script_name = "simstudy_script_checkpointed.R"
) {
    runs <- sort(unique(as.integer(runs)))

    if (length(runs) == 0L || anyNA(runs)) {
        stop("runs must contain at least one process ID.", call. = FALSE)
    }

    if (
        length(mem) != 1L ||
        !is.finite(mem) ||
        mem <= 0
    ) {
        stop("mem must be one positive number of GB.", call. = FALSE)
    }

    filename <- file.path(iridis_dir, script_name)
    time_text <- .checkpoint_format_slurm_time(time_each)

    lines <- c(
        "#!/bin/sh",
        "",
        "#SBATCH --nodes=1",
        paste0("#SBATCH --time=", time_text),
        paste0("#SBATCH --mem=", mem, "G"),
        "#SBATCH --error=storage/run-%A-%a.err",
        "#SBATCH --output=storage/run-%A-%a.out",
        paste0("#SBATCH --array=", paste(runs, collapse = ",")),
        "",
        "module load conda/python3",
        "conda activate tidyverse",
        "",
        paste("Rscript", r_script_name)
    )

    writeLines(lines, con = filename)
    Sys.chmod(filename, mode = "0755")

    invisible(filename)
}


write_R_script_iridis_checkpointed <- function(
    simulation_fun_name,
    iridis_dir,
    script_name = "simstudy_script_checkpointed.R"
) {
    filename <- file.path(iridis_dir, script_name)

    script <- paste0(
        "library(tidyverse)\n",
        "library(refund)\n",
        "devtools::load_all()\n\n",
        "sessionInfo()\n\n",
        "Sys.setenv(TZ = 'Europe/London')\n\n",
        "id <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))\n",
        "if (!is.finite(id)) stop('SLURM_ARRAY_TASK_ID is missing.')\n\n",
        "load('setup.Rda')\n",
        "load('subcases.Rda')\n\n",
        "simulation_fun <- get('",
        simulation_fun_name,
        "', mode = 'function')\n\n",
        "run_checkpointed_process_iridis(\n",
        "    id = id,\n",
        "    setup = setup,\n",
        "    subcases = subcases,\n",
        "    simulation_fun = simulation_fun,\n",
        "    iridis_dir = '.'\n",
        ")\n"
    )

    writeLines(script, con = filename)
    invisible(filename)
}


## Execute one checkpointed Slurm process.
##
## Every completed run is saved atomically before the next run begins.  On a
## repeated attempt, existing valid run files are skipped.
run_checkpointed_process_iridis <- function(
    id,
    setup,
    subcases,
    simulation_fun,
    iridis_dir = ".",
    simulation_arg_names = c("case", "subcase_id", "method", "seed"),
    gc_after_run = TRUE
) {
    .checkpoint_check_cols(
        setup,
        c("run_id", "process_id", simulation_arg_names)
    )

    setup_local <- setup |>
        dplyr::filter(process_id == id)

    if (nrow(setup_local) == 0L) {
        stop("No setup rows have process_id = ", id, ".", call. = FALSE)
    }

    if ("run_order" %in% names(setup_local)) {
        setup_local <- setup_local |>
            dplyr::arrange(run_order)
    }

    dir.create(
        file.path(iridis_dir, "storage", "runs"),
        recursive = TRUE,
        showWarnings = FALSE
    )
    dir.create(
        file.path(iridis_dir, "storage", "progress"),
        recursive = TRUE,
        showWarnings = FALSE
    )

    .checkpoint_append_event(
        iridis_dir,
        event = "PROCESS_START",
        process_id = id,
        message = paste0(nrow(setup_local), " assigned runs")
    )

    for (jj in seq_len(nrow(setup_local))) {
        setup_row <- setup_local[jj, , drop = FALSE]
        run_id <- setup_row$run_id[[1]]
        output_file <- .checkpoint_run_path(iridis_dir, run_id)

        if (.checkpoint_output_is_valid(output_file, run_id)) {
            .checkpoint_append_event(
                iridis_dir,
                event = "SKIP",
                process_id = id,
                setup_row = setup_row,
                message = "Valid checkpoint already exists"
            )
            next
        }

        if (file.exists(output_file)) {
            corrupt_file <- paste0(
                output_file,
                ".corrupt-",
                format(Sys.time(), "%Y%m%d-%H%M%S")
            )
            file.rename(output_file, corrupt_file)

            .checkpoint_append_event(
                iridis_dir,
                event = "CORRUPT",
                process_id = id,
                setup_row = setup_row,
                message = basename(corrupt_file)
            )
        }

        .checkpoint_append_event(
            iridis_dir,
            event = "START",
            process_id = id,
            setup_row = setup_row
        )

        run_args <- as.list(
            setup_row[1, simulation_arg_names, drop = FALSE]
        )
        run_args$subcases <- subcases

        worker_error <- NULL

        result <- tryCatch(
            do.call(simulation_fun, run_args),
            error = function(e) {
                worker_error <<- conditionMessage(e)
                .checkpoint_worker_error_result(
                    setup_row,
                    worker_error
                )
            }
        )

        if (!is.data.frame(result) || nrow(result) != 1L) {
            worker_error <- paste0(
                "simulation_fun returned ",
                if (is.data.frame(result)) nrow(result) else "a non-data-frame",
                "; exactly one row was expected."
            )
            result <- .checkpoint_worker_error_result(
                setup_row,
                worker_error
            )
        }

        result <- tibble::as_tibble(result) |>
            dplyr::mutate(
                run_id = as.character(run_id),
                process_id = as.integer(id),
                checkpoint_worker_status = if (
                    is.null(worker_error)
                ) "ok" else "error",
                checkpoint_worker_error = if (
                    is.null(worker_error)
                ) NA_character_ else worker_error,
                checkpoint_slurm_job_id = Sys.getenv(
                    "SLURM_JOB_ID",
                    unset = NA_character_
                ),
                checkpoint_slurm_array_job_id = Sys.getenv(
                    "SLURM_ARRAY_JOB_ID",
                    unset = NA_character_
                ),
                .before = 1
            )

        temporary_file <- paste0(
            output_file,
            ".tmp-",
            Sys.getpid()
        )

        saveRDS(result, temporary_file)

        if (!file.rename(temporary_file, output_file)) {
            unlink(temporary_file)
            stop(
                "Could not atomically move completed result to ",
                output_file,
                call. = FALSE
            )
        }

        .checkpoint_append_event(
            iridis_dir,
            event = if (is.null(worker_error)) "DONE" else "ERROR",
            process_id = id,
            setup_row = setup_row,
            message = if (is.null(worker_error)) "" else worker_error
        )

        rm(result)
        if (isTRUE(gc_after_run)) {
            invisible(gc())
        }
    }

    .checkpoint_append_event(
        iridis_dir,
        event = "PROCESS_DONE",
        process_id = id
    )

    invisible(TRUE)
}


## Create a new checkpointed Iridis run directory.
write_simstudy_iridis_checkpointed <- function(
    setup,
    subcases,
    simulation_fun,
    time_each,
    mem = 8,
    iridis_dir,
    simulation_fun_name = NULL
) {
    if (dir.exists(iridis_dir)) {
        stop(
            "Directory ",
            iridis_dir,
            " already exists. Please delete or rename it.",
            call. = FALSE
        )
    }

    setup <- add_checkpoint_run_ids(setup)

    dir.create(iridis_dir, recursive = TRUE)
    dir.create(
        file.path(iridis_dir, "storage", "runs"),
        recursive = TRUE
    )
    dir.create(
        file.path(iridis_dir, "storage", "progress"),
        recursive = TRUE
    )

    if (is.null(simulation_fun_name)) {
        simulation_fun_name <- deparse(substitute(simulation_fun))
    }

    if (
        length(simulation_fun_name) != 1L ||
        !grepl("^[A-Za-z.][A-Za-z0-9._]*$", simulation_fun_name)
    ) {
        stop(
            "simulation_fun_name must be the name of one function available ",
            "after devtools::load_all().",
            call. = FALSE
        )
    }

    save(setup, file = file.path(iridis_dir, "setup.Rda"))
    save(subcases, file = file.path(iridis_dir, "subcases.Rda"))

    write_R_script_iridis_checkpointed(
        simulation_fun_name,
        iridis_dir
    )

    process_ids <- sort(unique(setup$process_id))

    write_shell_script_iridis_checkpointed(
        runs = process_ids,
        time_each = time_each,
        mem = mem,
        iridis_dir = iridis_dir,
        script_name = "simstudy_run.slurm"
    )

    invisible(setup)
}


## Write an update Slurm script for specified checkpointed processes.
##
## Unlike write_update_script_iridis(), this function does not infer failed
## processes from missing process-level output files.  The caller supplies the
## process IDs, usually from find_timeout_ids_iridis().  Completed run-level
## checkpoints are skipped automatically when those processes are rerun.
write_update_script_iridis_checkpointed <- function(
    process_ids,
    time_each,
    mem = 8,
    iridis_dir,
    script_name = "simstudy_update.slurm"
) {
    load(file.path(iridis_dir, "setup.Rda"))

    process_ids <- sort(unique(as.integer(process_ids)))
    valid_process_ids <- sort(unique(setup$process_id))
    invalid <- setdiff(process_ids, valid_process_ids)

    if (length(process_ids) == 0L || anyNA(process_ids)) {
        stop(
            "process_ids must contain at least one valid process ID.",
            call. = FALSE
        )
    }

    if (length(invalid) > 0L) {
        stop(
            "Unknown process IDs: ",
            paste(invalid, collapse = ", "),
            call. = FALSE
        )
    }

    write_shell_script_iridis_checkpointed(
        runs = process_ids,
        time_each = time_each,
        mem = mem,
        iridis_dir = iridis_dir,
        script_name = script_name
    )

    invisible(process_ids)
}


.checkpoint_run_file_table <- function(iridis_dir) {
    runs_dir <- file.path(iridis_dir, "storage", "runs")

    if (!dir.exists(runs_dir)) {
        return(
            tibble::tibble(
                run_id = character(),
                file = character(),
                mtime = as.POSIXct(character())
            )
        )
    }

    files <- list.files(
        runs_dir,
        pattern = "^run-.*\\.rds$",
        full.names = TRUE
    )

    if (length(files) == 0L) {
        return(
            tibble::tibble(
                run_id = character(),
                file = character(),
                mtime = as.POSIXct(character())
            )
        )
    }

    basenames <- basename(files)
    run_id <- sub("^run-(.*)\\.rds$", "\\1", basenames)

    tibble::tibble(
        run_id = run_id,
        file = files,
        mtime = file.info(files)$mtime
    )
}


find_completed_run_ids_iridis_checkpointed <- function(
    iridis_dir,
    validate = TRUE
) {
    file_tab <- .checkpoint_run_file_table(iridis_dir)

    if (nrow(file_tab) == 0L) {
        return(character())
    }

    if (!isTRUE(validate)) {
        return(file_tab$run_id)
    }

    valid <- purrr::map2_lgl(
        file_tab$file,
        file_tab$run_id,
        .checkpoint_output_is_valid
    )

    file_tab$run_id[valid]
}


.checkpoint_empty_run_status <- function() {
    tibble::tibble(
        run_id = character(),
        process_id = integer(),
        case = character(),
        subcase_id = integer(),
        method = character(),
        seed = integer(),
        time = numeric(),
        time_sim = numeric(),
        time_fit = numeric(),
        time_pred = numeric(),
        time_gp = numeric(),
        sim_status = character(),
        fit_status = character(),
        pred_status = character(),
        gp_status = character(),
        k_status = character(),
        sim_error = character(),
        fit_error = character(),
        pred_error = character(),
        gp_error = character(),
        k_error = character(),
        checkpoint_worker_status = character(),
        checkpoint_worker_error = character(),
        has_CI = logical(),
        checkpoint_source_dir = character(),
        checkpoint_file = character(),
        checkpoint_mtime = as.POSIXct(character()),
        checkpoint_source_rank = integer()
    )
}


## Return the checkpoint files to use, resolving duplicates by permanent run_id.
##
## This is deliberately a file index, not the full simulation output.  It can be
## saved after finalisation and used later by analysis code which reads and
## summarises one checkpoint at a time.
checkpointed_result_file_index <- function(
    iridis_dirs,
    prefer = c("last", "latest_mtime"),
    warn_duplicates = TRUE
) {
    prefer <- match.arg(prefer)

    if (length(iridis_dirs) == 0L) {
        return(
            tibble::tibble(
                run_id = character(),
                file = character(),
                checkpoint_mtime = as.POSIXct(character()),
                checkpoint_source_dir = character(),
                checkpoint_source_rank = integer()
            )
        )
    }

    components <- purrr::map2(
        iridis_dirs,
        seq_along(iridis_dirs),
        function(iridis_dir, source_rank) {
            .checkpoint_run_file_table(iridis_dir) |>
                dplyr::transmute(
                    run_id = as.character(run_id),
                    file = file,
                    checkpoint_mtime = mtime,
                    checkpoint_source_dir = normalizePath(
                        iridis_dir,
                        mustWork = FALSE
                    ),
                    checkpoint_source_rank = as.integer(source_rank)
                )
        }
    )

    all_files <- dplyr::bind_rows(components)

    if (nrow(all_files) == 0L) {
        return(all_files)
    }

    duplicate_ids <- all_files |>
        dplyr::count(run_id, name = "n_records") |>
        dplyr::filter(n_records > 1L)

    if (isTRUE(warn_duplicates) && nrow(duplicate_ids) > 0L) {
        warning(
            nrow(duplicate_ids),
            " run IDs occur in more than one checkpoint directory; ",
            "the preferred checkpoint file will be retained.",
            call. = FALSE
        )
    }

    if (prefer == "last") {
        all_files <- all_files |>
            dplyr::arrange(
                run_id,
                checkpoint_source_rank,
                checkpoint_mtime
            )
    } else {
        all_files <- all_files |>
            dplyr::arrange(
                run_id,
                checkpoint_mtime,
                checkpoint_source_rank
            )
    }

    all_files |>
        dplyr::group_by(run_id) |>
        dplyr::slice_tail(n = 1L) |>
        dplyr::ungroup()
}


.checkpoint_read_run_status_file <- function(
    file,
    expected_run_id,
    checkpoint_source_dir = NA_character_,
    checkpoint_mtime = as.POSIXct(NA),
    checkpoint_source_rank = NA_integer_
) {
    result <- tryCatch(
        readRDS(file),
        error = function(e) {
            warning(
                "Could not read ",
                file,
                ": ",
                conditionMessage(e),
                call. = FALSE
            )
            NULL
        }
    )

    if (is.null(result)) {
        return(NULL)
    }

    if (!is.data.frame(result) || nrow(result) != 1L) {
        warning(
            file,
            " does not contain exactly one result row.",
            call. = FALSE
        )
        return(NULL)
    }

    result <- tibble::as_tibble(result)

    if (!("run_id" %in% names(result))) {
        result$run_id <- expected_run_id
    }

    if (!identical(
        as.character(result$run_id[[1]]),
        as.character(expected_run_id)
    )) {
        warning(
            "run_id inside ",
            file,
            " does not match its filename.",
            call. = FALSE
        )
        return(NULL)
    }

    has_CI <- FALSE

    if ("has_CI" %in% names(result)) {
        has_CI <- isTRUE(result$has_CI[[1]])
    } else if ("pred_data" %in% names(result)) {
        has_CI <- has_CI_predictions(result$pred_data[[1]])
    }

    needed_cols <- c(
        "run_id", "process_id", "case", "subcase_id", "method", "seed",
        "time", "time_sim", "time_fit", "time_pred", "time_gp",
        "sim_status", "fit_status", "pred_status", "gp_status", "k_status",
        "sim_error", "fit_error", "pred_error", "gp_error", "k_error",
        "checkpoint_worker_status", "checkpoint_worker_error"
    )

    result2 <- add_missing_cols(result, needed_cols)

    result2 |>
        dplyr::mutate(
            run_id = as.character(run_id),
            has_CI = has_CI,
            checkpoint_source_dir = checkpoint_source_dir,
            checkpoint_file = file,
            checkpoint_mtime = checkpoint_mtime,
            checkpoint_source_rank = as.integer(checkpoint_source_rank)
        ) |>
        dplyr::select(
            run_id,
            process_id,
            case,
            subcase_id,
            method,
            seed,
            time,
            time_sim,
            time_fit,
            time_pred,
            time_gp,
            sim_status,
            fit_status,
            pred_status,
            gp_status,
            k_status,
            sim_error,
            fit_error,
            pred_error,
            gp_error,
            k_error,
            checkpoint_worker_status,
            checkpoint_worker_error,
            has_CI,
            checkpoint_source_dir,
            checkpoint_file,
            checkpoint_mtime,
            checkpoint_source_rank
        )
}


## Read only scalar status fields from checkpoint files.
##
## Each checkpoint is still read from disk, but only one at a time and the large
## list columns such as pred_data, GP0 and GP_hat are discarded immediately.
## This is intended for finalisation/status construction, not metric analysis.
read_checkpointed_run_status <- function(
    checkpoint_index,
    verbose = FALSE
) {
    if (nrow(checkpoint_index) == 0L) {
        return(.checkpoint_empty_run_status())
    }

    if (isTRUE(verbose)) {
        message(
            "Reading scalar status fields from ",
            nrow(checkpoint_index),
            " checkpoint files"
        )
    }

    components <- purrr::pmap(
        list(
            checkpoint_index$file,
            checkpoint_index$run_id,
            checkpoint_index$checkpoint_source_dir,
            checkpoint_index$checkpoint_mtime,
            checkpoint_index$checkpoint_source_rank
        ),
        .checkpoint_read_run_status_file
    )

    output <- dplyr::bind_rows(components)

    if (nrow(output) == 0L) {
        return(.checkpoint_empty_run_status())
    }

    output
}


combine_checkpointed_run_status <- function(
    iridis_dirs,
    prefer = c("last", "latest_mtime"),
    warn_duplicates = TRUE,
    verbose = FALSE
) {
    checkpoint_index <- checkpointed_result_file_index(
        iridis_dirs = iridis_dirs,
        prefer = prefer,
        warn_duplicates = warn_duplicates
    )

    read_checkpointed_run_status(
        checkpoint_index = checkpoint_index,
        verbose = verbose
    )
}


read_simstudy_iridis_checkpointed <- function(
    iridis_dir,
    run_ids = NULL,
    drop = NULL,
    condition = NULL,
    verbose = FALSE
) {
    condition_quo <- rlang::enquo(condition)
    file_tab <- .checkpoint_run_file_table(iridis_dir)

    if (!is.null(run_ids)) {
        file_tab <- file_tab |>
            dplyr::filter(run_id %in% as.character(run_ids))
    }

    if (isTRUE(verbose)) {
        message("Reading ", nrow(file_tab), " checkpointed run files")
    }

    components <- purrr::map2(
        file_tab$file,
        file_tab$run_id,
        function(file, expected_run_id) {
            result <- tryCatch(
                readRDS(file),
                error = function(e) {
                    warning(
                        "Could not read ",
                        file,
                        ": ",
                        conditionMessage(e),
                        call. = FALSE
                    )
                    NULL
                }
            )

            if (is.null(result)) {
                return(NULL)
            }

            if (!is.data.frame(result) || nrow(result) != 1L) {
                warning(
                    file,
                    " does not contain exactly one result row.",
                    call. = FALSE
                )
                return(NULL)
            }

            result <- tibble::as_tibble(result)

            if (!("run_id" %in% names(result))) {
                result$run_id <- expected_run_id
            }

            if (!identical(
                as.character(result$run_id[[1]]),
                as.character(expected_run_id)
            )) {
                warning(
                    "run_id inside ",
                    file,
                    " does not match its filename.",
                    call. = FALSE
                )
                return(NULL)
            }

            result$checkpoint_source_dir <- normalizePath(
                iridis_dir,
                mustWork = FALSE
            )
            result$checkpoint_file <- file
            result$checkpoint_mtime <- file.info(file)$mtime

            result
        }
    )

    output <- dplyr::bind_rows(components)

    if (!rlang::quo_is_null(condition_quo) && nrow(output) > 0L) {
        output <- dplyr::filter(output, !!condition_quo)
    }

    if (!is.null(drop) && nrow(output) > 0L) {
        output <- dplyr::select(
            output,
            -dplyr::any_of(drop)
        )
    }

    output
}


## Combine primary and recovery directories by permanent run_id.
##
## With prefer = "last", later directories in iridis_dirs take precedence.
## With prefer = "latest_mtime", the newest checkpoint takes precedence.
combine_checkpointed_iridis_results <- function(
    iridis_dirs,
    prefer = c("last", "latest_mtime"),
    warn_duplicates = TRUE
) {
    prefer <- match.arg(prefer)

    if (length(iridis_dirs) == 0L) {
        return(tibble::tibble())
    }

    components <- purrr::map2(
        iridis_dirs,
        seq_along(iridis_dirs),
        function(iridis_dir, source_rank) {
            read_simstudy_iridis_checkpointed(iridis_dir) |>
                dplyr::mutate(
                    checkpoint_source_rank = as.integer(source_rank)
                )
        }
    )

    all_results <- dplyr::bind_rows(components)

    if (nrow(all_results) == 0L) {
        return(all_results)
    }

    duplicate_ids <- all_results |>
        dplyr::count(run_id, name = "n_records") |>
        dplyr::filter(n_records > 1L)

    if (isTRUE(warn_duplicates) && nrow(duplicate_ids) > 0L) {
        warning(
            nrow(duplicate_ids),
            " run IDs occur in more than one checkpoint directory; ",
            "the preferred result will be retained.",
            call. = FALSE
        )
    }

    if (prefer == "last") {
        all_results <- all_results |>
            dplyr::arrange(
                run_id,
                checkpoint_source_rank,
                checkpoint_mtime
            )
    } else {
        all_results <- all_results |>
            dplyr::arrange(
                run_id,
                checkpoint_mtime,
                checkpoint_source_rank
            )
    }

    all_results |>
        dplyr::group_by(run_id) |>
        dplyr::slice_tail(n = 1L) |>
        dplyr::ungroup()
}


find_missing_runs_iridis_checkpointed <- function(
    iridis_dir,
    setup = NULL,
    validate = TRUE
) {
    if (is.null(setup)) {
        load(file.path(iridis_dir, "setup.Rda"))
    }

    .checkpoint_check_cols(setup, "run_id")

    completed_ids <- find_completed_run_ids_iridis_checkpointed(
        iridis_dir,
        validate = validate
    )

    setup |>
        dplyr::filter(!as.character(run_id) %in% completed_ids)
}


find_processes_with_missing_runs_iridis_checkpointed <- function(
    iridis_dir,
    setup = NULL,
    validate = TRUE
) {
    missing <- find_missing_runs_iridis_checkpointed(
        iridis_dir,
        setup = setup,
        validate = validate
    )

    if (nrow(missing) == 0L) {
        return(
            tibble::tibble(
                process_id = integer(),
                n_missing_runs = integer()
            )
        )
    }

    missing |>
        dplyr::count(process_id, name = "n_missing_runs") |>
        dplyr::arrange(process_id)
}

read_iridis_progress_checkpointed <- function(
    iridis_dir,
    process_ids = NULL
) {
    progress_dir <- file.path(
        iridis_dir,
        "storage",
        "progress"
    )

    empty_progress <- tibble::tibble(
        process_id = integer(),
        timestamp = as.POSIXct(
            character(),
            tz = "UTC"
        ),
        event = character(),
        run_id = integer(),
        case = character(),
        subcase_id = integer(),
        method = character(),
        seed = integer(),
        message = character()
    )

    if (!dir.exists(progress_dir)) {
        return(empty_progress)
    }

    progress_files <- list.files(
        progress_dir,
        pattern = "^progress-[0-9]+[.]log$",
        full.names = TRUE
    )

    if (length(progress_files) == 0L) {
        return(empty_progress)
    }

    file_process_ids <- as.integer(
        sub(
            "^progress-([0-9]+)[.]log$",
            "\\1",
            basename(progress_files)
        )
    )

    if (!is.null(process_ids)) {
        keep <- file_process_ids %in% as.integer(process_ids)

        progress_files <- progress_files[keep]
        file_process_ids <- file_process_ids[keep]
    }

    if (length(progress_files) == 0L) {
        return(empty_progress)
    }

    components <- purrr::map2(
        progress_files,
        file_process_ids,
        function(file, process_id) {
            lines <- readLines(
                file,
                warn = FALSE
            )

            if (length(lines) == 0L) {
                return(empty_progress)
            }

            parsed <- purrr::map_dfr(
                lines,
                function(line) {
                    ## Expected examples:
                    ## 2026-... START run_id=12 case=2re subcase_id=1 method=AdaStruMM seed=4
                    ## 2026-... DONE run_id=12
                    ## 2026-... SKIP run_id=12
                    ## 2026-... ERROR run_id=12 message=...
                    parts <- strsplit(
                        line,
                        "\\s+"
                    )[[1]]

                    if (length(parts) < 3L) {
                        return(
                            tibble::tibble(
                                timestamp = as.POSIXct(
                                    NA_character_,
                                    tz = "UTC"
                                ),
                                event = NA_character_,
                                run_id = NA_integer_,
                                case = NA_character_,
                                subcase_id = NA_integer_,
                                method = NA_character_,
                                seed = NA_integer_,
                                message = line
                            )
                        )
                    }

                    timestamp_text <- paste(
                        parts[1L],
                        parts[2L]
                    )

                    event <- parts[3L]

                    key_value_parts <- parts[-seq_len(3L)]

                    key_values <- stats::setNames(
                        rep(
                            NA_character_,
                            length(key_value_parts)
                        ),
                        rep(
                            NA_character_,
                            length(key_value_parts)
                        )
                    )

                    if (length(key_value_parts) > 0L) {
                        keys <- sub(
                            "=.*$",
                            "",
                            key_value_parts
                        )

                        values <- sub(
                            "^[^=]*=",
                            "",
                            key_value_parts
                        )

                        key_values <- stats::setNames(
                            values,
                            keys
                        )
                    }

                    get_value <- function(name) {
                        if (name %in% names(key_values)) {
                            key_values[[name]]
                        } else {
                            NA_character_
                        }
                    }

                    tibble::tibble(
                        timestamp = as.POSIXct(
                            timestamp_text,
                            tz = "UTC"
                        ),
                        event = event,
                        run_id = suppressWarnings(
                            as.integer(
                                get_value("run_id")
                            )
                        ),
                        case = get_value("case"),
                        subcase_id = suppressWarnings(
                            as.integer(
                                get_value("subcase_id")
                            )
                        ),
                        method = get_value("method"),
                        seed = suppressWarnings(
                            as.integer(
                                get_value("seed")
                            )
                        ),
                        message = line
                    )
                }
            )

            parsed %>%
                dplyr::mutate(
                    process_id = process_id,
                    .before = 1
                )
        }
    )

    out <- dplyr::bind_rows(components)

    if (nrow(out) == 0L || !"process_id" %in% names(out)) {
        return(empty_progress)
    }

    out %>%
        dplyr::arrange(
            process_id,
            timestamp
        )
}


## Identify runs which were started but have no later DONE, ERROR or SKIP
## event and do not have a valid checkpoint.  These are the best candidates
## for the run active when a packed process failed.
find_active_runs_at_failure_iridis_checkpointed <- function(
    iridis_dir,
    process_ids = NULL
) {
    progress <- read_iridis_progress_checkpointed(
        iridis_dir,
        process_ids = process_ids
    )

    empty_active <- tibble::tibble(
        process_id = integer(),
        run_id = integer(),
        case = character(),
        subcase_id = integer(),
        method = character(),
        seed = integer(),
        started_at = as.POSIXct(
            character(),
            tz = "UTC"
        )
    )

    if (nrow(progress) == 0L) {
        return(empty_active)
    }

    starts <- progress %>%
        dplyr::filter(
            event == "START",
            !is.na(run_id)
        ) %>%
        dplyr::select(
            process_id,
            run_id,
            case,
            subcase_id,
            method,
            seed,
            started_at = timestamp
        )

    finishes <- progress %>%
        dplyr::filter(
            event %in% c(
                "DONE",
                "ERROR",
                "SKIP"
            ),
            !is.na(run_id)
        ) %>%
        dplyr::select(
            process_id,
            run_id
        ) %>%
        dplyr::distinct()

    starts %>%
        dplyr::anti_join(
            finishes,
            by = c(
                "process_id",
                "run_id"
            )
        ) %>%
        dplyr::group_by(
            process_id
        ) %>%
        dplyr::slice_max(
            started_at,
            n = 1,
            with_ties = FALSE
        ) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(
            process_id
        )
}


## Construct a recovery setup containing only incomplete runs.
##
## If process_ids is supplied, incomplete runs originally assigned to those
## processes are selected.  If run_ids is supplied, only those permanent runs
## are selected.  With isolate = TRUE, every remaining run receives its own
## new process_id, allowing definitive attribution of OOMs and timeouts.
make_isolated_recovery_setup <- function(
    setup,
    completed_run_ids = character(),
    process_ids = NULL,
    run_ids = NULL,
    isolate = TRUE
) {
    .checkpoint_check_cols(setup, c("run_id", "process_id"))

    selected <- setup

    if (!is.null(process_ids)) {
        selected <- selected |>
            dplyr::filter(process_id %in% process_ids)
    }

    if (!is.null(run_ids)) {
        selected <- selected |>
            dplyr::filter(as.character(run_id) %in% as.character(run_ids))
    }

    selected <- selected |>
        dplyr::filter(
            !as.character(run_id) %in% as.character(completed_run_ids)
        ) |>
        dplyr::mutate(
            original_process_id = process_id,
            original_run_order = if (
                "run_order" %in% names(selected)
            ) run_order else NA_integer_
        )

    if (isTRUE(isolate) && nrow(selected) > 0L) {
        selected <- selected |>
            dplyr::arrange(
                original_process_id,
                original_run_order,
                run_id
            ) |>
            dplyr::mutate(
                process_id = dplyr::row_number(),
                run_order = 1L
            )
    }

    selected
}


## Create an isolated checkpointed recovery directory from a primary run.
write_isolated_recovery_iridis_checkpointed <- function(
    source_iridis_dir,
    recovery_iridis_dir,
    simulation_fun,
    time_each,
    mem = 8,
    process_ids = NULL,
    run_ids = NULL
) {
    simulation_fun_name <- deparse(substitute(simulation_fun))

    load(file.path(source_iridis_dir, "setup.Rda"))
    load(file.path(source_iridis_dir, "subcases.Rda"))

    completed_run_ids <- find_completed_run_ids_iridis_checkpointed(
        source_iridis_dir,
        validate = TRUE
    )

    recovery_setup <- make_isolated_recovery_setup(
        setup = setup,
        completed_run_ids = completed_run_ids,
        process_ids = process_ids,
        run_ids = run_ids,
        isolate = TRUE
    )

    if (nrow(recovery_setup) == 0L) {
        stop(
            "There are no incomplete runs matching the requested recovery ",
            "set.",
            call. = FALSE
        )
    }

    write_simstudy_iridis_checkpointed(
        setup = recovery_setup,
        subcases = subcases,
        simulation_fun = simulation_fun,
        time_each = time_each,
        mem = mem,
        iridis_dir = recovery_iridis_dir,
        simulation_fun_name = simulation_fun_name
    )

    invisible(recovery_setup)
}


.checkpoint_validate_isolated_setup <- function(iridis_dir) {
    load(file.path(iridis_dir, "setup.Rda"))

    counts <- setup |>
        dplyr::count(process_id, name = "n_runs")

    if (any(counts$n_runs != 1L)) {
        stop(
            iridis_dir,
            " is not an isolated recovery directory: at least one process ",
            "contains more than one run.",
            call. = FALSE
        )
    }

    setup
}


## Definitive OOM run IDs from a one-run-per-process recovery directory.
find_oom_run_ids_iridis_checkpointed <- function(
    iridis_dir,
    min_mem_gb = 8,
    mem_limit_if_missing_gb = NULL
) {
    .checkpoint_validate_isolated_setup(iridis_dir)

    find_oom_ids_iridis(
        iridis_dir,
        min_mem_gb = min_mem_gb,
        mem_limit_if_missing_gb = mem_limit_if_missing_gb
    ) |>
        dplyr::distinct(run_id) |>
        dplyr::pull(run_id) |>
        as.character()
}


## Definitive timeout run IDs from a one-run-per-process recovery directory.
find_timeout_run_ids_iridis_checkpointed <- function(
    iridis_dir,
    min_time_minutes,
    time_limit_if_missing_minutes = NULL
) {
    .checkpoint_validate_isolated_setup(iridis_dir)

    find_timeout_ids_iridis(
        iridis_dir,
        min_time_minutes = min_time_minutes,
        time_limit_if_missing_minutes = time_limit_if_missing_minutes
    ) |>
        dplyr::distinct(run_id) |>
        dplyr::pull(run_id) |>
        as.character()
}


## Run-level status construction for checkpointed primary and recovery runs.
##
## oom_run_ids and timeout_run_ids should contain only definitively attributed
## failures, normally obtained from isolated recovery directories.  Runs still
## awaiting recovery may be supplied in pending_run_ids.
make_complete_run_status_checkpointed <- function(
    setup,
    simstudy,
    oom_run_ids = NULL,
    timeout_run_ids = NULL,
    pending_run_ids = NULL
) {
    run_keys <- c("case", "subcase_id", "method", "seed")

    .checkpoint_check_cols(setup, c("run_id", run_keys))
    .checkpoint_check_cols(simstudy, "run_id")

    if (anyDuplicated(setup$run_id)) {
        stop("setup contains duplicated run_id values.", call. = FALSE)
    }

    if (anyDuplicated(simstudy$run_id)) {
        stop(
            "simstudy contains duplicated run_id values. Combine recovery ",
            "directories first.",
            call. = FALSE
        )
    }

    setup_keys <- setup |>
        dplyr::select(run_id, dplyr::all_of(run_keys))

    result_keys <- simstudy |>
        dplyr::select(
            run_id,
            dplyr::any_of(run_keys)
        )

    common_keys <- intersect(run_keys, names(result_keys))

    if (length(common_keys) > 0L && nrow(result_keys) > 0L) {
        key_check <- setup_keys |>
            dplyr::inner_join(
                result_keys,
                by = "run_id",
                suffix = c("_setup", "_result")
            )

        mismatch <- rep(FALSE, nrow(key_check))

        for (key in common_keys) {
            left <- key_check[[paste0(key, "_setup")]]
            right <- key_check[[paste0(key, "_result")]]
            mismatch <- mismatch |
                as.character(left) != as.character(right)
        }

        if (any(mismatch, na.rm = TRUE)) {
            print(key_check[mismatch, , drop = FALSE])
            stop(
                "Run keys in simstudy do not match setup for some run_id ",
                "values.",
                call. = FALSE
            )
        }
    }

    ## The finalisation stage should pass a slim run-status table containing
    ## has_CI, rather than a full simstudy table containing pred_data.  For
    ## backwards compatibility, compute has_CI here if pred_data is present.
    if (!("has_CI" %in% names(simstudy))) {
        if ("pred_data" %in% names(simstudy)) {
            simstudy$has_CI <- purrr::map_lgl(
                simstudy$pred_data,
                has_CI_predictions
            )
        } else {
            simstudy$has_CI <- FALSE
        }
    }

    needed_cols <- c(
        "run_id",
        "time", "time_sim", "time_fit", "time_pred", "time_gp",
        "sim_status", "fit_status", "pred_status", "gp_status", "k_status",
        "sim_error", "fit_error", "pred_error", "gp_error", "k_error",
        "checkpoint_worker_status", "checkpoint_worker_error",
        "process_id", "has_CI"
    )

    simstudy2 <- add_missing_cols(simstudy, needed_cols)

    run_status <- simstudy2 |>
        dplyr::select(
            dplyr::all_of(needed_cols)
        )

    setup |>
        dplyr::left_join(
            run_status,
            by = "run_id",
            suffix = c("_setup", "_result")
        ) |>
        dplyr::mutate(
            process_id = dplyr::coalesce(
                process_id_result,
                process_id_setup
            ),

            missing_from_simstudy =
                is.na(sim_status) &
                is.na(checkpoint_worker_status),

            worker_fail =
                !missing_from_simstudy &
                dplyr::coalesce(
                    checkpoint_worker_status == "error",
                    FALSE
                ),

            oom_fail =
                missing_from_simstudy &
                as.character(run_id) %in% as.character(oom_run_ids),

            timeout_fail =
                missing_from_simstudy &
                !oom_fail &
                as.character(run_id) %in% as.character(timeout_run_ids),

            pending_recovery =
                missing_from_simstudy &
                !oom_fail &
                !timeout_fail &
                as.character(run_id) %in% as.character(pending_run_ids),

            unknown_missing =
                missing_from_simstudy &
                !oom_fail &
                !timeout_fail,

            sim_fail =
                !missing_from_simstudy &
                (
                    worker_fail |
                    is.na(sim_status) |
                    sim_status != "ok"
                ),

            fit_fail =
                !missing_from_simstudy &
                !dplyr::coalesce(worker_fail, FALSE) &
                sim_status == "ok" &
                fit_status != "ok",

            pred_fail =
                !missing_from_simstudy &
                !dplyr::coalesce(worker_fail, FALSE) &
                sim_status == "ok" &
                fit_status == "ok" &
                pred_status != "ok",

            gp_fail =
                !missing_from_simstudy &
                !dplyr::coalesce(worker_fail, FALSE) &
                sim_status == "ok" &
                fit_status == "ok" &
                pred_status == "ok" &
                gp_status != "ok",

            k_fail =
                !missing_from_simstudy &
                !dplyr::coalesce(worker_fail, FALSE) &
                sim_status == "ok" &
                fit_status == "ok" &
                k_status != "ok",

            ci_fail =
                !missing_from_simstudy &
                !dplyr::coalesce(worker_fail, FALSE) &
                sim_status == "ok" &
                fit_status == "ok" &
                pred_status == "ok" &
                !dplyr::coalesce(has_CI, FALSE)
        ) |>
        dplyr::select(
            -dplyr::any_of(c("process_id_setup", "process_id_result"))
        )
}
