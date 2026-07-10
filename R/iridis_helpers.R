write_R_script_iridis <- function(simulation_fun_name, iridis_dir) {
    filename <- file.path(iridis_dir, "simstudy_script.R")
    cat(
"library(tidyverse)

library(refund)
devtools::load_all()

sessionInfo()

Sys.setenv(TZ = \"Europe/London\")

id <- as.integer(Sys.getenv(\"SLURM_ARRAY_TASK_ID\"))

load(\"setup.Rda\")
load(\"subcases.Rda\")

setup_local <- setup %>% filter(process_id == id)

simstudy_local <- pmap_dfr(setup_local, ", simulation_fun_name, ", subcases = subcases)

save(simstudy_local, file = file.path(\"storage\", paste(id, \".Rout\", sep = \"\")))",
file = filename,
sep = "")
}

write_shell_script_iridis <- function(runs, time_each, mem, iridis_dir, script_name) {
    filename <- file.path(iridis_dir, script_name)
    cat(
"#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=00:", time_each, ":00
#SBATCH --mem=", mem, "G
#SBATCH --error=storage/run-%A-%a.err
#SBATCH --output=storage/run-%A-%a.out
#SBATCH --array=", paste(runs, collapse = ","),"

module load conda/python3
conda activate tidyverse

Rscript simstudy_script.R $SLURM_ARRAY_TASK_ID",
file = filename,
sep = "")
    Sys.chmod(filename)
}


write_simstudy_iridis <- function(setup, subcases, simulation_fun, time_each, mem = 4,
                                  iridis_dir) {
    if(dir.exists(iridis_dir))
        stop("Directory ", iridis_dir, " already exists. Please delete or rename it.")
    else {
        dir.create(iridis_dir)
        storage_dir <- file.path(iridis_dir, "storage")
        dir.create(storage_dir)
    }
    
    simulation_fun_name <- deparse(substitute(simulation_fun))
    ## check if setup has process ids, if not add one per row
    if(!("process_id" %in% names(setup)))
        setup <- setup %>% mutate(process_id = 1:nrow(setup))
    save(setup, file = file.path(iridis_dir, "setup.Rda"))
    save(subcases, file = file.path(iridis_dir, "subcases.Rda"))
    n_processes <- get_n_processes(setup)
    write_R_script_iridis(simulation_fun_name, iridis_dir)
    write_shell_script_iridis(1:n_processes, time_each, mem, iridis_dir, "simstudy_run.slurm")
}

read_simstudy_local <- function(id, iridis_dir, drop, summary_fun,
                                condition, ...) {
    filename <- file.path(iridis_dir, "storage", paste(id, ".Rout", sep = ""))
    if(file.exists(filename)) {
        load(filename)
        result <- simstudy_local %>% select(-!!drop)
        if(!is.null(condition)) {
            result <- filter(result, eval(condition))
        }
        if(!is.null(summary_fun)) {
            result <- summary_fun(result, ...)
        }
    } else {
        warning("Case ", id, " is missing", call. = FALSE)
        result <- NULL
    }
    
    result
}

read_simstudy_iridis <- function(iridis_dir, ids = NULL, drop = NULL, summary_fun = NULL,
                                 condition = NULL, verbose = FALSE, ...) {
    load(file.path(iridis_dir, "setup.Rda"))
    if(!("process_id" %in% names(setup)))
        setup <- setup %>% mutate(process_id = 1:nrow(setup))

    if(!is.null(condition)) {
        setup <- setup %>% filter(eval(condition))
    }
    if(is.null(ids)) {
        ids <- unique(setup$process_id)
    }
    if(verbose) {
        cat("reading ids", ids, "\n")
    }
    components <- lapply(ids, read_simstudy_local,
                         iridis_dir = iridis_dir, drop = drop,
                         summary_fun = summary_fun,
                         condition = condition, ...)
    output <- bind_rows(components)
    
}



find_process_ids <- function(allowed_time, time_each, back) {
    if(any(is.na(allowed_time)))
        stop("do not have timing information for some setups")
    ## use large multiplier as knapsack::binpacking needs weights to be integers
    allowed_time_int <- ceiling(1000 * allowed_time)
    time_each_int <- 1000 * time_each
    time_ranks <- rank(-allowed_time_int, ties.method = "random")
    bins <- knapsack::binpacking(sort(allowed_time_int, decreasing = TRUE),
                                 cap = time_each_int, back = back)
    bins$xbins[time_ranks]                  
}


## back = -1 gives an exact solution to the bin packing problem
## increase to 0 or positive integer to give a quicker approximate solution
find_iridis_setup_from_feasibility <- function(
    setup,
    simstudy_feasibility,
    time_each = 240,
    buffer = 1.5,
    back = 0
) {
    set.seed(1)

    timing_tab <- simstudy_feasibility %>%
        dplyr::filter(
            sim_status == "ok",
            fit_status == "ok",
            pred_status == "ok",
            gp_status == "ok",
            is.finite(time)
        ) %>%
        dplyr::group_by(
            case,
            subcase_id,
            method
        ) %>%
        dplyr::summarise(
            n_successful_trials = dplyr::n(),
            median_time = stats::median(time),
            max_time = max(time),
            .groups = "drop"
        ) %>%
        dplyr::mutate(
            allowed_time = buffer * max_time
        )

    setup_with_times <- setup %>%
        dplyr::left_join(
            timing_tab,
            by = c(
                "case",
                "subcase_id",
                "method"
            )
        )

    missing_times <- setup_with_times %>%
        dplyr::filter(
            !is.finite(allowed_time)
        )

    if (nrow(missing_times) > 0) {
        stop(
            "No successful feasibility timing is available for ",
            nrow(missing_times),
            " full-simulation rows."
        )
    }

    min_time_each <- max(
        setup_with_times$allowed_time
    )

    if (min_time_each > time_each) {
        stop(
            "time_each must be at least ",
            ceiling(min_time_each),
            " minutes"
        )
    }

    setup_with_times %>%
        dplyr::mutate(
            process_id = find_process_ids(
                allowed_time,
                time_each,
                back
            )
        )
}


get_n_processes <- function(setup) {
    length(unique(setup$process_id))
}

is_id_missing <- function(id, iridis_dir) {
    filename <- file.path(iridis_dir, "storage", paste(id, ".Rout", sep = ""))
    !file.exists(filename)
}

find_missing_ids <- function(iridis_dir) {
    load(file.path(iridis_dir, "setup.Rda"))
    ids <- unique(setup$process_id)
    ids[sapply(ids, is_id_missing, iridis_dir = iridis_dir)]
}


parse_iridis_log_filename <- function(files) {
    base <- basename(files)

    m <- stringr::str_match(
        base,
        "^run-([0-9]+)-([0-9]+)\\.(out|err)$"
    )

    tibble::tibble(
        file = files,
        array_job_id = as.integer(m[, 2]),
        process_id = as.integer(m[, 3]),
        type = m[, 4]
    ) %>%
        dplyr::filter(!is.na(array_job_id), !is.na(process_id))
}


parse_memory_to_gb <- function(x) {
    if (length(x) == 0 || is.na(x)) {
        return(NA_real_)
    }

    m <- stringr::str_match(
        x,
        "([0-9.]+)\\s*(KB|MB|GB|TB)"
    )

    if (is.na(m[1, 1])) {
        return(NA_real_)
    }

    value <- as.numeric(m[1, 2])
    unit <- m[1, 3]

    dplyr::case_when(
        unit == "KB" ~ value / 1024^2,
        unit == "MB" ~ value / 1024,
        unit == "GB" ~ value,
        unit == "TB" ~ value * 1024,
        TRUE ~ NA_real_
    )
}

parse_slurm_duration_minutes <- function(x) {
    if (length(x) == 0 || is.na(x)) {
        return(NA_real_)
    }

    x <- trimws(x)

    ## Supports HH:MM:SS and D-HH:MM:SS.
    days <- 0

    if (grepl("-", x, fixed = TRUE)) {
        pieces <- strsplit(x, "-", fixed = TRUE)[[1]]

        if (length(pieces) != 2) {
            return(NA_real_)
        }

        days <- suppressWarnings(as.numeric(pieces[1]))
        x <- pieces[2]
    }

    pieces <- suppressWarnings(
        as.numeric(strsplit(x, ":", fixed = TRUE)[[1]])
    )

    if (anyNA(pieces) || is.na(days)) {
        return(NA_real_)
    }

    duration_minutes <- switch(
        as.character(length(pieces)),
        "3" = pieces[1] * 60 + pieces[2] + pieces[3] / 60,
        "2" = pieces[1] + pieces[2] / 60,
        "1" = pieces[1] / 60,
        NA_real_
    )

    days * 24 * 60 + duration_minutes
}

parse_iridis_out_file <- function(file) {
    if (length(file) != 1 || is.na(file) || !file.exists(file)) {
        return(
            tibble::tibble(
                file = file,
                state = NA_character_,
                mem_used_gb = NA_real_,
                mem_limit_gb = NA_real_,
                elapsed_minutes = NA_real_,
                time_limit_minutes = NA_real_
            )
        )
    }

    lines <- readLines(file, warn = FALSE)

    state_line <- lines[grepl("^State:", lines)]
    mem_eff_line <- lines[grepl("^Memory Efficiency:", lines)]
    mem_used_line <- lines[grepl("^Memory Utilized:", lines)]
    elapsed_line <- lines[grepl("^Elapsed time\\s*:", lines)]

    state <- if (length(state_line) > 0) {
        sub("^State:\\s*", "", state_line[length(state_line)])
    } else {
        NA_character_
    }

    mem_used <- if (length(mem_used_line) > 0) {
        sub(
            "^Memory Utilized:\\s*",
            "",
            mem_used_line[length(mem_used_line)]
        )
    } else {
        NA_character_
    }

    mem_limit <- if (length(mem_eff_line) > 0) {
        stringr::str_match(
            mem_eff_line[length(mem_eff_line)],
            "of\\s+([0-9.]+\\s*(KB|MB|GB|TB))"
        )[, 2]
    } else {
        NA_character_
    }

    elapsed <- if (length(elapsed_line) > 0) {
        stringr::str_match(
            elapsed_line[length(elapsed_line)],
            "^Elapsed time\\s*:\\s*([^ ]+)"
        )[, 2]
    } else {
        NA_character_
    }

    time_limit <- if (length(elapsed_line) > 0) {
        stringr::str_match(
            elapsed_line[length(elapsed_line)],
            "Timelimit=([^\\)]+)"
        )[, 2]
    } else {
        NA_character_
    }

    tibble::tibble(
        file = file,
        state = state,
        mem_used_gb = parse_memory_to_gb(mem_used),
        mem_limit_gb = parse_memory_to_gb(mem_limit),
        elapsed_minutes = parse_slurm_duration_minutes(elapsed),
        time_limit_minutes = parse_slurm_duration_minutes(time_limit)
    )
}

parse_iridis_err_file <- function(file) {
    if (length(file) != 1 || is.na(file) || !file.exists(file)) {
        return(
            tibble::tibble(
                err_file = file,
                err_has_oom_kill = FALSE,
                err_has_killed_rscript = FALSE,
                err_has_timeout = FALSE,
                err_text = NA_character_
            )
        )
    }

    lines <- readLines(file, warn = FALSE)

    has_oom_kill <- any(
        grepl(
            "oom_kill|OOM Killed|Out Of Memory|out of memory",
            lines,
            ignore.case = TRUE
        )
    )

    has_killed_rscript <- any(
        grepl("Killed\\s+Rscript", lines)
    )

    has_timeout <- any(
        grepl(
            "DUE TO TIME LIMIT|TIME LIMIT|TIMEOUT",
            lines,
            ignore.case = TRUE
        )
    )

    tibble::tibble(
        err_file = file,
        err_has_oom_kill = has_oom_kill,
        err_has_killed_rscript = has_killed_rscript,
        err_has_timeout = has_timeout,
        err_text = paste(lines, collapse = "\n")
    )
}

find_latest_iridis_log_files <- function(iridis_dir) {
    storage_dir <- file.path(iridis_dir, "storage")

    files <- list.files(
        storage_dir,
        pattern = "^run-[0-9]+-[0-9]+\\.(out|err)$",
        full.names = TRUE
    )

    parse_iridis_log_filename(files) %>%
        tidyr::pivot_wider(
            id_cols = c(array_job_id, process_id),
            names_from = type,
            values_from = file,
            names_prefix = "file_"
        ) %>%
        group_by(process_id) %>%
        slice_max(
            order_by = array_job_id,
            n = 1,
            with_ties = FALSE
        ) %>%
        ungroup()
}


find_oom_ids_iridis <- function(iridis_dir, min_mem_gb = 8,
                                mem_limit_if_missing_gb = NULL,
                                suspicious_mem_gb = 0.1) {
    latest_logs <- find_latest_iridis_log_files(iridis_dir)

    out_tab <- purrr::map_dfr(
        latest_logs$file_out,
        parse_iridis_out_file
    ) %>%
        dplyr::rename(out_file = file)

    err_tab <- purrr::map_dfr(
        latest_logs$file_err,
        parse_iridis_err_file
    )

    oom_tab <- latest_logs %>%
        dplyr::left_join(out_tab, by = c("file_out" = "out_file")) %>%
        dplyr::left_join(err_tab, by = c("file_err" = "err_file")) %>%
        dplyr::mutate(
            out_state_oom = grepl("^OUT_OF_MEMORY", state),
            err_oom = err_has_oom_kill | err_has_killed_rscript,
            confirmed_oom = out_state_oom | err_oom,

            mem_limit_gb_reported = mem_limit_gb,
            mem_limit_gb = dplyr::if_else(
                is.na(mem_limit_gb) & confirmed_oom &
                    !is.null(mem_limit_if_missing_gb),
                as.numeric(mem_limit_if_missing_gb),
                mem_limit_gb
            ),
            mem_limit_source = dplyr::case_when(
                !is.na(mem_limit_gb_reported) ~ "slurm_out",
                is.na(mem_limit_gb_reported) & confirmed_oom &
                    !is.null(mem_limit_if_missing_gb) ~ "fallback_argument",
                TRUE ~ NA_character_
            ),

            sufficient_mem_limit =
                !is.na(mem_limit_gb) & mem_limit_gb >= min_mem_gb,

            suspicious_low_mem_accounting =
                !is.na(mem_used_gb) & mem_used_gb < suspicious_mem_gb
        ) %>%
        dplyr::filter(
            confirmed_oom,
            sufficient_mem_limit
        )

    load(file.path(iridis_dir, "setup.Rda"))

    if (!("process_id" %in% names(setup))) {
        setup <- setup %>%
            dplyr::mutate(process_id = dplyr::row_number())
    }

    oom_tab %>%
        dplyr::left_join(setup, by = "process_id") %>%
        dplyr::arrange(process_id)
}

find_timeout_ids_iridis <- function(
    iridis_dir,
    min_time_minutes = 240,
    time_limit_if_missing_minutes = NULL
) {
    latest_logs <- find_latest_iridis_log_files(iridis_dir)

    out_tab <- purrr::map_dfr(
        latest_logs$file_out,
        parse_iridis_out_file
    ) %>%
        dplyr::rename(out_file = file)

    err_tab <- purrr::map_dfr(
        latest_logs$file_err,
        parse_iridis_err_file
    )

    fallback_time <- if (is.null(time_limit_if_missing_minutes)) {
        NA_real_
    } else {
        as.numeric(time_limit_if_missing_minutes)
    }

    timeout_tab <- latest_logs %>%
        dplyr::left_join(
            out_tab,
            by = c("file_out" = "out_file")
        ) %>%
        dplyr::left_join(
            err_tab,
            by = c("file_err" = "err_file")
        ) %>%
        dplyr::mutate(
            out_state_timeout = dplyr::coalesce(
                grepl("^TIMEOUT", state),
                FALSE
            ),

            err_timeout = dplyr::coalesce(
                err_has_timeout,
                FALSE
            ),

            confirmed_timeout =
                out_state_timeout | err_timeout,

            time_limit_minutes_reported =
                time_limit_minutes,

            time_limit_minutes = dplyr::case_when(
                !is.na(time_limit_minutes_reported) ~
                    time_limit_minutes_reported,

                confirmed_timeout & !is.na(fallback_time) ~
                    fallback_time,

                TRUE ~ NA_real_
            ),

            time_limit_source = dplyr::case_when(
                !is.na(time_limit_minutes_reported) ~
                    "slurm_out",

                confirmed_timeout & !is.na(fallback_time) ~
                    "fallback_argument",

                TRUE ~ NA_character_
            ),

            sufficient_time_limit =
                !is.na(time_limit_minutes) &
                time_limit_minutes >= min_time_minutes
        ) %>%
        dplyr::filter(
            confirmed_timeout,
            sufficient_time_limit
        )

    load(file.path(iridis_dir, "setup.Rda"))

    if (!("process_id" %in% names(setup))) {
        setup <- setup %>%
            dplyr::mutate(
                process_id = dplyr::row_number()
            )
    }

    timeout_tab %>%
        dplyr::left_join(setup, by = "process_id") %>%
        dplyr::arrange(process_id)
}

write_update_script_iridis <- function(simulation_fun, time_each, mem = 4,
                                       iridis_dir,
                                       extra_ids = NULL,
                                       exclude_ids = NULL,
                                       exclude_oom = FALSE,
                                       oom_min_mem_gb = 8) {
    missing_ids <- find_missing_ids(iridis_dir)

    if (exclude_oom) {
        oom_ids <- find_oom_ids_iridis(
            iridis_dir,
            min_mem_gb = oom_min_mem_gb
        ) %>%
            dplyr::pull(process_id)

        exclude_ids <- unique(c(exclude_ids, oom_ids))
    }

    ids <- sort(unique(c(missing_ids, extra_ids)))
    ids <- setdiff(ids, exclude_ids)

    if (length(ids) == 0) {
        stop("No missing or extra process IDs to run after exclusions.")
    }

    write_shell_script_iridis(
        ids,
        time_each,
        mem,
        iridis_dir,
        "simstudy_update.slurm"
    )

    invisible(ids)
}
