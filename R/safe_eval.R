safe_eval <- function(expr) {
    warnings <- character()
    error <- NA_character_
    value <- NULL
    has_error <- FALSE

    time <- system.time({
        value <- withCallingHandlers(
            tryCatch(
                expr,
                error = function(e) {
                    has_error <<- TRUE
                    error <<- conditionMessage(e)
                    NULL
                }
            ),
            warning = function(w) {
                warnings <<- c(warnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            }
        )
    })

    status <- if (has_error) "error" else "ok"

    list(
        value = value,
        status = status,
        error = error,
        warnings = warnings,
        n_warnings = length(warnings),
        time = unname(time["elapsed"]) / 60
    )
}
