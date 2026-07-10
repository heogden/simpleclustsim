as_error_chr <- function(x) {
    if (is.null(x) || length(x) == 0) {
        NA_character_
    } else {
        as.character(x[[1]])
    }
}

run_simstudy_each <- function(case, subcase_id, method, seed, subcases, ...) {
    arg_names <- c("case", "subcase_id", "method", "seed")
    args <- c(case, subcase_id, method, seed)
    names(args) <- arg_names

    cat(paste(arg_names, args, sep = " = "), sep = ", ")
    cat("\n")

    total_time <- system.time({

        sim_res <- safe_eval(
            simulate(case, subcases[subcase_id, ], seed)
        )

        data_full <- sim_res$value

        data <- NULL
        pred_data <- NULL
        GP0 <- NULL
        GP0_emp <- NULL
        GP_hat <- list(Unavailable = NA_real_)
        mod <- NULL
        k_hat <- NA_real_

        fit_res <- list(
            value = NULL,
            status = "skipped",
            error = NA_character_,
            warnings = character(),
            n_warnings = 0,
            time = NA_real_
        )

        pred_res <- list(
            value = NULL,
            status = "skipped",
            error = NA_character_,
            warnings = character(),
            n_warnings = 0,
            time = NA_real_
        )

        gp_res <- list(
            value = NULL,
            status = "skipped",
            error = NA_character_,
            warnings = character(),
            n_warnings = 0,
            time = NA_real_
        )

        k_res <- list(
            value = NA,
            status = "skipped",
            error = NA_character_,
            warnings = character(),
            n_warnings = 0,
            time = NA_real_
        )

        k_hat <- if (k_res$status == "ok" && !is.null(k_res$value)) {
                     as.numeric(k_res$value)
                 } else {
                     NA_real_
                 }
        
        if (sim_res$status == "ok") {
            data <- data_full$data
            GP0 <- data_full$GP0
            GP0_emp <- data_full$GP0_emp

            fit_res <- safe_eval(
                fit_mod(
                    data,
                    method,
                    k = find_k0(case),
                    nbasis = find_nbasis(case, method),
                    seed = seed,
                    z_fun = data_full$z_fun,
                    oracle_spec = data_full$oracle_spec,
                    true_u = data_full$true_u
                )
            )

            mod <- fit_res$value

            if (!is.null(mod)) {
                pred_res <- safe_eval(
                    data_full$pred_data %>%
                        mutate(
                            mu_c_hat = predict_eta(
                                mod,
                                !!method,
                                c,
                                x,
                                interval = TRUE
                            )
                        )
                )

                pred_data <- pred_res$value

                x_poss <- unique(data_full$pred_data$x)
                c_poss <- unique(data_full$pred_data$c)

                gp_res <- safe_eval(
                    get_GP_hat(method, mod, x_poss, c_poss)
                )

                GP_hat <- if (gp_res$status == "ok" && !is.null(gp_res$value)) {
                              gp_res$value
                          } else {
                              list(Unavailable = NA_real_)
                          }
                

                k_res <- safe_eval(
                    get_k_hat(mod, method)
                )
                
                k_hat <- if (k_res$status == "ok" &&
                             !is.null(k_res$value) &&
                             length(k_res$value) > 0) {
                             as.numeric(k_res$value[1])
                         } else {
                             NA_real_
                         }
            }
        }
    })

    tibble(
        case = case,
        subcase_id = subcase_id,
        method = method,
        seed = seed,

        data = list(data),
        pred_data = list(pred_data),
        GP0 = list(GP0),
        GP0_emp = list(GP0_emp),
        GP_hat = list(GP_hat),
        k_hat = k_hat,

        time = unname(total_time["elapsed"]) / 60,
        time_sim = sim_res$time,
        time_fit = fit_res$time,
        time_pred = pred_res$time,
        time_gp = gp_res$time,

        sim_status = sim_res$status,
        fit_status = fit_res$status,
        pred_status = pred_res$status,
        gp_status = gp_res$status,
        k_status = k_res$status,

        sim_error = as_error_chr(sim_res$error),
        fit_error = as_error_chr(fit_res$error),
        pred_error = as_error_chr(pred_res$error),
        gp_error = as_error_chr(gp_res$error),
        k_error = as_error_chr(k_res$error),

        sim_warnings = list(sim_res$warnings),
        fit_warnings = list(fit_res$warnings),
        pred_warnings = list(pred_res$warnings),
        gp_warnings = list(gp_res$warnings),
        k_warnings = list(k_res$warnings),

        n_sim_warnings = sim_res$n_warnings,
        n_fit_warnings = fit_res$n_warnings,
        n_pred_warnings = pred_res$n_warnings,
        n_gp_warnings = gp_res$n_warnings,
        n_k_warnings = k_res$n_warnings
    )
}
