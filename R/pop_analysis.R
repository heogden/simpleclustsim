get_GP_emp_pred <- function(pred_data, estimate = FALSE) {
    if(estimate)
        if(length(colnames(pred_data$mu_c_hat)) == 0)
            return(NA)
       
    if(estimate)
        pred_data <- pred_data %>% select(-mu_c) %>% mutate(mu_c = mu_c_hat$estimate)
    
    pa_data <- pred_data %>%
        group_by(x) %>%
        summarise(mu = mean(mu_c))
    
    m <- pa_data$mu
    
    r_c_data <- pred_data %>%
        left_join(pa_data, by = "x") %>%
        mutate(r_c = mu_c - mu) %>%
        select(x, c, r_c)

    
    x_poss <- unique(pred_data$x)

    r_c_data_s <- r_c_data %>%
        rename(s = x, r_c_s = r_c)

    r_c_data_t <- r_c_data %>%
        rename(t = x, r_c_t = r_c)
    

    cov_data <- crossing(s = unique(pred_data$x),
                         t = unique(pred_data$x),
                         c = unique(pred_data$c)) %>%
        left_join(r_c_data_s, by = c("s", "c")) %>%
        left_join(r_c_data_t, by = c("t", "c")) %>%
        mutate(r_c_prod = r_c_s * r_c_t) %>%
        group_by(s, t) %>%
        summarise(cov = mean(r_c_prod))

    list(x_poss = x_poss,
         m = pa_data$mu,
         C = matrix(cov_data$cov, nrow = length(x_poss)))
}


get_GP_hat <- function(method, mod, x, c_poss) {
    if(is.null(mod))
        return(NA)

    pred_fun <- function(mod, c, x) {
        predict_eta(mod, method, c, x, interval = FALSE)
    }
    
    package <- find_package(method)
    switch(package,
           "fdapace" = simpleclust::get_GP_fdapace(mod, x),
           "adastrumm" = simpleclust::get_GP_adastrumm(mod, x),
           "face" = simpleclust::get_GP_face(mod, x),
           "bayesFPCA" = simpleclust::get_GP_bayesFPCA(mod, x),
           "mgcv_gamm" = simpleclust::get_GP_gamm(mod, x),
           "oracle_known_gp" = list(Individual = NA_real_),
           "oracle_known_sitar" = list(Individual = NA_real_),
           list(Empirical = simpleclust::get_GP_emp(pred_fun, c_poss, mod)(x)))
}

sqrt_degen <- function(A) {
    tol <- 1e-8
    eA <- eigen(A)
    r <- sum(eA$values > tol)
    Q <- eA$vectors[, 1:r, drop = FALSE]
    sqrt_Lambda <- diag(sqrt(eA$values[1:r]), nrow = r, ncol = r)
    Q %*% sqrt_Lambda %*% t(Q)
}




find_d_m_bar_GP <- function(GP0, GP_hat) {
    if(length(GP_hat) == 0 || any(is.na(GP_hat)))
        return(NA)

    result <- NA
    try(result <- mean((GP0$m - GP_hat$m)^2))
    result
}

find_d_C_bar <- function(C, C_hat) {
    N <- nrow(C)
    sqrt_C <- sqrt_degen(C)
    A <- sqrt_C %*% C_hat %*% sqrt_C
    sqrt_A <- sqrt_degen(A)
    B <- C_hat + C - 2 * sqrt_A
    tr_B <- sum(diag(B))
    tr_B / N
}

find_d_C_bar_GP <- function(GP0, GP_hat) {
    if(length(GP_hat) == 0 || any(is.na(GP_hat)))
        return(NA)

    result <- NA
    try(result <- find_d_C_bar(GP0$C, GP_hat$C))
    result
}



