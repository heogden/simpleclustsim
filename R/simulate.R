find_Cov_re_from_sitar <- function(mod) {
    sigma <- mod$sigma

    s <- summary(mod$modelStruct)
    s1 <- s[[1]]
    s11 <- s1[[1]]


    ss11 <- summary(s11)
    sds <- attr(ss11, "stdDev")
    sigma_re <- sds * sigma

    ll <- lower.tri(ss11)
    corr_re <- ss11[ll]
    
    sigma_a <- sigma_re[1]
    sigma_b <- sigma_re[2]
    sigma_c <- sigma_re[3]
    r_ab <- corr_re[1]
    r_ac <- corr_re[2]
    r_bc <- corr_re[3]
    S_ab <- sigma_a * sigma_b * r_ab
    S_ac <- sigma_a * sigma_c * r_ac
    S_bc <- sigma_b * sigma_c * r_bc
    matrix(c(sigma_a^2, S_ab, S_ac,
             S_ab, sigma_b^2, S_bc,
             S_ac, S_bc, sigma_c^2),
           nrow = 3, ncol = 3)
}


find_z_fun_from_GP <- function(GP, k_max = 10) {
    C <- GP$C
    x_poss <- GP$x_poss
    dx <- GP$x_poss[2] - GP$x_poss[1]
    A <- C * dx
    ee <- eigen(A, symmetric = TRUE)
    lambda <- ee$values[1:k_max]
    z_unnorm <- ee$vectors[, 1:k_max, drop = FALSE]
    z_x_poss <- apply(z_unnorm, 2, function(v) v / sqrt(sum(v^2) * dx))

    z_fun <- function(x) {
        zi_list <- lapply(1:ncol(z_x_poss), function(i) {
            stats::splinefun(x_poss, z_x_poss[,i])(x)
        })
        Reduce(cbind, zi_list)
    }
    z_fun
}

simulate_sitar <- function(seed, h, x_range, mean_re, Cov_re, sigma, n_clusters, n_obs_per_cluster, GP0) {
    set.seed(seed)
 
    c <- rep(1:n_clusters, each = n_obs_per_cluster)
    j <- rep(1:n_obs_per_cluster, times = n_clusters)

    gap <- (x_range[2] - x_range[1]) / n_obs_per_cluster
    x1 <- runif(n_clusters, min = x_range[1], max = x_range[1] + gap)
    x <- rep(x1, each = n_obs_per_cluster) + (j - 1) * gap
    
    u <- mvtnorm::rmvnorm(n_clusters, mean = mean_re, sigma = Cov_re)
    alpha <- u[,1]
    beta <- u[,2]
    gamma <- u[,3]

    true_u <- tibble::tibble(
                          c = seq_len(n_clusters),
                          alpha = alpha,
                          beta = beta,
                          gamma = gamma
                      )

    mu_c <- function(x, c) {
        alpha[c] + h((x - beta[c])/exp(-gamma[c]))
    }

    if(length(GP0) > 0) {
        x_poss <- GP0$x_poss
    } else {
        x_poss <- seq(min(x), max(x), length.out = 100)
    }
    
    
    pred_data <- tidyr::crossing(x = x_poss,
                          c = 1:n_clusters) %>%
        dplyr::mutate(mu_c = mu_c(x, c))
   
    mu <- alpha[c] + h((x - beta[c])/exp(-gamma[c]))
    epsilon <- rnorm(length(mu), sd = sigma)
    
    y <- mu + epsilon

    data <- tibble(c = c,
                   x = x,
                   y = y,
                   mu = mu)

    GP0_emp <- get_GP_emp_pred(pred_data, estimate = FALSE)

    if(length(GP0) > 0)
        z_fun <- find_z_fun_from_GP(GP0)
    else
        z_fun <- NULL
    
    list(
        data = data,
        pred_data = pred_data,
        GP0 = GP0,
        GP0_emp = GP0_emp,
        z_fun = z_fun,
        oracle_spec = simpleclust::make_oracle_sitar_spec(
                                       h = h,
                                       mean_re = mean_re,
                                       Cov_re = Cov_re,
                                       sigma = sigma
                                   ),
        true_u = true_u
    )
    
}

init_Berkeley <- function() {
    berkeley <- sitar::berkeley
    ff <- na.omit(berkeley[berkeley$sex == 2 & berkeley$age >= 8 & berkeley$age <= 18, 
                           c('id', 'age', 'height')])
    fh1 <- sitar::sitar(x = age, y = height, id = id, data = ff, df = 5)

    theta_hat <- nlme::fixef(fh1)
    s <- as.numeric(theta_hat[1:5])
    
    knots <-  c(-3, -1, 1, 3)
    bounds <- c(-5.4, 5.4)

    h <- function(x_star) {
        as.numeric(tcrossprod(s, splines::ns(x_star,k=knots,B=bounds)))
    }


    mean_re <- as.numeric(theta_hat[6:8])
    Cov_re <- find_Cov_re_from_sitar(fh1)

    data_full <- simulate_sitar(1, h = h, x_range = c(-5, 5), mean_re = mean_re,
                                Cov_re = Cov_re, sigma = fh1$sigma,
                                n_clusters = 5000,
                                n_obs_per_cluster = 1,
                                GP0 = NULL)

    list(GP0 =  data_full$GP0_emp,
         h = h,
         x_range = c(-5, 5),
         mean_re = mean_re,
         Cov_re = Cov_re,
         sigma = fh1$sigma)
}


simulate_Berkeley <- function(seed, n_clusters, n_obs_per_cluster) {
    init <- init_Berkeley()
    
    simdata <- simulate_sitar(seed, h = init$h,
                              x_range = init$x_range,
                              mean_re = init$mean_re,
                              Cov_re = init$Cov_re,
                              sigma = init$sigma,
                              n_clusters = n_clusters,
                              n_obs_per_cluster = n_obs_per_cluster,
                              GP0 = init$GP0)
    
    #simdata$data$x <- simdata$data$x + 13 # back to original x scale
    #simdata$pred_data$x <- simdata$pred_data$x + 13
    #simdata$GP0_emp$x_poss <- simdata$GP0_emp$x_poss + 13

    simdata
}



simulate_2re <- function(seed, n_clusters, n_obs_per_cluster, sigma) {
    set.seed(seed)
    c <- rep(1:n_clusters, each = n_obs_per_cluster)
    x <- runif(length(c), 0, 3*pi)

    u1 <- rnorm(n_clusters)
    u2 <- rnorm(n_clusters)

    true_u <- tibble::tibble(
                          c = seq_len(n_clusters),
                          u1 = u1,
                          u2 = u2
                      )
    
    mu_c <- function(x, c) {
        (1 + u1[c]) * (x/2 + sin(x)) + u2[c]
    }

    pred_data <- tidyr::crossing(x = seq(min(x), max(x), length.out = 100),
                                 c = 1:n_clusters) %>%
        dplyr::mutate(mu_c = mu_c(x, c))
   
    mu <- (1 + u1[c]) * (x/2 + sin(x)) + u2[c]
    epsilon <- rnorm(length(mu), sd = sigma)
    
    y <- mu + epsilon

    data <- tibble(c = c,
                   x = x,
                   y = y,
                   mu = mu)

    x_poss <- unique(pred_data$x)

    f0 <- function(x) { x/2 + sin(x) }
    h <- function(x) { x/2 + sin(x) }
    f <- function(x) {
        cbind(1, h(x))
    }
    
    GP0 <- simpleclust::get_GP_fpc(f0, f)(x_poss)
    GP0_emp <- get_GP_emp_pred(pred_data, estimate = FALSE)

    oracle_spec <- simpleclust::make_oracle_gp_from_fpc(
                                    f0 = f0,
                                    f = f,
                                    sigma = sigma,
                                    K = 2
                                )

    list(
        data = data,
        pred_data = pred_data,
        GP0 = GP0,
        GP0_emp = GP0_emp,
        z_fun = f,
        oracle_spec = oracle_spec,
        true_u = true_u
    )
}


simulate_3re <- function(seed, n_clusters, n_obs_per_cluster, sigma) {
    set.seed(seed)
    c <- rep(1:n_clusters, each = n_obs_per_cluster)
    x <- runif(length(c), 0, 3*pi)

    u1 <- rnorm(n_clusters)
    u2 <- rnorm(n_clusters)
    u3 <- rnorm(n_clusters)

    true_u <- tibble::tibble(
                          c = seq_len(n_clusters),
                          u1 = u1,
                          u2 = u2,
                          u3 = u3
                      )

    mu_c <- function(x, c) {
        (1 + u1[c]) * (x/2 + sin(x)) + u2[c] + u3[c] * cos(x)
    }

    pred_data <- tidyr::crossing(x = seq(min(x), max(x), length.out = 100),
                                 c = 1:n_clusters) %>%
        dplyr::mutate(mu_c = mu_c(x, c))
   
    mu <- (1 + u1[c]) * (x/2 + sin(x)) + u2[c] + u3[c] * cos(x)
    epsilon <- rnorm(length(mu), sd = sigma)
    
    y <- mu + epsilon

    data <- tibble(c = c,
                   x = x,
                   y = y,
                   mu = mu)

    x_poss <- unique(pred_data$x)

    f0 <- function(x) { x/2 + sin(x) }
    h <- function(x) { x/2 + sin(x) }
    f <- function(x) {
        cbind(1, h(x), cos(x))
    }
    
    GP0 <- simpleclust::get_GP_fpc(f0, f)(x_poss)
    GP0_emp <- get_GP_emp_pred(pred_data, estimate = FALSE)

    oracle_spec <- simpleclust::make_oracle_gp_from_fpc(
                                    f0 = f0,
                                    f = f,
                                    sigma = sigma,
                                    K = 3
                                )

    list(
        data = data,
        pred_data = pred_data,
        GP0 = GP0,
        GP0_emp = GP0_emp,
        z_fun = f,
        oracle_spec = oracle_spec,
        true_u = true_u
    )
}


simulate <- function(case, subcase, seed) {
    required_cols <- c("n_clusters", "n_obs_per_cluster")
    missing_cols <- setdiff(required_cols, names(subcase))

    if (length(missing_cols) > 0) {
        stop(
            "subcase must contain columns ",
            paste(required_cols, collapse = ", "),
            "; missing ",
            paste(missing_cols, collapse = ", ")
        )
    }
    switch(case,
           "sitar" = simulate_Berkeley(seed,
                                       n_clusters = subcase$n_clusters,
                                       n_obs_per_cluster = subcase$n_obs_per_cluster),
           "2re" = simulate_2re(seed, n_clusters = subcase$n_clusters,
                                n_obs_per_cluster = subcase$n_obs_per_cluster,
                                sigma = 0.1),
           "3re" = simulate_3re(seed, n_clusters = subcase$n_clusters,
                                n_obs_per_cluster = subcase$n_obs_per_cluster,
                                sigma = 0.1),
           stop("case ", case, " not found"))
}


find_k0 <- function(case) {
    switch(case,
           "2re" = 2,
           "3re" = 3,
           "sitar" = NA,
           stop("case ", case, " not found"))
}


