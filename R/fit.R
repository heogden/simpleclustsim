#' could use switch to give different options by case
find_nbasis <- function(case, method) {
    if(method == "SITAR")
        nbasis <- 5
    else
        nbasis <- 10
    nbasis
}



fit_mod <- function(data, method, k = NA, nbasis = 10, seed = NULL,
                    z_fun = NULL, oracle_spec = NULL, true_u = NULL) {
    spec <- get_method_spec(method)

    if (!is.null(spec)) {
        return(
            fit_mod_from_spec(
                data = data,
                spec = spec,
                k = k,
                nbasis = nbasis,
                seed = seed,
                z_fun = z_fun,
                oracle_spec = oracle_spec
            )
        )
    }
    switch(method,
           "AdaStruMM" = simpleclust::fit_mod_adastrumm(data, nbasis = nbasis),
           "LMM-RI" = simpleclust::fit_mod_lmer(data, random_slope = FALSE),
           "LMM-RS" = simpleclust::fit_mod_lmer(data, random_slope = TRUE),
           "HGAM-GS" = simpleclust::fit_mod_mgcv(data, nbasis = nbasis, bam = TRUE),
           "SITAR" = simpleclust::fit_mod_sitar(data, nbasis = 5),
           "GAMM-RI" = simpleclust::fit_mod_gamm(data,
                                                 nbasis = nbasis,
                                                 z_fun = function(x){matrix(1, ncol = 1, nrow = length(x))},
                                                 re_structure = "correlated",
                                                 control = nlme::lmeControl(maxIter = 1000,
                                                                            msMaxIter = 1000,
                                                                            niterEM = 200)),
           "GAMM-RS" = simpleclust::fit_mod_gamm(data,
                                                 nbasis = nbasis,
                                                 z_fun = function(x){cbind(1, x)},
                                                 re_structure = "correlated",
                                                 control = nlme::lmeControl(maxIter = 500,
                                                                            msMaxIter = 500,
                                                                            niterEM = 200)),
           "Oracle-MM" = simpleclust::fit_mod_gamm(data,
                                                   nbasis = nbasis,
                                                   z_fun = z_fun,
                                                   re_structure = "correlated",
                                                   control = nlme::lmeControl(maxIter = 1000,
                                                                              msMaxIter = 1000,
                                                                              niterEM = 200)),
           "Oracle-MM-D" = simpleclust::fit_mod_gamm(data,
                                                     nbasis = nbasis,
                                                     z_fun = z_fun,
                                                     re_structure = "independent",
                                                     control = nlme::lmeControl(maxIter = 1000,
                                                                                msMaxIter = 1000,
                                                                                niterEM = 200)),
           
           "FACE" = simpleclust::fit_mod_face(data, knots = nbasis - 3,
                                              pve = 0.99),
           "Oracle-GP" = simpleclust::fit_mod_oracle_known_gp(
                                          data = data,
                                          oracle_spec = oracle_spec
                                      ),
           "Oracle-SITAR" = simpleclust::fit_mod_oracle_known_sitar(
                                             data = data,
                                             oracle_spec = oracle_spec,
                                             true_u = true_u,
                                             use_true_u_start = TRUE,
                                             n_draws = 5000,
                                             seed = seed,
                                             proposal_scale = c(0.5, 0.75, 1, 1.5, 2),
                                             mode_logpost_drop = 30
                                         ),
           stop("method ", method, " not found"))
}

find_package <- function(method) {
    spec <- get_method_spec(method)

    if (!is.null(spec)) {
        return(spec$family)
    }
    
    switch(method,
           "Di-Oracle" = "refund",
           "Di-95" = "refund",
           "Di-99" = "refund",
           "Goldsmith-Oracle" = "refund",
           "Goldsmith-95" = "refund",
           "Goldsmith-99" = "refund",
           "AdaStruMM" = "adastrumm",
           "LMM-RI" = "lme4",
           "LMM-RS" = "lme4",
           "HGAM-GS" = "mgcv",
           "SITAR" = "sitar",
           "GAMM-RI" = "mgcv_gamm",
           "GAMM-RS" = "mgcv_gamm",
           "Oracle-GP" = "oracle_known_gp",
           "Oracle-SITAR" = "oracle_known_sitar",
           "Oracle-MM" = "mgcv_gamm",
           "Oracle-MM-D" = "mgcv_gamm",
           "FACE" = "face",
           stop("method ", method, " not found"))     

}


predict_eta <- function(mod, method, c, x, interval = TRUE) {

    package <- find_package(method)
    switch(package,
           "fdapace" = simpleclust::predict_eta_fdapace(mod, c, x, interval),
           "adastrumm" = simpleclust::predict_eta_adastrumm(mod, c, x, interval),
           "refund" = simpleclust::predict_eta_refund(mod, c, x, interval),
           "lme4" = simpleclust::predict_eta_lmer(mod, c, x, interval),
           "mgcv" = simpleclust::predict_eta_mgcv(mod, c, x, interval),
           "mgcv_gamm" = simpleclust::predict_eta_gamm(mod, c, x, interval),
           "sitar" = simpleclust::predict_eta_sitar(mod, c, x, interval),
           "bayesFPCA" = simpleclust::predict_eta_bayesFPCA(mod, c, x, interval),
           "face" = simpleclust::predict_eta_face(mod, c, x, interval),
           "oracle_known_gp" = simpleclust::predict_eta_oracle_known_gp(
                                                mod, c, x, interval
                                            ),
           "oracle_known_sitar" = simpleclust::predict_eta_oracle_known_sitar(
                                                   mod, c, x, interval
                                               ),
           "local_gam" = simpleclust::predict_eta_local_gam(mod, c, x, interval),
           stop("prediction for package ", package, " not found"))
}

get_k_hat <- function(mod, method) {
    package <- find_package(method)
    switch(package,
           "fdapace" = mod$selectK,
           "adastrumm" = mod$k,
           "refund" = if (!is.null(mod$npc)) as.numeric(mod$npc) else NA_real_,
           "lme4" = ncol(coef(mod)$c),
           "mgcv" = NA,
           "mgcv_gamm" = length(mod$z_names),
           "sitar" = NA,
           "bayesFPCA" = mod$L,
           "face" = length(mod$mod$eigenvalues),
           "oracle_known_gp" = mod$oracle_spec$K,
           "local_gam" = NA_real_,
           stop("k_hat for package ", package, " not found"))
    
}

