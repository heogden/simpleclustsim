make_pace_optns <- function(method_select_k = c("BIC", "AIC", "FVE", "Oracle"),
                            fve_threshold = NULL,
                            bandwidth = c("GCV", "default")) {
    method_select_k <- match.arg(method_select_k)
    bandwidth <- match.arg(bandwidth)

    optns <- list(dataType = "Sparse")

    if (bandwidth == "GCV") {
        optns$methodBwMu <- "GCV"
        optns$methodBwCov <- "GCV"
    }

    if (method_select_k != "Oracle") {
        optns$methodSelectK <- method_select_k
    }

    if (method_select_k == "FVE") {
        if (is.null(fve_threshold)) {
            stop("fve_threshold must be supplied when method_select_k = 'FVE'")
        }
        optns$FVEthreshold <- fve_threshold
    }

    optns
}


make_pace_spec <- function(method_select_k,
                           fve_threshold = NULL,
                           bandwidth = c("GCV", "default")) {
    bandwidth <- match.arg(bandwidth)

    list(
        family = "fdapace",
        oracle_k = method_select_k == "Oracle",
        optns = make_pace_optns(
            method_select_k = method_select_k,
            fve_threshold = fve_threshold,
            bandwidth = bandwidth
        )
    )
}

make_bayesFPCA_spec <- function(method_choose_L,
                                PVEthreshold = NULL,
                                Lmax = NULL,
                                oracle_L = FALSE,
                                lambda_L = 2,
                                use_nbasis_as_K = TRUE) {
    optns <- list(
        methodChooseL = method_choose_L,
        seed = NULL,
        K = NULL
    )

    if (!is.null(PVEthreshold)) {
        optns$PVEthreshold <- PVEthreshold
    }

    if (!is.null(Lmax)) {
        optns$Lmax <- Lmax
    }

    if (!is.null(lambda_L)) {
        optns$lambda_L <- lambda_L
    }

    list(
        family = "bayesFPCA",
        optns = optns,
        oracle_L = oracle_L,
        use_nbasis_as_K = use_nbasis_as_K
    )
}

make_face_spec <- function(pve = 0.99) {
    list(
        family = "face",
        pve = pve
    )
}

make_refund_spec <- function(refund_method = c("fpca.sc", "ccb.fpc"),
                             pve = 0.99,
                             oracle_k = FALSE) {
    refund_method <- match.arg(refund_method)

    list(
        family = "refund",
        refund_method = refund_method,
        pve = pve,
        oracle_k = oracle_k
    )
}

make_local_gam_spec <- function(k = NULL, bs = "cr", method = "REML") {
    list(
        family = "local_gam",
        k = k,
        bs = bs,
        method = method
    )
}

method_specs <- function() {
    list(
        ## Canonical methods used in the full simulations
        "PACE" = make_pace_spec(
            method_select_k = "FVE",
            fve_threshold = 0.99,
            bandwidth = "GCV"
        ),

        "FACE" = make_face_spec(
            pve = 0.99
        ),

        "bayesFPCA" = make_bayesFPCA_spec(
            method_choose_L = "fixed",
            oracle_L = TRUE
        ),

        "Local-GAM" = make_local_gam_spec(),

        ## Existing aliases, so old simulation scripts still work
        "PACE-BIC" = make_pace_spec("BIC", bandwidth = "GCV"),
        "PACE-AIC" = make_pace_spec("AIC", bandwidth = "GCV"),
        "PACE-95" = make_pace_spec("FVE", fve_threshold = 0.95, bandwidth = "GCV"),
        "PACE-99" = make_pace_spec("FVE", fve_threshold = 0.99, bandwidth = "GCV"),
        "PACE-Oracle" = make_pace_spec("Oracle", bandwidth = "GCV"),

        ## New pilot variants
        "PACE-BIC-GCV" = make_pace_spec("BIC", bandwidth = "GCV"),
        "PACE-AIC-GCV" = make_pace_spec("AIC", bandwidth = "GCV"),
        "PACE-FVE95-GCV" = make_pace_spec("FVE", fve_threshold = 0.95, bandwidth = "GCV"),
        "PACE-FVE99-GCV" = make_pace_spec("FVE", fve_threshold = 0.99, bandwidth = "GCV"),

        "PACE-BIC-default" = make_pace_spec("BIC", bandwidth = "default"),
        "PACE-AIC-default" = make_pace_spec("AIC", bandwidth = "default"),
        "PACE-FVE95-default" = make_pace_spec("FVE", fve_threshold = 0.95, bandwidth = "default"),
        "PACE-FVE99-default" = make_pace_spec("FVE", fve_threshold = 0.99, bandwidth = "default"),

        "PACE-Oracle-GCV" = make_pace_spec("Oracle", bandwidth = "GCV"),
        "PACE-Oracle-default" = make_pace_spec("Oracle", bandwidth = "default"),
        
        ## Existing aliases, but fixed
        "bayesFPCA-95" = make_bayesFPCA_spec(
            method_choose_L = "PVE",
            PVEthreshold = 95,
            Lmax = 15
        ),

        "bayesFPCA-99" = make_bayesFPCA_spec(
            method_choose_L = "PVE",
            PVEthreshold = 99,
            Lmax = 15
        ),

        "bayesFPCA-Oracle" = make_bayesFPCA_spec(
            method_choose_L = "fixed",
            oracle_L = TRUE
        ),

        ## Clearer pilot names
        "bayesFPCA-PVE95-L15" = make_bayesFPCA_spec(
            method_choose_L = "PVE",
            PVEthreshold = 95,
            Lmax = 15
        ),

        "bayesFPCA-PVE99-L15" = make_bayesFPCA_spec(
            method_choose_L = "PVE",
            PVEthreshold = 99,
            Lmax = 15
        ),

        "bayesFPCA-model-choice-L15" = make_bayesFPCA_spec(
            method_choose_L = "model_choice",
            Lmax = 15,
            lambda_L = 2
        ),
        
        "FACE-99" = make_face_spec(pve = 0.99),
        "FACE-95" = make_face_spec(pve = 0.95),

        "Di-Oracle" = make_refund_spec(
            refund_method = "fpca.sc",
            oracle_k = TRUE
        ),
        "Di-95" = make_refund_spec(
            refund_method = "fpca.sc",
            pve = 0.95
        ),
        "Di-99" = make_refund_spec(
            refund_method = "fpca.sc",
            pve = 0.99
        ),
        "Goldsmith-95" = make_refund_spec(
            refund_method = "ccb.fpc",
            pve = 0.95
        ),
        "Goldsmith-99" = make_refund_spec(
            refund_method = "ccb.fpc",
            pve = 0.99
        )
    )
}


get_method_spec <- function(method) {
    specs <- method_specs()

    if (method %in% names(specs)) {
        specs[[method]]
    } else {
        NULL
    }
}


final_methods <- function(case) {
    general_methods <- c(
        "AdaStruMM",
        "Local-GAM",
        "HGAM-GS",
        "PACE",
        "FACE"
    )

    switch(
        case,

        "2re" = c(
            general_methods,
            "bayesFPCA",
            "Oracle-GP"
        ),

        "3re" = c(
            general_methods,
            "bayesFPCA",
            "Oracle-GP"
        ),

        "sitar" = c(
            general_methods,
            "SITAR",
            "Oracle-SITAR"
        ),

        stop("No final method list defined for case ", case)
    )
}

method_labels <- function() {
    c(
        "AdaStruMM" = "AdaStruMM",
        "Local-GAM" = "Local-GAM",
        "HGAM-GS" = "HGAM-GS",
        "PACE" = "PACE",
        "FACE" = "FACE",
        "bayesFPCA" = "bayesFPCA (oracle K)",
        "SITAR" = "SITAR",
        "Oracle-GP" = "Individual oracle bound",
        "Oracle-SITAR" = "Individual oracle bound"
    )
}

method_label <- function(method) {
    labels <- method_labels()

    out <- unname(labels[method])
    out[is.na(out)] <- method[is.na(out)]

    out
}

fit_mod_from_spec <- function(data, spec, k = NA, nbasis = 10, seed = NULL,
                              z_fun = NULL, oracle_spec = NULL) {
    switch(
        spec$family,

        "fdapace" = {
            k_fdapace <- if (isTRUE(spec$oracle_k)) k else NA

            simpleclust::fit_mod_fdapace(
                data = data,
                k = k_fdapace,
                optns = spec$optns
            )
        },

        "bayesFPCA" = {
            optns <- spec$optns

            optns$seed <- seed

            if (isTRUE(spec$use_nbasis_as_K)) {
                optns$K <- nbasis
            }

            if (isTRUE(spec$oracle_L)) {
                optns$L <- k
            }

            simpleclust::fit_mod_bayesFPCA(
                data = data,
                optns = optns
            )
        },

        "face" = {
            simpleclust::fit_mod_face(
                             data = data,
                             knots = nbasis - 3,
                             pve = spec$pve
                         )
        },

        "refund" = {
            k_refund <- if (isTRUE(spec$oracle_k)) k else NA
            
            simpleclust::fit_mod_refund(
                             data = data,
                             method = spec$refund_method,
                             k = k_refund,
                             pve = spec$pve,
                             nbasis = nbasis
                         )
        },

        "local_gam" = {
            k_local <- if (is.null(spec$k)) nbasis else spec$k

            simpleclust::fit_mod_local_gam(
                             data = data,
                             k = k_local,
                             bs = spec$bs,
                             method = spec$method
                         )
        },

        stop("No fitting method defined for family ", spec$family)
    )
}
