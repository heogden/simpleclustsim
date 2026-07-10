project_root <- here::here()
devtools::load_all(project_root)

case_name <- "3re"

analysis_dir <- here::here(
    "instances",
    case_name,
    "analysis",
    "output"
)

figure_dir <- here::here(
    "instances",
    case_name,
    "analysis",
    "output",
    "figures"
)

dir.create(
    figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

analysis_results <- readRDS(
    file.path(
        analysis_dir,
        "analysis_results.rds"
    )
)

if (!identical(analysis_results$case, case_name)) {
    stop(
        "Loaded analysis_results for case ",
        analysis_results$case,
        ", but this script is for case ",
        case_name,
        ".",
        call. = FALSE
    )
}

plots <- make_case_simstudy_plots(
    analysis_results = analysis_results,
    legend = FALSE,
    tex = TRUE,
    show_interval = TRUE,
    log_time = FALSE,
    plot_exclusions_file = here::here(
        "instances",
        case_name,
        "analysis",
        "output",
        "plot_exclusions.csv"
    )
)

save_case_simstudy_plots(
    plots = plots,
    case_name = case_name,
    output_dir = figure_dir,
    formats = c("tikz", "pdf")
)
