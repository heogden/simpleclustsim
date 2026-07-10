## Reproduce paper figure panels from saved simulation summaries.
##
## This script does not use the raw checkpointed simulation outputs. It reads
## the processed CSV summaries in instances/<case>/analysis/output/ and
## regenerates the TikZ figure panels used for the paper.
##
## From the repository root, run:
##
##     Rscript scripts/reproduce_paper_figures.R

required_packages <- c("devtools", "readr")

missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_packages) > 0L) {
    stop(
        "Please install the following packages before running this script: ",
        paste(missing_packages, collapse = ", "),
        ".",
        call. = FALSE
    )
}

find_project_root <- function() {
    root <- normalizePath(getwd(), mustWork = TRUE)

    if (
        !file.exists(file.path(root, "DESCRIPTION")) ||
        !dir.exists(file.path(root, "instances"))
    ) {
        stop(
            "Run this script from the simpleclustsim repository root.",
            call. = FALSE
        )
    }

    root
}

make_analysis_results_from_csv <- function(case_name, project_root) {
    analysis_output_dir <- file.path(
        project_root,
        "instances",
        case_name,
        "analysis",
        "output"
    )

    method_plot_file <- file.path(analysis_output_dir, "paper_plot_data.csv")
    oracle_file <- file.path(analysis_output_dir, "oracle_reference_data.csv")

    missing_files <- c(method_plot_file, oracle_file)[
        !file.exists(c(method_plot_file, oracle_file))
    ]

    if (length(missing_files) > 0L) {
        stop(
            "Missing saved summary file(s) for case ",
            case_name,
            ": ",
            paste(missing_files, collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    list(
        case = case_name,
        method_plot_data = readr::read_csv(
            method_plot_file,
            show_col_types = FALSE
        ),
        oracle_reference_data = readr::read_csv(
            oracle_file,
            show_col_types = FALSE
        )
    )
}

reproduce_case_figures <- function(case_name, project_root, formats = c("tikz")) {
    message("Reproducing figures for case: ", case_name)

    analysis_results <- make_analysis_results_from_csv(
        case_name = case_name,
        project_root = project_root
    )

    figure_dir <- file.path(
        project_root,
        "instances",
        case_name,
        "analysis",
        "output",
        "figures"
    )

    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    plot_exclusions_file <- file.path(
        project_root,
        "instances",
        case_name,
        "analysis",
        "output",
        "plot_exclusions_regenerated.csv"
    )

    plots <- make_case_simstudy_plots(
        analysis_results = analysis_results,
        legend = FALSE,
        tex = TRUE,
        show_interval = TRUE,
        log_time = FALSE,
        plot_exclusions_file = plot_exclusions_file
    )

    save_case_simstudy_plots(
        plots = plots,
        case_name = case_name,
        output_dir = figure_dir,
        formats = formats
    )
}

project_root <- find_project_root()

devtools::load_all(project_root, quiet = TRUE)

case_names <- c("2re", "3re", "sitar")

invisible(
    lapply(
        case_names,
        reproduce_case_figures,
        project_root = project_root,
        formats = c("tikz")
    )
)

message(
    "Done. Regenerated figure files are in ",
    file.path("instances", "<case>", "analysis", "output", "figures"),
    "."
)
