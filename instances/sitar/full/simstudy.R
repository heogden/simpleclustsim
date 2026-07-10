## Change `stage` as the Iridis workflow progresses:
##   "prepare" -> "after_4h" -> "after_8h" ->
##   "after_recovery" -> "finalise"
##
## Some stages may be skipped when the preceding stage reports that no
## update or recovery run is required.

project_root <- here::here()
devtools::load_all(project_root)

case_to_run <- "sitar"
stage <- "finalise"

source(
    here::here(
        "instances",
        case_to_run,
        "design.R"
    )
)

if (!identical(case_name, case_to_run)) {
    stop(
        "The sourced design has case_name = ",
        case_name,
        ", but case_to_run = ",
        case_to_run,
        ".",
        call. = FALSE
    )
}

run_full_simstudy_stage(
    case_name = case_name,
    subcases = subcases,
    applicable_method_cells =
        applicable_method_cells,
    simulation_fun = run_simstudy_each,
    simulation_fun_name = "run_simstudy_each",
    stage = stage,
    initial_walltime_minutes = 240,
    extended_walltime_minutes = 480,
    mem_gb = 8,
    binpacking_back = 0,
    pack_by = "method",
    project_root = project_root
)
