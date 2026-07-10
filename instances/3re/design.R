library(tidyverse)
devtools::load_all()

case_name <- "3re"

subcases <- tidyr::crossing(
    n_clusters = c(50, 100, 200, 300, 400, 500),
    n_obs_per_cluster = c(2, 3, 5, 10)
) %>%
    dplyr::mutate(
        subcase_id = dplyr::row_number()
    )

methods <- final_methods(case_name)

applicable_method_cells <- tidyr::crossing(
    case = case_name,
    subcase_id = subcases$subcase_id,
    method = methods
) %>%
    dplyr::left_join(
        subcases,
        by = "subcase_id"
    ) %>%
    dplyr::filter(
        !(
            method == "Local-GAM" &
                n_obs_per_cluster == 2
        )
    ) %>%
    dplyr::select(
        case,
        subcase_id,
        method
    )
