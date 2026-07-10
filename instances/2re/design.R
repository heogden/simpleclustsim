library(tidyverse)
devtools::load_all()

case_name <- "2re"

subcases <- crossing(
    n_clusters = c(50, 100, 200, 300, 400, 500),
    n_obs_per_cluster = c(2, 3, 5, 10)
) %>%
    mutate(
        subcase_id = dplyr::row_number()
    )

methods <- final_methods(case_name)

applicable_method_cells <- crossing(
    case = case_name,
    subcase_id = subcases$subcase_id,
    method = methods
) %>%
    left_join(
        subcases,
        by = "subcase_id"
    ) %>%
    filter(
        !(method == "Local-GAM" &
              n_obs_per_cluster == 2)
    ) %>%
    select(
        case,
        subcase_id,
        method
    )
