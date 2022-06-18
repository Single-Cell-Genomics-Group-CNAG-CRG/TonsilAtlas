gc_dict <- list()

hard <- "{gct}/{robj_dir}/candidates_NFkB_subset_expression_HARD.rds" %>%
  glue::glue() %>%
  here::here() %>%
  readRDS(.)

gc_dict[["HARD"]] <- hard

mid <- "{gct}/{robj_dir}/candidates_NFkB_subset_expression_MID.rds" %>%
  glue::glue() %>%
  here::here() %>%
  readRDS(.)

gc_dict[["MID"]] <- mid

candidates <- "{gct}/{robj_dir}/candidates_NFkB.rds" %>%
    glue::glue() %>%
    here::here() %>%
    readRDS(.)

gc_dict[["candidate-targets"]] <- candidates

all_t <- "{gct}/{robj_dir}/all_targets_NFkB.rds" %>%
    glue::glue() %>%
    here::here() %>%
    readRDS(.)

gc_dict[["all-targets"]] <- all_t

all_inter <- "{gct}/{robj_dir}/all_targets_NFkB_intersection_upregulated_and_closestFeatures.rds" %>%
  glue::glue() %>%
  here::here() %>%
  readRDS(.)

gc_dict[["all-targets-intersect"]] <- all_inter

gc_vec <- unique(unlist(gc_dict))

