# create_new_dm_grid <- function(var, model_fit){
#
#   output <- list()
#
#   for(s in seq_along(model_fit@TMB$obj$tmb.data$dm[[var]])){
#     df <- model_fit@TMB$obj$tmb.data$dm[[var]][[s]]$X
#     fr <- model_fit@TMB$obj$tmb.data$dm[[var]][[s]]$fr
#     frm <- model_fit@MR_settings$frms[[var]][[s]]
#
#     if(attr(gregexpr("\\|", frm)[[1]], 'match.length') != (-1)){
#
#       trms <- all.vars(formula(frm))
#       factor_terms <- sapply(trms, function(var) !is.numeric(fr))
#       continuous_terms <- sapply(trms, function(var) is.numeric(fr))
#
#       reTrms <- names(model_fit@TMB$obj$tmb.data$dm[[var]][[s]]$reTrms$nl)
#       factor_terms <- names(factor_terms[!(factor_terms%in%c(reTrms,continuous_terms))])
#       factor_id <- which(factor_terms%in%trms, arr.ind = TRUE)
#       factor_cols <- which(attr(df, "assign") %in% factor_id,
#                            arr.ind = TRUE)
#
#       continuous_id <- which(continuous_terms%in%trms, arr.ind = TRUE)
#       continuous_cols <- which(attr(df, "assign") %in% continuous_id,
#                                arr.ind = TRUE)
#
#       reTrms <- t(as.matrix(model_fit@TMB$obj$tmb.data$dm[[var]][[1]]$reTrms$Zt))
#
#
#     }else{
#
#       trms <- attr(terms(formula(model_fit@MR_settings$frms[[var]][[s]])),"term.labels")
#       has_intercept <- attr(terms(formula(frm)), "intercept") == 1
#       factor_terms <- names(attr(df, "contrasts"))
#       if(has_intercept){
#         continuous_terms <- c("(Intercept)",trms[!(trms%in%factor_terms)])
#       }else{
#         continuous_terms <- trms[!(trms%in%factor_terms)]
#       }
#
#       factor_id <- which(trms%in%factor_terms, arr.ind = TRUE)
#       continuous_id <- which(trms%in%continuous_terms, arr.ind = TRUE)
#
#       factor_cols <- which(attr(df, "assign") %in% factor_id,
#                            arr.ind = TRUE)
#       continuous_cols <- which(attr(df, "assign") %in% continuous_id,
#                                arr.ind = TRUE)
#
#       poly_basis <- list()
#
#       for(i in seq_along(continuous_terms)){
#         tmp_var <- extract_poly_vars(continuous_terms[i])
#         if(attr(gregexpr("poly",continuous_terms[i])[[1]], "useBytes") & continuous_terms[i]!="(Intercept)"){
#           poly_basis[[tmp_var]]$poly <- create_poly_basis(fr[,tmp_var],2)
#           poly_basis[[tmp_var]]$range <- seq(min(fr[,tmp_var]),max(fr[,tmp_var]),length.out = 10)
#           poly_basis[[tmp_var]]$new_transform <- poly_basis[[tmp_var]]$poly$transform(poly_basis[[tmp_var]]$range)
#           continuous_terms[i] <- tmp_var
#         }
#       }
#     }
#
#     # Get the row numbers of distinct rows
#     fr_tmp <- as.data.frame(fr)
#     distinct_row_numbers <- fr_tmp %>%
#       distinct(across(all_of(c(factor_terms))),.keep_all = TRUE) %>%  # Keep all columns
#       pull(id)# Extract the row numbers
#     fr_tmp <- fr_tmp %>%
#       filter(id %in% distinct_row_numbers) %>%
#       select(factor_terms)
#
#     if(has_intercept){
#       fr_tmp <- data.frame("(Intercept)" = 1, fr_tmp)
#     }
#
#     if(sum(continuous_terms!="(Intercept)")>0){
#       ranges_df <- extract_ranges_to_df(poly_basis)
#       # Your initial data frames
#       fr_tmp_expanded <- expand_grid(data.frame(fr_tmp), loc = unique(fit@data$loc))
#
#       # Expand by each column in ranges_df
#       new_grid <- fr_tmp_expanded
#       for (col in colnames(ranges_df)) {
#         new_grid <- expand_grid(new_grid, !!sym(col) := ranges_df[[col]])
#       }
#
#     }else{
#       new_grid <- expand_grid(data.frame(fr_tmp), loc = unique(fit@data$loc))
#     }
#     new_dm <- model.matrix(formula(frm), new_grid)
#
#     for(i in seq_along(continuous_terms[!(continuous_terms%in%"(Intercept)")])) {
#       print(i)
#       new_dm[,attr(df, "assign") %in% continuous_id[i]] <- poly_basis[[continuous_terms[i]]]$poly$transform(as.vector(new_grid[,c(continuous_terms[i])][[1]]))
#     }
#
#     #parameter estimates
#     pars <- unlist(c(model_fit@TMB$rep[[paste0(var,"_s")]],model_fit@TMB$rep[[paste0(var,'_re_s')]]))
#
#     est <- new_dm %*% pars
#     est <- exp(est)/(exp(1) + exp(est))
#     new_grid$est <- est
#   }
#
#   return(list(new_grid = new_grid, new_dm = new_dm))
# }
#
# extract_poly_vars <- function(poly_term) {
#   # Use regex to extract the variable inside poly()
#   var_name <- gsub("^poly\\(([^,]+),.*$", "\\1", poly_term)
#   return(var_name)
# }
#
# create_poly_basis <- function(variable, degree) {
#   # Store the original basis
#   basis <- poly(variable, degree = degree)
#   list(
#     basis = basis,
#     transform = function(new_variable) predict(basis, new_variable)
#   )
# }
#
# extract_ranges_to_df <- function(poly_basis) {
#   # Extract the names of the sublists
#   sublist_names <- names(poly_basis)
#
#   # Extract the range from each sublist
#   ranges <- lapply(sublist_names, function(name) {
#     if (!is.null(poly_basis[[name]]$range)) {
#       poly_basis[[name]]$range
#     } else {
#       # Handle missing ranges (e.g., return NA if 'range' is absent)
#       NA
#     }
#   })
#
#   # Combine the ranges into a data frame
#   ranges_df <- as.data.frame(ranges)
#
#   # Set column names to the names of the sublists
#   colnames(ranges_df) <- sublist_names
#
#   return(ranges_df)
# }
#
# create_poly_basis_list <- function(continuous_terms, fr) {
#   poly_basis <- list()
#
#   for (term in continuous_terms) {
#     if (grepl("poly", term) && term != "(Intercept)") {
#       var_name <- extract_poly_vars(term)
#       poly_basis[[var_name]] <- list(
#         poly = create_poly_basis(fr[, var_name], degree = 2),
#         range = seq(min(fr[, var_name]), max(fr[, var_name]), length.out = 10)
#       )
#       poly_basis[[var_name]]$new_transform <- poly_basis[[var_name]]$poly$transform(poly_basis[[var_name]]$range)
#     }
#   }
#
#   return(poly_basis)
# }
#
# create_new_grid <- function(fr_tmp, loc_values, poly_basis, continuous_terms) {
#   if (length(continuous_terms[continuous_terms != "(Intercept)"]) > 0) {
#     ranges_df <- extract_ranges_to_df(poly_basis)
#     expanded_grid <- expand_grid(fr_tmp, loc = loc_values)
#
#     for (col in colnames(ranges_df)) {
#       expanded_grid <- expand_grid(expanded_grid, !!sym(col) := ranges_df[[col]])
#     }
#
#     return(expanded_grid)
#   } else {
#     return(expand_grid(fr_tmp, loc = loc_values))
#   }
# }
#
# test <- create_new_dm_grid(var = "lam",fit)
#
# test$new_grid %>%
#   filter(loc == "as.Smolt") %>%
#   # filter(Year == 2002) %>%
#   ggplot(aes(x = fY, y = est)) +
#   # facet_wrap(~ReleaseSite) +
#   geom_point()
