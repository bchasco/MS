create_poly_basis <- function(object){
  poly_transform <- list()
  for(i in names(object@MR_settings$frms)){
    poly_transform[[i]] <- list()
    state_frms <- object@MR_settings$frms[[i]]
    for(s in names(state_frms)){
      terms <- attr(terms(formula(state_frms[[s]])),"term.labels")
      for(trm in terms){
        if(grepl("poly", trm)){
          var <- gsub("^poly\\(([^,]+),.*$", "\\1", trm)
          df <- gsub(".*,(.*)\\)", "\\1", trm)  # Regex to extract
          df <- as.numeric(trimws(df))                          # Remove any whitespace
          poly_transform[[i]][[s]][[var]] <- function(new_variable) predict(poly(object@data[[var]],df), new_variable)
        }
      }
    }
  }
  return(poly_transform)
}

extract_poly_vars <- function(poly_term) {
  # Use regex to extract the variable inside poly()
  if(attr(gregexpr("poly",poly_term)[[1]],"match.length")>(-1)){
    var_name <- gsub("^poly\\(([^,]+),.*$", "\\1", poly_term)
    return(var_name)
  }
}

extract_unique_vars <- function(frm_list) {
  # Collapse all formulas into a single string
  all_vars <- paste(unlist(frm_list), collapse = "+")

  # Remove unwanted parts like "~", "1", "-1"
  all_vars <- gsub("~|\\b1\\b|-1", "", all_vars)

  # Split by any of "+", "|", "/", or ":" using regex
  terms <- unlist(strsplit(all_vars, "\\+|\\||\\/|:"))

  # Trim whitespace from each term
  terms <- trimws(terms)

  # Process each term to handle polynomial terms and clean variables
  clean_terms <- unlist(lapply(terms, function(term) {
    if (grepl("poly\\(", term)) {
      # Extract variable inside poly() using regex
      return(gsub(".*poly\\(([^,]+).*", "\\1", term))
    }
    # Return the term as is, stripping unwanted parentheses or symbols
    return(gsub("\\(|\\)|\\|", "", term))
  }))

  # Remove duplicates and empty entries
  unique_terms <- unique(clean_terms)
  unique_terms <- unique_terms[unique_terms != ""]

  # Return the cleaned vector
  return(unique_terms)
}

expand_grid <- function(data, cols, include_time = TRUE, time.length = 20, state_included = c("yr","unk")) {
  # Check the input columns
  if (!all(cols %in% colnames(data))) {
    stop("The specified columns are not in the data.")
  }

  data <- data %>%
    filter(state %in% state_included)

  # Extract unique combinations of fY and ReleaseSite
  unique_combinations <- data %>%
    distinct(across(all_of(cols))) %>%
    arrange(across(all_of(cols)))

  if(include_time){
    # q <- quantile(na.omit(tmb.data$data$time),c(0.025,0.975))
    # time_range <- seq(max(0.5,q[1]),
    #                   q[2],
    #                   length.out = time.length)
    time_range <- unique(quantile(na.omit(data$time[data$loc=="as.Smolt"]), probs = seq(0.05,0.95,0.1)))
    # print(time_range)
  }

  # Extract unique locations (loc)
  unique_locs <- unique(data$loc)

  # Expand the grid for each combination
  expanded_data <- unique_combinations %>%
    rowwise() %>%
    mutate(expanded = list(
      expand.grid(
        loc = unique_locs,
        state = state_included,
        time = if (include_time) time_range else NA
      )
    )) %>%
    unnest(expanded)

  add_group_id <- function(data, dynamic_cols) {
    # Ensure dynamic_cols exist in the data
    if (!all(dynamic_cols %in% colnames(data))) {
      stop("Some columns in 'dynamic_cols' are not in the data.")
    }

    # Define grouping columns, always including 'state' and 'time'
    grouping_cols <- c(dynamic_cols, "state", "time")

    # Group by the specified columns and add a unique group ID
    data <- data %>%
      group_by(across(all_of(grouping_cols))) %>%
      mutate(id = cur_group_id()) %>%
      ungroup()  # Ungroup after assigning IDs for safe downstream operations

    return(data)
  }

  # Example usage
  expanded_data_with_ids <- add_group_id(expanded_data, dynamic_cols = cols) %>%
    arrange(id,state)


  return(expanded_data_with_ids)
}
