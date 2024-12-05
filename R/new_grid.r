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

create_expanded_grid <- function(object, breaks = 10) {
  search_vars <- extract_unique_vars(object@MR_settings$frms)
  var_names <- search_vars
  print(var_names)
  # ii <- 1
  expanded_data <- lapply(search_vars, function(var) {
    #Get rid of the poly for expanding the grid
    if(attr(gregexpr("poly",var)[[1]],"match.length")>(-1)){
      var <- extract_poly_vars(var)
    }

    if(is.character(object@data[[var]]) || (is.factor(object@data[[var]]) && is.character(as.character(object@data[[var]])))){
      unique(object@data[[var]])
    } else if (is.numeric(object@data[[var]])) {
      # Expand using a sequence for continuous variables
      seq(min(object@data[[var]], na.rm = TRUE), max(object@data[[var]], na.rm = TRUE), length.out = breaks)
    } else {
      stop(paste("Variable type for", var, "is not supported."))
    }
    # ii <- ii + 1
  })

  expanded_data[["loc"]] <- unique(object@data$loc)
  expanded_data[["state"]] <- unique(object@data$state)[1]

  for(i in seq_along(var_names)){
    if(attr(gregexpr("poly",var_names[i])[[1]],"match.length")>(-1)){
      var_names[i] <- extract_poly_vars(var_names[i])
    }
  }

  names(expanded_data) <- c(var_names,"loc","state")  # Assign variable names

  tmp_cols <- names(expanded_data)[!(names(expanded_data)%in%c("loc"))]
  # Create the expanded grid
  expanded_grid <- expand.grid(expanded_data)
  expanded_grid <- expanded_grid %>%
    arrange(across(all_of(tmp_cols))) %>%
    group_by(across(all_of(tmp_cols))) %>%
    mutate(id = cur_group_id())

  return(expanded_grid)
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
  all_vars <- paste(unlist(frm_list), collapse = " ")
  all_vars <- gsub("~|1|-1", "", all_vars)

  # Split by "+" to get individual components
  terms <- trimws(unlist(strsplit(all_vars, "\\+")))

  # Process each term to handle polynomial terms and clean variables
  clean_terms <- unlist(lapply(terms, function(term) {
    if (grepl("poly\\(", term)) {
      # Extract variable inside poly() using regex
      gsub(".*poly\\(([^,]+).*", "\\1", term)
    }
    # Return the term as is
    if(length(grep(':',term)>0)){
      strsplit(gsub("\\(|\\)|\\|", "", term), "\\:|\\/")
    }else{
      gsub("\\(|\\)|\\|", "\\1", term)
    }
  }))

  # Remove duplicates and empty entries
  unique_terms <- unique(clean_terms)
  unique_terms <- unique_terms[unique_terms != ""]

  # Return the cleaned vector
  return(unique_terms)
}
