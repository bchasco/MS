myvars <- c("ReleaseSite","loc","state","Species","Year", "travel.days","adult.days")
x <- read.csv("data/combinedData.csv") %>%
  mutate(tag = 1,
         Length = round(Length/5)*5,
         travels.days = round(as.numeric(as.character(travel.days))/7),
         adult.days = round(as.numeric(as.character(adult.days))/7)) %>%
  pivot_longer(cols = c(tag,as.Smolt,As.Adult.ballard), names_to = "loc", values_to = "state") %>%
  group_by(PITCode) %>%
  summarise(
    pattern = paste(apply(across(all_of(myvars)), 1, paste, collapse = "\t"), collapse = "# ")
  ) %>%  # Now group by the unique cum_time patterns
  group_by(pattern) %>%
  # # Summarize or count the Tag Codes that have the same cum_time pattern
  summarise(tag_codes = list(PITCode),
            count = n(), # Number of Tag Codes with this cum_time pattern
            .groups = 'drop')


ch_list <- x %>%
  rowwise() %>%
  # Split cum_time_pattern by commas, then by spaces, and convert values to numeric
  mutate(matrix = list(
    # Split and then bind rows, ensuring we keep them in a matrix format
    matrix(do.call(rbind, strsplit(strsplit(pattern, "# ")[[1]], "\t"))) %>%
      apply(2, as.character)  # Apply as.numeric to each column
    # matrix(strsplit(strsplit(cum_time_pattern, ", ")[[1]], " "),) #%>%
    # apply(2, as.character)  # Apply as.numeric to each column
  )) %>%
  # Extract the first column as cum_time and the remaining columns as covariates (covars)
  mutate(#cum_time = list(matrix[, 1]),         # First column is cum_time
    X = list(matrix(matrix[, ]))) %>%      # Remaining columns are covariates (covars)
  # Create a list where the first element is count, the second is cum_time, and the third is covars
  mutate(list_element = list(list(n = count, X = X))) %>%
  # Extract the new list column
  pull(list_element)


# Apply the transformation to each element in the list, adding an 'id' column
ch_list <- lapply(seq_along(ch_list), function(i) {
  # Access the current element
  x <- ch_list[[i]]

  # Create a matrix from X
  m <- matrix(x$X, ncol = length(myvars), byrow = FALSE)

  # Convert the matrix to a data frame
  df <- as.data.frame(m, stringsAsFactors = FALSE)

  # Assign the column names to the data frame
  names(df) <- myvars

  # Add 'n' as a column in the data frame
  df$n <- x$n

  # Add the index as 'id' in the data frame
  df$id <- i

  # Replace X with the new data frame
  x$X <- df

  # Return the modified list element
  return(x)
})

LkWA <- do.call(rbind, lapply(ch_list, function(x) (x$X)))
LkWA[LkWA=='NA'] <- NA
LkWA <- LkWA %>%
  mutate(time = ifelse(loc == "tag", 0, ifelse(loc == "as.Smolt", travel.days, adult.days))) %>%
  mutate(time = na_if(time, "#N/A"))# %>%
# filter(!is.na(days))

save(LkWA, file="data/LkWA_time.rda")
