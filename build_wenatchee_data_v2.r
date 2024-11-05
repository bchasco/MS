library(tidyr)
library(dplyr)

# Load and clean data
sheetNames <- c("2021 ChiwCompHist", "2022 ChiwCompHist", "Updated2023 ChiwComHist", "2021 NasCompHist", "2022 NasCompHist", "2023 NasCompHist")
locs <- c(1, 2, 4, 6, 7, 8, 9)  # Locations to include
dimnames <- data.frame(Caploc = c("Trib","First_Trap","Tumwater","Wen","RIS_RIA","MCJ","JDJ","BON","TWX_EST"),
                       dimID = c(1, 2, 3, 4, 5, 6, 7, 8, 9))
recode <- read.csv("data/location.csv")

# Read and combine sheets
sheets <- lapply(sheetNames, function(sheet) {
  data <- readxl::read_xlsx('data/Upper Wenatchee 11W Remote Tagging_AllTags.xlsx', sheet = sheet)
  names(data) <- names(sheets[[1]])
  data
})
dat <- bind_rows(sheets)

# Data transformation
ch <- dat %>%
  group_by(`Tag Code`) %>%
  mutate(init_year = min(`Event Year YYYY`),
         event_day = lubridate::yday(`Event Date MMDDYYYY`),
         stage = ifelse(event_day <= 181, "yr", "sub")) %>%
  filter((`Event Year YYYY` - init_year) <= 1) %>%  # Exclude adults
  distinct(`Tag Code`, `Event Site Code Value`, event_day, .keep_all = TRUE) %>%
  mutate(loc_code = recode$Code[match(`Event Site Code Value`, recode$Location)],
         IDnum = recode$Idnum[match(`Event Site Code Value`, recode$Location)],
         cum_time = as.numeric(as.Date(`Event Date MMDDYYYY`) - as.Date(first(`Event Date MMDDYYYY`))))

# Filtering and preparing capture histories
ch_complete <- ch %>%
  filter(IDnum %in% locs) %>%
  left_join(dimnames, by = c("IDnum" = "dimID")) %>%
  mutate(loc = ifelse(!is.na(Caploc), Caploc, IDnum),
         loc = factor(loc, levels = dimnames$Caploc[locs])) %>%
  group_by(`Tag Code`, loc) %>%
  arrange(desc(cum_time)) %>%
  slice(1) %>%
  arrange(`Tag Code`, loc)

# Create summary patterns
myvars <- c('loc','cum_time','init_year','event_day', 'stage')

library(tidyr)
library(dplyr)

# Assuming 'ch_complete' has been created as in the earlier steps

# Group and summarize the data based on the unique cum_time patterns
ch_complete_summary <- ch_complete %>%
  group_by(`Tag Code`, loc, init_year, cum_time, event_day) %>%
  summarise(stage = first(stage)) %>% # Remove duplicates and return the stage
  group_by(`Tag Code`) %>%
  complete(loc = dimnames$Caploc[locs]) %>% # Expand the data to include all locations
  mutate(loc = factor(loc, levels = dimnames$Caploc[locs])) %>%
  group_by(`Tag Code`, loc) %>%
  arrange(desc(cum_time)) %>%
  slice(1) %>% # Keep the first row (max cum_time or NA)
  arrange(`Tag Code`, loc) %>%
  group_by(`Tag Code`) %>%
  summarise(
    cum_time_pattern = paste(apply(across(all_of(myvars)), 1, paste, collapse = "\t"), collapse = ", ")
  ) %>%
  group_by(cum_time_pattern) %>%
  summarise(
    tag_codes = list(`Tag Code`),
    count = n(), # Correctly count the number of unique Tag Codes with the same cum_time_pattern
    .groups = 'drop'
  )

# Create the list structure with matrices and the correct 'n' values
ch_list2 <- ch_complete_summary %>%
  rowwise() %>%
  mutate(matrix = list(
    matrix(do.call(rbind, strsplit(strsplit(cum_time_pattern, ", ")[[1]], "\t"))) %>%
      apply(2, as.character)  # Ensure matrix values are in the correct format
  )) %>%
  mutate(X = list(matrix(matrix[, ]))) %>%
  mutate(list_element = list(list(n = count, X = X))) %>%
  pull(list_element)

# Assign 'myCols' for the data frame
myCols <- c('loc', 'cum_time', 'init_year', 'event_day', 'stageID')

# Apply the transformation to each element in 'ch_list', adding 'id' and correct 'n'
ch_list2 <- lapply(seq_along(ch_list), function(i) {
  x <- ch_list[[i]]
  m <- matrix(x$X, ncol = length(myCols), byrow = FALSE)
  df <- as.data.frame(m, stringsAsFactors = FALSE)
  names(df) <- myCols
  df$n <- x$n  # Correctly assign 'n' as the count of unique Tag Codes
  df$id <- i
  x$X <- df
  return(x)
})

# Combine all the data frames
combined_df <- do.call(rbind, lapply(ch_list, function(x) x$X))

# Inspect the result
head(combined_df)
