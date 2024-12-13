library(tidyr)
library(dplyr)

rm(list=ls())

# #Build simple data
# locs <- c('CHImark','CHI','WEN','MCJ','JDJ','BON')
states <- c("sub", "yr", "unk")

sheetNames <- c(

  'EntUpperRST2024YRLCompHist'
  ,'EntLowerRST2024YRLCompHist'
  ,"MethLowerRST2024YRLCompHist"
  ,"TwispRST2024YRLCompHist"
  ,"EntLowerRST2023SubsCompHist"
)

sheets <- list()
for(i in seq_along(sheetNames)){
  sheets[[i]] <- readxl::read_xlsx('data/Entiat_Methow_23Subs24Yearlings_11W.xlsx', sheet = sheetNames[i])
  names(sheets[[i]]) <- names(sheets[[1]])
}

combined_data <-
  do.call('rbind',sheets)

#Compress all of the individual locations into a set of nine locations based on the following lookup table
recode <- read.csv("data/hatchery_lookup.csv", header = TRUE)


#parse out the event site information
#This throws an error because one of the entries has 2 "-" instead of 1. Doesn't affect the result.
combined_data$`Event Site Code Value` <- do.call(rbind,strsplit(combined_data$`Event Site Name`," - "))[,1]


#These are locations to use in the analysis
locs <- c(2,3,4,8) #Location to include, Leaveout Tumawter and RIS/RIA
dimnames <- data.frame("Caploc" = c("Trib","First_Trap","Lwr_Main_TRAP","RRJ","RIS_RIA","MCJ","JDJ","BON","TWX_EST"),
                       "dimID" = c(1,2,3,4,5,6,7,8,9))

#This is a flat file with a line for each mark or recap of a fish


#This is a flat file with a line for each mark or recap of a fish
#Transform a lot of the data to get dates and stages
ch <- combined_data[,] %>%
  # slice(-grep("LWB",dat$`Event Site Code Value`)) %>% #remove barge data
  group_by(`Tag Code`) %>%
  mutate(init_year = min(`Event Year YYYY`)) %>%
  mutate(event_year = `Event Year YYYY`) %>%
  filter((event_year - init_year)<=1) %>% #remove adults
  mutate(`Event Date MMDDYYYY` = lubridate::ymd(`Event Date MMDDYYYY`)) %>%
  mutate(event_day = lubridate::yday(`Event Date MMDDYYYY`)) %>%
  mutate(stage = ifelse(event_day<=181,"yr","sub")) %>% #<July 1 is a yearling
  group_by(`Tag Code`,
           init_year,
           `Event Site Code Value`) %>%
  mutate(`Event Type Name` = factor(`Event Type Name`, levels = c("Mark","Observations"))) %>%
  arrange(`Tag Code`,
          event_year,
          event_day) %>%
  distinct(`Tag Code`,
           `Event Year YYYY`,
           `Event Site Code Value`,
           event_day, .keep_all = TRUE)

# Create the lookup table
species_lookup <- ch %>%
  distinct(`Tag Code`, `Event Species Name`) %>%
  group_by(`Tag Code`) %>%
  slice(1) %>%
  ungroup()

#Group the different locations into coarser geographic groups based on the following look up
ch$loc_code <- recode$Code[match(ch$`Event Site Name`,recode$`Event.Site.Name`)]
ch$IDnum <- recode$Idnum[match(ch$`Event Site Name`,recode$`Event.Site.Name`)]
ch$init_species <- species_lookup$`Event Species Name`[match(ch$`Tag Code`,species_lookup$`Tag Code`)]


ch <- ch %>%
  # left_join(species_lookup,by = c("Tag Code" = "Tag Code")) %>%
  # mutate(init_species = `Event Species Name.x`) %>%
  group_by(`Tag Code`) %>%
  mutate(init_stage = first(stage),
         init_loc = first(`Event Site Code Value`),
         init_date = first(`Event Date MMDDYYYY`)) %>%
  mutate(init_week = lubridate::week(init_date)) %>%
  # mutate(stage = ifelse(init_stage == "sub", "yr", stage)) %>%
  mutate(cum_time = julian(as.Date(`Event Date MMDDYYYY`),
                           origin = as.Date("1970-01-01")) - julian(as.Date(init_date), origin = as.Date("1970-01-01"))) %>%
  filter(cum_time <= 365)


#Remove recaptures not represented by the location lookup above.
myvars <- c('loc','stage', 'Release Site', "init_week","cum_time","Event Species Name")

# myvars <- c('loc','tagSite', 'stage')
ch_complete <- ch[,] %>%
  group_by(`Tag Code`) %>%
  mutate(init_code = first(loc_code), init_num = first(IDnum)) %>%
  arrange(init_num) %>%
  mutate(tagSite = first(`Event Site Code Value`)) %>%
  filter(IDnum!="")%>%
  mutate(stage = factor(stage, levels = states)) %>%
  mutate(state = as.integer(stage)) %>%
  droplevels() %>%
  filter(init_num%in%(locs)) %>% #be suer you don't take capture histories for fish tagged upstream
  select(`Tag Code`, stage, IDnum, loc_code, init_week, cum_time, `Event Species Name`) %>%
  left_join(dimnames, by = c("IDnum" = "dimID")) %>% #This adds the capture locations from dimnames
  mutate(loc = ifelse(!is.na(Caploc), Caploc, IDnum)) %>%
  filter(loc %in% dimnames$Caploc[dimnames$dimID%in%locs]) %>%
  mutate(loc = factor(loc, levels = dimnames$Caploc[locs])) %>%
  group_by(`Tag Code`,loc) %>%
  mutate(stage = first(stage)) %>%
  droplevels() %>%
  # filter(loc %in% c("Trib", "First_Trap", "WEN", "MCJ", "JDJ", "BON", "TWX_EST")) %>%
  group_by(`Tag Code`, loc, `Event Species Name`, stage, init_week) %>%
  summarise(cum_time = first(cum_time))  %>% #REmove duplicates and return the stage when a fish was first detected at a location
  # summarise(stage = first(stage))  %>% #REmove duplicates and return the stage when a fish was first detected at a location
  # distinct(`Tag Code`, loc, `Event Species Name`, init_week) %>%
  group_by(`Tag Code`) %>%

  # Expand the data to include all locations for each fish
  complete(loc = dimnames$Caploc[locs]) %>%
  mutate(loc = factor(loc, levels = dimnames$Caploc[locs])) %>% #You have to reorder the locations, AGAIN!!

  # Remove duplicates by keeping only the max cum_time for each loc
  group_by(`Tag Code`, loc, `Event Species Name`, cum_time, init_week) %>%
  arrange(`Tag Code`, loc) %>%
  group_by(`Tag Code`) %>%
  mutate(`Release Site` = first(`Event Species Name`)) %>%
  # Collapse the cum_time values across locations into a vector or string
  # summarise(cum_time_pattern = paste(loc, cum_time,init_year,event_day, stageID, collapse = ", ")) %>%
  summarise(
    pattern = paste(apply(across(all_of(myvars)), 1, paste, collapse = "\t"), collapse = "# ")
  ) %>%  # Now group by the unique cum_time patterns
  group_by(pattern) %>%
  #
  # # Summarize or count the Tag Codes that have the same cum_time pattern
  summarise(tag_codes = list(`Tag Code`),
            count = n(), # Number of Tag Codes with this cum_time pattern
            .groups = 'drop')



ch_list <- ch_complete %>%
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

Ent <- do.call(rbind, lapply(ch_list, function(x) (x$X)))
Ent[Ent=='NA'] <- NA

save(Ent, file="data/Ent.rda")
