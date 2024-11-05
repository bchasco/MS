library(tidyr)
library(dplyr)

# rm(list=ls())


# #Build simple data
locs <- c('CHImark','CHI','WEN','MCJ','JDJ','BON')
states <- c("Sub", "Yr", "Dead")

#Read in all of the sheets in Tom's data set
sheetNames <- c("2021 ChiwCompHist"
                ,"2022 ChiwCompHist"
                ,"Updated2023 ChiwComHist" #had to change a couple of the header names, Event Site Code Value	Event Capture Method Name
                ,"2021 NasCompHist"
                ,"2022 NasCompHist"
                ,"2023 NasCompHist"
                ,"NasRST2021SubsCompHist"
                ,"NasRST2022YCompHist"
                ,"NasRST2022SubsCompHist"
                ,"NasRST2023YCompHist"
                ,"NasRST2023SubsCompHist"
                ,"NasRST2023YCompHist"
                ,"ChiwRST2021SubsCompHist"
                ,"ChiwRST2022YCompHist"
                ,"ChiwRST2022SubsCompHist"
                ,"ChiwRST2023YCompHist"
                ,"ChiwRST2023SubsCompHist"
                ,"ChiwRST2023YCompHist"
                ,"LWenRST2022YCompHist"
                ,"LWenRST2023YCompHist"
                ,"LWenRST2024YCompHist"
)

sheets <- list()
for(i in seq_along(sheetNames)){
  sheets[[i]] <- readxl::read_xlsx('data/Upper Wenatchee 11W Remote Tagging_AllTags.xlsx', sheet = sheetNames[i])
  names(sheets[[i]]) <- names(sheets[[1]])
}

#Combine all of the datasheet
dat <-
  do.call('rbind',sheets)

#Compress all of the individual locations into a set of nine locations based on the following lookup table
recode <- read.csv("data/location.csv")

#These are locations to use in the analysis
locs <- c(1,2,4,6,7,8,9) #Location to include, Leaveout Tumawter and RIS/RIA
dimnames <- data.frame("Caploc" = c("Trib","First_Trap","Tumwater","Wen","RIS_RIA","MCJ","JDJ","BON","TWX_EST"),
                       "dimID" = c(1,2,3,4,5,6,7,8,9))

#This is a flat file with a line for each mark or recap of a fish
#Transform a lot of the data to get dates and stages
ch <- dat[,] %>%
  # slice(-grep("LWB",dat$`Event Site Code Value`)) %>% #remove barge data
  group_by(`Tag Code`) %>%
  mutate(init_year = min(`Event Year YYYY`)) %>%
  mutate(event_year = `Event Year YYYY`) %>%
  filter((event_year - init_year)<=1) %>% #remove adults
  mutate(event_day = lubridate::yday(`Event Date MMDDYYYY`)) %>%
  mutate(stage = ifelse(event_day<=181,"yr","sub")) %>% #<July 1 is a yearling
  group_by(`Tag Code`,init_year,`Event Site Code Value`) %>%
  mutate(max_event_day = max(event_day)) %>%  #find the max detection for a site if there are dups
  filter(event_day == max_event_day) %>%
  arrange(`Tag Code`,event_year,event_day) %>%
  distinct(`Tag Code`,`Event Year YYYY`,`Event Site Code Value`, event_day, .keep_all = TRUE) %>%
  group_by(`Tag Code`) %>%
  mutate(init_stage = first(stage), init_loc = first(`Event Site Code Value`), init_date = first(`Event Date MMDDYYYY`)) %>%
  mutate(stage = ifelse(init_stage == "yr", "yr", stage)) %>%
  mutate(cum_time = julian(as.Date(`Event Date MMDDYYYY`), origin = as.Date("1970-01-01")) - julian(as.Date(init_date), origin = as.Date("1970-01-01"))) %>%
  filter(cum_time <= 365)

#Group the different locations into coarser geographic groups based on the following look up
ch$loc_code <- recode$Code[match(ch$`Event Site Code Value`,recode$Location)]
ch$IDnum <- recode$Idnum[match(ch$`Event Site Code Value`,recode$Location)]


#Remove recaptures not represented by the location lookup above.
myvars <- c('loc','tagSite', 'stage')
ch_complete <- ch[,] %>%
  group_by(`Tag Code`) %>%
  mutate(init_code = first(loc_code), init_num = first(IDnum)) %>%
  arrange(init_num) %>%
  mutate(tagSite = first(`Event Site Code Value`)) %>%
  filter(IDnum!="")%>%
  mutate(stageID = ifelse(stage == "sub", 1, 2)) %>%
  filter(init_num%in%(locs)) %>% #be suer you don't take capture histories for fish tagged upstream
  select(`Tag Code`, stage, stageID,IDnum, loc_code, tagSite) %>%
  left_join(dimnames, by = c("IDnum" = "dimID")) %>% #This adds the capture locations from dimnames
  mutate(loc = ifelse(!is.na(Caploc), Caploc, IDnum)) %>%
  filter(loc%in%dimnames$Caploc[dimnames$dimID%in%locs]) %>%
  mutate(loc = factor(loc, levels = dimnames$Caploc[locs])) %>%
  group_by(`Tag Code`, loc, tagSite) %>%
  summarise(stage = first(stage))  %>% #REmove duplicates and return the stage when a fish was first detected at a location
  group_by(`Tag Code`) %>%
  # Expand the data to include all locations for each fish
  complete(loc = dimnames$Caploc[locs]) %>%
  mutate(loc = factor(loc, levels = dimnames$Caploc[locs])) %>% #You have to reorder the locations, AGAIN!!
  # Remove duplicates by keeping only the max cum_time for each loc
  group_by(`Tag Code`, loc) %>%
  slice(1) %>%                 # Keep the first row (max cum_time or NA)
  # Arrange the final output
  arrange(`Tag Code`, loc) %>%
  group_by(`Tag Code`) %>%
  # Collapse the cum_time values across locations into a vector or string
  # summarise(cum_time_pattern = paste(loc, cum_time,init_year,event_day, stageID, collapse = ", ")) %>%
  summarise(
    cum_time_pattern = paste(apply(across(all_of(myvars)), 1, paste, collapse = "\t"), collapse = "# ")
  ) %>%  # Now group by the unique cum_time patterns
  group_by(cum_time_pattern) %>%
  # Summarize or count the Tag Codes that have the same cum_time pattern
  summarise(tag_codes = list(`Tag Code`),
            count = n(), # Number of Tag Codes with this cum_time pattern
            .groups = 'drop')



ch_list <- ch_complete %>%
  rowwise() %>%
  # Split cum_time_pattern by commas, then by spaces, and convert values to numeric
  mutate(matrix = list(
    # Split and then bind rows, ensuring we keep them in a matrix format
    matrix(do.call(rbind, strsplit(strsplit(cum_time_pattern, "# ")[[1]], "\t"))) %>%
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

# Assuming ch_list is your original list and myCols is defined
myCols <- c('loc','tagSite', 'stageID')

# Apply the transformation to each element in the list, adding an 'id' column
ch_list <- lapply(seq_along(ch_list), function(i) {
  # Access the current element
  x <- ch_list[[i]]

  # Create a matrix from X
  m <- matrix(x$X, ncol = length(myCols), byrow = FALSE)

  # Convert the matrix to a data frame
  df <- as.data.frame(m, stringsAsFactors = FALSE)

  # Assign the column names to the data frame
  names(df) <- myCols

  # Add 'n' as a column in the data frame
  df$n <- x$n

  # Add the index as 'id' in the data frame
  df$id <- i

  # Replace X with the new data frame
  x$X <- df

  # Return the modified list element
  return(x)
})

WEN <- do.call(rbind, lapply(ch_list, function(x) (x$X)))
WEN[WEN=='NA'] <- NA
WEN <- WEN %>%
  group_by(id) %>%
  mutate(tagSite = first(na.omit(tagSite)))


save(WEN, file="data/WEN.rda")
