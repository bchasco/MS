library(dplyr)
library(tidyr)
library(ggplot2)
# rm(list = ls())
# load("data/NASON_allFish_size.rda")
load("data/NASON_allFish.rda")
locs <- c("Trib", "WEN","BON","TWX_EST")
last_site <- "TWX_EST"

WENh <- WENh %>%
  mutate(state = ifelse(is.na(stage),"unk","yr")) %>%
  mutate(loc = factor(loc, levels = c("Trib","First_Trap","RIS_RIA", "WEN","MCJ","JDJ","BON","TWX_EST"))) %>%
  mutate(tag_site = "Trib",
         last_site = locs[length(locs)]) %>%
  filter(loc %in% locs) #%>%
    # mutate(len = as.numeric(ind_len)/100) %>%
    # mutate(`sp` = ifelse(grepl("Chinook",`Release Site`),"Chinook", ifelse(grepl("Coho",`Release Site`),"Coho", "Steelhead"))) %>%
    # filter(sp != "Coho")# %>%
# mutate(`Release Site` = ifelse(grepl("Chinook",`Release Site`),ifelse(grepl("HXH",`Release Site`),"Chinook" ,"Chinook" ),
#                                ifelse(grepl("Steelhead",`Release Site`),ifelse(grepl("HXH",`Release Site`),"Steelhead" ,"Steelhead" ), "Coho")))


names(WENh)[3] <- "ReleaseSite"



# Define MR_settings
# Get unique values for each of the first three columns
unique_locs <- unique(WENh$loc)
unique_releases <- unique(WENh$ReleaseSite)

# Use expand.grid to create all combinations
expanded_grid <- data.frame(expand.grid(
  # len = mean(WENh$len),
  loc = unique_locs,
  # sp = unique(WENh$sp),
  current_state = unique(WENh$state),
  next_state = unique(WENh$state),
  ReleaseSite = unique_releases,
  init_site = "Trib",
  last_iste = last_site
))
expanded_grid <- expanded_grid %>%
  rowwise() %>%
  # mutate(same_sp = grepl(sp, ReleaseSite)) %>%
  # filter(same_sp == TRUE) %>%
  group_by(across(-contains("state"))) %>%
  mutate(id = cur_group_id()) %>%
  ungroup() %>%
  arrange(id) %>%
  mutate(tag_site = "Trib", last_site = locs[length(locs)]) %>%
  filter(loc != tag_site & loc != last_site) %>%
  droplevels()


# WEN_data <- WENh %>%
#   mutate(state = ifelse(is.na(stage), "unk", "yr")) %>%
#   mutate(state = factor(state, levels = c("yr","unk"))) %>%
#   mutate(loc = factor(loc, levels = c("Trib","First_Trap","RIS_RIA", "WEN","MCJ","JDJ","BON","TWX_EST"))) %>%
#   mutate(tag_site = "Trib",
#          last_site = locs[length(locs)]) %>%
#   group_by(state) %>%
#   droplevels()
#

MR_settings <- list(state = 'LH',
                    frms = list(phi = ' ~ -1 + loc:ReleaseSite ', #survival
                                p = ' ~ -1 + loc', #detection
                                lam = ' ~ 1', #nuisance parameter
                                eta = ' ~ 1'), #transition probability
                    pdm = list(),
                    pred_data = expanded_grid,
                    mod = "MS")

# Create an instance of the tmb_list class
input <- new("tmb_list",
                   MR_settings = MR_settings,
                   data = WENh)


#take the input and put it into a tmb model
fit <- MStmb(input)
sd <- RTMB::sdreport(fit@TMB$obj)

# round((fit@TMB$rep$gamma_pred[,,1,]),3)
#
# WENh %>%
#   group_by(ReleaseSite, loc) %>%
#   mutate(loc = factor(loc, levels = locs)) %>%
#   arrange(loc) %>%
#   na.omit() %>%
#   summarize(n = sum(n)) %>%
#   pivot_wider(names_from = ReleaseSite, values_from = n)
#
# if(MR_settings$mod=="MS"){
#
#   print(round(plogis(fit@TMB$rep$gamma_pred[1,1,,]),3))
#   print(round(plogis(fit@TMB$rep$omega[1,1,]),3))
# }else{
#   print(round((plogis(fit@TMB$opt$par)),3))
#   print(round((fit@TMB$rep$omega[1,1,]),3))
# }
# print(sd)
#
source("ggplot_survival_plot.r")

x <- cbind(fit@MR_settings$pred_data , round(fit@TMB$rep$pred_x,3), round(fit@TMB$rep$cum_pred_x,3))
x[x$ReleaseSite==levels(x$ReleaseSite)[4] & x$current_state == x$next_state,]
