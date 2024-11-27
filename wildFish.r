library(dplyr)
library(tidyr)
library(ggplot2)
rm(list = ls())
load("data/wildFish.rda")
WENh <- WENh %>%
  group_by(id) %>%
  mutate(init_state = first(stage)) %>%
  mutate( stage = ifelse(init_state == "yr" & stage == "sub", "yr", stage))

# load("data/NASON_allFish.rda")
# WENh$stage[!is.na(WENh$stage)] <- "yr"


locs <- c("Trib","WEN","BON","TWX_EST")
last_site <- "TWX_EST"
states <- c("sub","yr", "unk")

WENh <- WENh %>%
  # mutate(stage = ifelse("sub","yr",stage)) #%>%
  mutate(state = ifelse(is.na(stage),"unk",stage)) %>%
  mutate(state = factor(state, levels = states)) %>%
  droplevels() %>%
  mutate(loc = factor(loc, levels = c("Trib","First_Trap","RIS_RIA", "WEN","MCJ","JDJ","BON","TWX_EST"))) %>%
  mutate(tag_site = "Trib",
         last_site = locs[length(locs)]) %>%
  mutate(locyr = case_when(
    loc == "BON" ~ "WEN",
    TRUE ~ as.character(loc)  # Convert `loc` to character here
  )) %>%
  mutate(loceta = case_when(
    loc == "TWX_EST" ~ "BON",
    TRUE ~ as.character(loc)  # Convert `loc` to character here
  )) %>%
  mutate(locyr = factor(locyr)) %>%
  filter(loc %in% locs) %>%
  droplevels() %>%
  droplevels()



names(WENh)[3] <- "ReleaseSite"

# WENh <- WENh %>%
  # filter(ReleaseSite == "Updated2023 ChiwComHist")
  # filter(ReleaseSite = "Chinook_LwrWENRST_W") %>%
  # mutate(ReleaseSite = "Chiwawa_Nason")

MR_settings <- list(state = 'LH',
                    frms = list(
                                phi = list(sub = ' ~ -1 + loc:ReleaseSite'
                                           ,yr = ' ~ -1 + loc'), #survival
                                p = list(sub = ' ~ -1 + loc',
                                         yr = ' ~ -1 + loc'), #detection
                                lam = ' ~ 1', #nuisance parameter
                                eta = ' ~ -1 + locyr'), #transition probability
                    mod = "MS")

# Create an instance of the tmb_list class
input <- new("tmb_list",
                   MR_settings = MR_settings,
                   data = WENh)


#take the input and put it into a tmb model
fit <- MStmb(input)
# Assuming `fit@TMB$opt$par` is your named vector
params <- fit@TMB$opt$par
# Split the vector into a list by names
grouped_params <- split(params, names(params))

sd <- RTMB::sdreport(fit@TMB$obj)
# Check the result
print(grouped_params)
source("ggplot_pred_p.r")
source("ggplot_pred_phi.r")
g <- ggpubr::ggarrange(phi, p, nrow = 2, align = "v")
# print(g)
# source("ggplot_cum_pred_x.r")
