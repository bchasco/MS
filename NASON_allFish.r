library(dplyr)
library(tidyr)
library(ggplot2)
rm(list = ls())
load("data/NASON_allFish.rda")
WENh$stage[!is.na(WENh$stage)] <- "yr"


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
  mutate(`sp` = ifelse(grepl("Chinook",`Release Site`),"Chinook", ifelse(grepl("Coho",`Release Site`),"Coho", "Steelhead"))) %>%
  mutate(sp = factor(sp, levels = c("Chinook", "Coho", "Steelhead"))) %>%
  mutate(`sp2` = ifelse(grepl("Chinook",`Release Site`), ifelse(grepl("HXH",`Release Site`),"Chinook HXH", "Chinook WXW"),
                                                               ifelse(grepl("Coho",`Release Site`),"Coho",
                                                                      ifelse(grepl("HXH",`Release Site`),"Steelhead HXH", "Steelhead WXW")))) %>%
  mutate(sp2 = factor(sp2, levels = c("Chinook HXH","Chinook WXW", "Coho", "Steelhead HXH", "Steelhead WXW"))) %>%
  mutate(`sp4` = ifelse(grepl("Steelhead",`Release Site`), "Steelhead",
                        ifelse(grepl("Chinook",`Release Site`),
                               ifelse(grepl("LFNH",`Release Site`),"Chinook LNFH", "Chinook Chiwawa and Nason"),
                               "Coho"))) %>%
  mutate(sp4 = factor(sp4, levels = c("Steelhead","Chinook LNFH", "Chinook Chiwawa and Nason", "Coho"))) %>%
  mutate(sp3 = `Release Site`) %>%
  mutate(sp3 = factor(sp3, levels = levels(`Release Site`))) %>%
  mutate(`Release Site` = factor(`Release Site`)) %>%
  droplevels()



names(WENh)[3] <- "ReleaseSite"


MR_settings <- list(state = 'LH',
                    frms = list(
                                phi = list(yr = ' ~ -1 + loc:sp4'), #survival
                                p = list(yr = ' ~ -1 + loc:sp'), #detection
                                lam = ' ~ 1', #nuisance parameter
                                eta = ' ~ -1'), #transition probability
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
# source("ggplot_pred_p.r")
# source("ggplot_pred_phi.r")
# g <- ggpubr::ggarrange(phi, p, nrow = 2, align = "v")
# print(g)
# source("ggplot_cum_pred_x.r")
