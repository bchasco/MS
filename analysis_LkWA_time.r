library(dplyr)
library(tidyr)
library(ggplot2)
# rm(list = ls())
load("data/LkWA_time.rda")
#
# load("data/LkWA_length.rda")


init_site = "tag"
last_site = "As.Adult.ballard"
locs <- c("tag", "as.Smolt","As.Adult.ballard")
LkWA$loc <- factor(LkWA$loc , levels = locs)
# LkWA$Length <- as.numeric(LkWA$Length)
LkWA <- LkWA %>%
  mutate(time = as.numeric(as.character(time))) %>%
  mutate(state = ifelse(state == 0, "unk", "yr")) %>%
  mutate(state = factor(state, levels = c("yr", "unk"))) %>%
  mutate(fY = factor(Year, levels = sort(unique(Year)))) %>%
  mutate(fSp = ifelse(Species == "1", "Chinook", "Coho")) %>%
  mutate(fSp = factor(fSp, levels = c("Chinook", "Coho"))) %>%
  # mutate(wk = as.numeric(as.factor(ReleaseWeek)) ) %>%
  # mutate(fwk = factor(wk, levels = sort(unique(wk)))) %>%
  # mutate(Length = as.numeric(as.character(Length))) %>%
  mutate(Year = as.numeric(as.character(Year))) %>%
  # filter(Length < 150 & Length > 70) %>%
  filter(Year < 2022, Year >= 2015) %>%
  filter(fSp == "Chinook") %>%
  filter(ReleaseSite %in% c("LWCEDR","LWBEAR")) %>%
  # filter(!is.na(Length)) %>%
  droplevels()


MR_settings <- list(state = 'state',
                    frms = list(
                                phi = list(yr = ' ~ -1 + ReleaseSite + (1|fY)'), #survival
                                p = list(yr = ' ~ 1 + (1|fY)'), #detection probability
                                lam = list(yr = ' ~ -1 + fY'), #nuisance parameter
                                eta = list(yr = ' ~ -1')), #transition probability
                    mod = "MSt")


# Create an instance of the tmb_list class
input <- new("tmb_list",
                   MR_settings = MR_settings,
                   data = LkWA)

#take the input and put it into a tmb model
fit <- MStmb(input)

# Assuming `fit@TMB$opt$par` is your named vector
params <- fit@TMB$opt$par

# Split the vector into a list by names
grouped_params <- split(params, names(params))


