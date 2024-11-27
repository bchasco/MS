library(dplyr)
library(tidyr)
library(ggplot2)
# rm(list = ls())
load("data/wildFish_sub.rda")
wild_sub <- WENh

locs <- c("Trib", "WEN",  "Downstream")
eta_locs <- c("Trib",  "WEN", "Downstream")
last_site <- "Downstream"
states <- c("sub","yr", "unk")
# init_site <- 'Trib'

wild_sub <- wild_sub %>%
  mutate(state = ifelse(is.na(stage),"unk",stage)) %>%
  mutate(state = factor(state, levels = states)) %>%
  mutate(loc2 = ifelse(loc %in% locs, loc, "Downstream")) %>%
  mutate(eta_locs = ifelse(loc %in% eta_locs, loc, "Downstream")) %>%
  group_by(id) %>%
  mutate(stage2 = ifelse(loc %in% locs, stage,
                         ifelse(stage%in%c("sub","yr"), "yr", stage))) %>%
  group_by(id,loc2) %>%
  select(-loc) %>%
  distinct() %>%
  group_by(id, loc2) %>%
  filter(!is.na(stage2) | n() == 1) %>% # Keep non-NA stage2 or single-row groups
  arrange(desc(stage2 == "yr")) %>% # Prioritize "yr" in stage2
  slice(1) %>% # Keep the first row per group
  ungroup() %>%
  mutate(loc2 = factor(loc2, levels = locs)) %>%
  group_by(id, loc2) %>%
  arrange(id,loc2) %>%
  filter(loc2 %in% locs) %>%
  droplevels() %>%
  mutate(state = ifelse(is.na(stage2),"unk",stage2)) %>%
  mutate(state = factor(state, levels = c("sub","yr","unk"))) %>%
  # mutate(state = ifelse(is.na(stage2),"unk","sub")) %>%
  # mutate(state = factor(state, levels = c("sub","unk"))) %>%
  mutate(id = factor(id)) %>%
  mutate(loc = loc2) %>%
  ungroup()

load("data/wildFish_yr.rda")
wild_yr <- WENh
wild_yr$stage[wild_yr$stage == "sub"] <- "yr"


wild_yr <- wild_yr %>%
  mutate(state = ifelse(is.na(stage),"unk",stage)) %>%
  mutate(state = factor(state, levels = states)) %>%
  mutate(loc2 = ifelse(loc %in% locs, loc, "Downstream")) %>%
  mutate(eta_locs = ifelse(loc %in% eta_locs, loc, "Downstream")) %>%
  group_by(id) %>%
  mutate(stage2 = ifelse(loc %in% locs, stage,
                         ifelse(stage%in%c("sub","yr"), "yr", stage))) %>%
  group_by(id,loc2) %>%
  select(-loc) %>%
  distinct() %>%
  group_by(id, loc2) %>%
  filter(!is.na(stage2) | n() == 1) %>% # Keep non-NA stage2 or single-row groups
  arrange(desc(stage2 == "yr")) %>% # Prioritize "yr" in stage2
  slice(1) %>% # Keep the first row per group
  ungroup() %>%
  mutate(loc2 = factor(loc2, levels = locs)) %>%
  group_by(id, loc2) %>%
  arrange(id,loc2) %>%
  filter(loc2 %in% locs) %>%
  droplevels() %>%
  mutate(state = ifelse(is.na(stage2),"unk",stage2)) %>%
  mutate(state = factor(state, levels = c("yr","unk"))) %>%
  # mutate(state = ifelse(is.na(stage2),"unk","sub")) %>%
  # mutate(state = factor(state, levels = c("sub","unk"))) %>%
  mutate(id = factor(id)) %>%
  mutate(loc = loc2) %>%
  ungroup()

names(wild_yr)[names(wild_yr)=='Release Site'] <- "ReleaseSite"
wild_yr$id <- factor(as.numeric(wild_yr$id) + max(as.numeric(wild_sub$id)))

names(wild_sub)[names(wild_sub)=='Release Site'] <- "ReleaseSite"

wild <- rbind(wild_sub) %>%
  droplevels()

wild$state[wild$state!="unk"] <- "sub"
wild <- wild %>%
  droplevels()
MR_settings <- list(state = 'state',
                    frms = list(
                                p = list(sub_yr = ' ~ 1'), #survival
                                phi = list(sub_yr = ' ~ 1'), #detection
                                lam = list(sub_yr = ' ~ 1 '), #nuisance parameter
                                eta = list(sub_yr = ' ~ -1')), #transition probability
                    mod = "CJS")


# Create an instance of the tmb_list class
input <- new("tmb_list",
                   MR_settings = MR_settings,
                   data = wild)

#take the input and put it into a tmb model
fit <- MStmb(input)

# Assuming `fit@TMB$opt$par` is your named vector
params <- fit@TMB$opt$par

# Split the vector into a list by names
grouped_params <- split(params, names(params))
print(fit@TMB$opt$par)

# sd <- RTMB::sdreport(fit@TMB$obj)
# source("HMM_plot_phi.r")
# source("HMM_plot_p.r")

