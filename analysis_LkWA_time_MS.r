library(dplyr)
library(tidyr)
library(ggplot2)
rm(list = ls())
# load("data/LkWA.rda")
load("data/LkWA_time.rda")
# load("data/LkWA_time_releasewk_length.rda")




init_site = "tag"
last_site = "As.Adult.ballard"
locs <- c("tag", "as.Smolt","As.Adult.ballard")
LkWA$loc <- factor(LkWA$loc , levels = locs)

LkWA <- LkWA %>%
  # mutate(time := ifelse(!("time" %in% names(.)), 1, time)) %>%
  mutate(time = as.numeric(as.character(time))) %>%
  # mutate(state = ) %>%
  mutate(state = ifelse(state == 0, "unk",
                        ifelse(loc == "as.Smolt" & state != 0 & time < 40, "yr",
                               ifelse(loc == "as.Smolt" & state != 0 & time >= 40,"yr_plus","yr")))) %>%
  mutate(state = factor(state, levels = c("yr", "yr_plus","unk"))) %>%
  group_by(id) %>%
  mutate(init_state = first(state)) %>%
  mutate(fY = factor(Year, levels = sort(unique(Year)))) %>%
  mutate(fSp = ifelse(Species == "1", "Chinook", "Coho")) %>%
  mutate(fSp = factor(fSp, levels = c("Chinook", "Coho"))) %>%
  filter(Year < 2022 & Year>2001 ) %>%
  mutate(Year  = as.numeric(as.character(Year))) %>%
  filter(fSp %in% c("Chinook")) %>%
  filter(ReleaseSite %in% c("LWCEDR")) %>%
  droplevels()


MR_settings <- list(state = 'state',
                    frms = list(
                                phi = list(yr = ' ~ 1 + (1|fY)',
                                           yr_plus = ' ~ 1'), #survival
                                p = list(yr = ' ~ 1 + (1|fY)'), #detection probability
                                lam = list(yr = ' ~ 1'), #nuisance parameter
                                eta = list(yr = ' ~ 1'), #transition probability
                                tau = list(yr = ' ~ -1 + loc + (1|fY)',
                                           yr_plus = '~ 1') #travel time
                                ), #transition probability
                    dv = list(state = "yr",
                              loc = "as.Smolt"
                    ),
                    mod = "MSt")


# Create an instance of the tmb_list class
input <- new("tmb_list",
                   MR_settings = MR_settings,
                   data = LkWA)


#take the input and put it into a tmb model
fit <- MStmb(input)
print(fit@TMB$AICc)

sd <- RTMB::sdreport(fit@TMB$obj)
# save(fit, file = "fit_Length_wk_year_no_re_time_year_wk_no_re.rda")
