library(dplyr)
library(tidyr)
library(ggplot2)
rm(list = ls())
load("data/LkWA_time.rda")
#
# load("data/LkWA_length.rda")


init_site = "tag"
last_site = "As.Adult.ballard"
locs <- c("tag", "as.Smolt","As.Adult.ballard")
LkWA$loc <- factor(LkWA$loc , levels = locs)
# LkWA$Length <- as.numeric(LkWA$Length)
LkWA <- LkWA %>%
  mutate(state = ifelse(state == 0, "unk", "yr")) %>%
  mutate(state = factor(state, levels = c("yr", "unk"))) %>%
  filter(Year < 2022) %>%
  mutate(fY = factor(Year, levels = sort(unique(Year)))) %>%
  mutate(fSp = ifelse(Species == "1", "Chinook", "Coho")) %>%
  mutate(fSp = factor(fSp, levels = c("Chinook", "Coho"))) %>%
  # mutate(wk = as.numeric(as.factor(ReleaseWeek)) ) %>%
  # mutate(fwk = factor(wk, levels = sort(unique(wk)))) %>%
  # mutate(L2 = round(as.numeric(as.factor(Length))/5)*5) %>%
  filter(fSp == "Chinook") %>%
  # filter(ReleaseSite %in% c("LWCEDR","LWBEAR")) %>%
  mutate(Year = as.numeric(as.character(Year))) %>%
  droplevels() #%>%
  # filter(Length < 150 & Length > 70)

LkWA <- na.omit(LkWA)
# LkWA$fL <- factor(round(LkWA$Length/5)*5, levels = sort(unique(round(LkWA$Length/5)*5)))

MR_settings <- list(state = 'state',
                    frms = list(
                                phi = list(yr = ' ~ 1'), #survival
                                p = list(yr = ' ~ 1'), #detection probability
                                lam = list(yr = ' ~ 1'), #nuisance parameter
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


sd <- RTMB::sdreport(fit@TMB$obj)
# source("HMM_plot_phi.r")
# source("HMM_plot_p.r")

# fr <- var_output("p", fit) #
# fr %>%
#   ggplot(aes(x = fY, y = est)) +
#   geom_point() +
#   # facet_wrap(~ReleaseSite, ncol = 1) +
#   theme_bw()


