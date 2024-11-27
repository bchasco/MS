library(dplyr)
library(tidyr)
library(ggplot2)
rm(list = ls())
load("data/LkWA_length.rda")


init_site = "tag"
last_site = "As.Adult.ballard"
locs <- c("tag", "as.Smolt","As.Adult.ballard")
LkWA$loc <- factor(LkWA$loc , levels = locs)
LkWA$Length <- as.numeric(LkWA$Length)
LkWA <- LkWA %>%
  mutate(state = ifelse(state == 0, "unk", "yr")) %>%
  mutate(state = factor(state, levels = c("yr", "unk"))) %>%
  filter(Year < 2022) %>%
  mutate(L2 = Length + rnorm(length(Length),0,10)) %>%
  mutate(fY = factor(Year, levels = sort(unique(Year)))) %>%
  mutate(fSp = ifelse(Species == "1", "Chinook", "Coho")) %>%
  mutate(fSp = factor(fSp, levels = c("Chinook", "Coho"))) %>%
  filter(fSp == "Chinook") %>%
  filter(ReleaseSite %in% c("LWCEDR","LWBEAR")) %>%
  mutate(Year = as.numeric(Year)) %>%
  droplevels() #%>%
  # filter(Length < 150 & Length > 70)

LkWA <- na.omit(LkWA)
LkWA$fL <- factor(round(LkWA$Length/5)*5, levels = sort(unique(round(LkWA$Length/5)*5)))

groups <- LkWA %>%
  group_by(ReleaseSite, fSp) %>%
  summarize(id  = min(as.integer(id)))

MR_settings <- list(state = 'state',
                    frms = list(
                                p = list(yr = ' ~ 1 '), #detection probability
                                phi = list(yr = ' ~ -1 + ReleaseSite + poly(Length,2) + fY'), #survival
                                lam = list(yr = ' ~ -1 + fY '), #nuisance parameter
                                eta = list(yr = ' ~ -1')), #transition probability
                    mod = "MS")


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
# print(plogis(fit@TMB$opt$par))
# print(exp(fit@TMB$opt$par)/sum(exp(1)+ exp(fit@TMB$opt$par)))


sd <- RTMB::sdreport(fit@TMB$obj)
# source("HMM_plot_phi.r")
# source("HMM_plot_p.r")

fr <- var_output("p", fit) #
fr %>%
  ggplot(aes(x = fY, y = est)) +
  geom_point() +
  # facet_wrap(~ReleaseSite, ncol = 1) +
  theme_bw()


