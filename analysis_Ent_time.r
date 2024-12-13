library(dplyr)
library(tidyr)
library(ggplot2)
# rm(list = ls())
load("data/Ent.rda")


init_site = "Lwr_Main_TRAP"
last_site = "BON"
locs <- c("First_Trap","Lwr_Main_TRAP","RRJ","BON")
Ent$loc <- factor(Ent$loc , levels = locs)
# LkWA$Length <- as.numeric(LkWA$Length)
Ent <- Ent %>%
  mutate(time = as.numeric(as.character(cum_time))) %>%
  mutate(state = ifelse(is.na(stage), "unk", stage)) %>%
  mutate(state = factor(state, levels = c("sub", "yr", "unk"))) %>%
  # mutate(wk = as.numeric(as.factor(init_week)) ) %>%
  group_by(id) %>%
  mutate(wk = min(na.omit(as.numeric(init_week)))) %>%
  ungroup() %>%
  ungroup() %>%
  mutate(logwk = log(wk)) %>%
  mutate(fwk = factor(wk, levels = sort(unique(wk)))) %>%
  group_by(id) %>%
  mutate(ReleaseSite = na.omit(unique(`Event Species Name`))) %>%
  ungroup() %>%
  droplevels()

Ent_mat <- Ent %>%
  group_by(state=="RRJ")

MR_settings <- list(state = 'state',
                    frms = list(
                                phi = list(sub = ' ~ 1',
                                           yr = '~ 1'), #survival
                                p = list(yr = ' ~ 1 '), #detection probability
                                lam = list(yr = ' ~ 1  '), #nuisance parameter
                                eta = list(yr = ' ~ 1'), #transition probability
                                time = list(yr = ' ~ -1 + ReleaseSite + fwk') #travel time
                                ), #transition probability

                    mod = "MSt")


# Create an instance of the tmb_list class
input <- new("tmb_list",
                   MR_settings = MR_settings,
                   data = Ent)

#take the input and put it into a tmb model
fit <- MStmb(input)
print(fit@TMB$AICc)
df <- plot_data(fit)

# Plot
ggplot(df, aes(x = arrival_day, y = value, color = est, group = est)) +
  geom_point(data = df %>% filter(est == "obs_n"), size = 3) +  # Points for obs_n
  geom_line(data = df %>% filter(est == "pred_n"), size = 1) +  # Lines for pred_n
  facet_wrap(~ ReleaseSite, ncol = 1, scales = "free_y") +
  labs(
    title = "Observation and Prediction by Release Site",
    x = "Arrival Day",
    y = "Value",
    color = "Estimate Type"
  ) +
  theme_minimal()



