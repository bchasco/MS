state <- "yr"
# loc <- "as.Smolt"
loc <- "As.Adult.ballard"

state_id <- which(levels(tmb.data$data$state)%in%state, arr.ind = TRUE)


summary_df <- data.frame(fit@TMB$obj$tmb.data$data[fit@TMB$obj$tmb.data$data$loc %in%
                                                     c("as.Smolt","As.Adult.ballard"),],
                         est = fit@TMB$rep$tau_pred[state_id,]) %>%
  filter(loc == !!loc & state == !!state) %>%
  group_by(loc, state, wk) %>%
  summarize(mean_time = sum((time*n)/sum(n)), .groups = "drop") %>%
  ungroup()

xx <- data.frame(fit@TMB$obj$tmb.data$data[fit@TMB$obj$tmb.data$data$loc %in%
                                             c("as.Smolt","As.Adult.ballard"),],
                 est = fit@TMB$rep$tau_pred[state_id,]) %>%
  filter(loc == !!loc & state == !!state) %>%
  left_join(summary_df, by = c("wk", "loc", "state")) %>%
  ggplot(aes(x = as.factor(wk), y = time)) +
  geom_boxplot() +
  stat_summary(aes(group = 1), colour = "black",
               geom = "line", fun.y = "median") +   # geom_point(aes(x = wk, y = mean_time, color = "red", size = 5)) +
  # geom_point(aes(x = as.factor(wk), y = mean_time, color = "blue", size = 5)) +
  geom_point(aes(x = as.factor(wk), y = est, color = "green", size = 5)) +
  ylab("Travel time in weeks") +
  xlab("Tagging date") +
  theme_bw()
  # facet_grid(~wk)

  print(xx)

