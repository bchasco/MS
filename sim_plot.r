draw <- data.frame(sim = rmultinom(1,sum(fit@TMB$rep$nn),fit@TMB$rep$pp/sum(fit@TMB$rep$pp)),
      obs = fit@TMB$rep$nn) %>%
  group_by(obs) %>%
  summarize(mean_sim = mean(sim))
p <- draw %>%
  ggplot(aes(x = obs, y = obs)) +
  geom_point() +
  geom_point(aes(x = obs, y = mean_sim, color = "red")) +
  xlim(0, 50) +
  ylim(0, 50)

print(p)

# plot(draw$obs, draw$sim)

