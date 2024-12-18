
# err <- as.list(sd,"Std. Error", report=TRUE)

fit@TMB$obj$tmb.data$data %>%
  filter(!is.na(time) & time!=0) %>%
  filter(state == "yr") %>%
  filter(loc == "as.Smolt") %>%
  filter(time>40) %>%
  group_by(fwk) %>%
  summarise(mean = mean(time)) %>%
  ggplot(aes(x = fwk, y = mean)) +
  geom_point()

xx <- cbind(fit@TMB$obj$tmb.data$pred.grid[fit@TMB$obj$tmb.data$pred.grid$loc==fit@TMB$obj$tmb.data$MR_settings$dv$loc,],
      est = fit@TMB$rep$tau__pred
      # ,err = err$gamma_pred
      ) %>%
  # filter(wk == 10, Length == 90) %>%
  filter(state == "yr")

print(plogis(xx$est))

xx %>%
  ggplot(aes(x = as.numeric(as.character(wk))
             ,y = plogis(est)
             ,color = as.factor(fY)
             )) +
  geom_line() +
  geom_ribbon(aes(ymax = plogis(est + 1.96 * err),
                  ymin = plogis(est - 1.96 * err)
                  ,fill = as.factor(fY)
  ),
  color = NA,
  alpha = 0.2) +
  # facet_grid(as.factor(wk)~as.factor(Length)) +
  theme_classic()

print(xx)


