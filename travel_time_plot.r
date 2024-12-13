
err <- as.list(sd,"Std. Error", report=TRUE)

xx <- cbind(fit@TMB$obj$tmb.data$pred.grid[fit@TMB$obj$tmb.data$pred.grid$loc=="as.Smolt",],
      est = fit@TMB$rep$omega_pred,
      err = err$omega_pred) %>%
  filter(state == "yr") %>%
  ggplot(aes(x = as.numeric(wk), y = plogis(est))) +
  geom_line() +
  geom_ribbon(aes(ymax = plogis(est + 1.96 * err),
                  ymin = plogis(est - 1.96 * err),
                  # fill = as.factor(fY)
                  ),
              color = NA,
              alpha = 0.2) +
  facet_wrap(~fY, ncol = 4) +
  theme_classic()

print(xx)


