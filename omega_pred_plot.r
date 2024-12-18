
err <- as.list(sd,"Std. Error", report=TRUE)

xx <- cbind(fit@TMB$obj$tmb.data$pred.grid[fit@TMB$obj$tmb.data$pred.grid$loc=="as.Smolt",],
      est = fit@TMB$rep$omega_pred,
      err = err$omega_pred) %>%
  filter(wk == 10, Length == 90) %>%
  filter(state == "yr")

print(round(plogis(xx$est),2))

xx %>%
  ggplot(aes(x = as.numeric(as.character(fY)), y = plogis(est))
         # , color = ReleaseSite
         ) +
  geom_line() +
  geom_ribbon(aes(ymax = plogis(est + 1.96 * err),
                  ymin = plogis(est - 1.96 * err),
                  fill = as.factor(ReleaseSite)
                  ),
              color = NA,
              alpha = 0.2) +
  facet_grid(as.factor(wk)~as.factor(Length)) +
  theme_classic()

print(xx)


