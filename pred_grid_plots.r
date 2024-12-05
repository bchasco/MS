est <- as.list(sd,"Estimate", report = TRUE)
est_sd <- as.list(sd,"Std. Error", report = TRUE)

df <- data.frame(fit@TMB$obj$tmb.data$pred.grid,
           reshape2::melt(fit@TMB$rep$gamma_pred),
           gam_est = reshape2::melt(est$gamma_pred)$value,
           gam_sd = reshape2::melt(est_sd$gamma_pred)$value,
           om_est = reshape2::melt(est$omega_pred)$value,
           om_sd = reshape2::melt(est_sd$omega_pred)$value) %>%
  filter(complete.cases(.)) %>%
  filter(loc == "as.Smolt")

  p <- df %>%
  # filter(fY %in% c("2001","2002")) %>%
  # mutate(fwk = factor((wk))) %>%
  ggplot(aes(x = loc, y = plogis(gam_est),
             color = "ReleaseSite")) +
  geom_point() +
  geom_errorbar(aes(ymin = plogis(gam_est - 1.96 * gam_sd),
                  ymax = plogis(gam_est + 1.96 * gam_sd),
                  # fill = ReleaseSite,
                  # group = fY,
                  # color = Length
                  ),
              alpha = 0.9,
              width = 0.2) +
  # facet_wrap(~ReleaseSite, ncol = 4) +
  ylim(0,1)+
  theme_bw()

print(p)

p <- df %>%
  filter(loc == "as.Smolt") %>%
  # mutate(fY = as.numeric(as.character(fY))) %>%
  ggplot(aes(x = loc,
             y = plogis(om_est))) +
  geom_point() +
  geom_errorbar(aes(ymin = plogis(om_est - 1.96 * om_sd),
                  ymax = plogis(om_est + 1.96 * om_sd),
                  # fill = ReleaseSite,
                  # group = fY,
                  # color = Length
  ),
  alpha = 0.9,
  width = 0.2) +
  # facet_wrap(~ReleaseSite, ncol = 5) +
  theme_bw()

print(p)
