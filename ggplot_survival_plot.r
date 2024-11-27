nm <- dimnames(fit@TMB$rep$cum_gamma_pred)
gamma_pred_est <- as.array(as.list(sd,"Est",report = TRUE), drop = FALSE)$gamma_pred
gamma_pred_sd <- as.array(as.list(sd,"Std. Error",report = TRUE), drop = FALSE)$gamma_pred
dimnames(gamma_pred_est) <- nm
surv <- reshape2::melt(gamma_pred_est[,,,], drop = FALSE)
surv$sd <- reshape2::melt(gamma_pred_sd[,,,], drop = FALSE)$value
if(!("WEN" %in% levels(surv$Locations))){
  surv$Locations <- factor(surv$Locations, levels = c(locs[1],"WEN",locs[2:length(locs)]))
  surv$Locations[nrow(surv)] <- "WEN"
  surv$current_state[nrow(surv)] <- "yr"
  surv$next_state[nrow(surv)] <- "yr"
  surv$value[nrow(surv)] <- NA
  surv$sd[nrow(surv)] <- NA

}

p <- as.data.frame(surv) %>%
  filter(next_state=="yr",
         current_state=="yr") %>%
  filter(Locations!="Trib" & Locations!=last_site) %>%
  droplevels() %>%
  # tibble::rownames_to_column(var = "Location") %>%
  # mutate(Locations = factor(Locations, levels = locs)) %>%
  ggplot2::ggplot(aes(y = plogis(value), x = Locations, group = rel_site, color = rel_site)) +
  theme_bw() +
  # facet_wrap(~sp) +
  ylab('Cumulative survival') +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = plogis(value - 1.96 * sd), ymax = plogis(value + 1.96 * sd)),
                position = position_dodge(width = 0.5), width = 0.25) +
    theme_classic() +
  xlab("Survival") +
  theme(text = element_text(size = 16)) +
  ylim(0,1)

print(p)
