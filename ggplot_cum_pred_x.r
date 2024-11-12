surv <- fit@MR_settings$pred_data
surv$cum_est <- as.array(as.list(sd,"Est",report = TRUE), drop = FALSE)$cum_pred_x
surv$cum_sd <- as.array(as.list(sd,"Std. Error",report = TRUE), drop = FALSE)$cum_pred_x
if(!("WEN" %in% levels(surv$loc))){
  surv$loc <- factor(surv$loc, levels = c(locs[1],"WEN",locs[2:length(locs)]))
  surv$loc[nrow(surv)] <- "WEN"
  surv$current_state[nrow(surv)] <- "yr"
  surv$next_state[nrow(surv)] <- "yr"
  surv$cum_est[nrow(surv)] <- NA
  surv$cum_sd[nrow(surv)] <- NA

}

p <- as.data.frame(surv) %>%
  filter(next_state=="yr",
         current_state=="yr") %>%
  # filter(loc!="Trib" & loc!=last_site) %>%
  droplevels() %>%
  # tibble::rownames_to_column(var = "Location") %>%
  # mutate(Locations = factor(Locations, levels = locs)) %>%
  ggplot2::ggplot(aes(y = plogis(cum_est), x = loc, group = ReleaseSite, color = ReleaseSite)) +
  theme_bw() +
  # facet_wrap(~sp) +
  ylab('Cumulative survival') +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = plogis(cum_est - 1.96 * cum_sd), ymax = plogis(cum_est + 1.96 * cum_sd)),
                position = position_dodge(width = 0.5), width = 0.25) +
  theme_classic() +
  xlab("Detection location") +
  theme(text = element_text(size = 16)) +
  ylim(0,1)

print(p)
