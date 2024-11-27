det <- fit@TMB$obj$tmb.data$dm$grid_p
group <- names(det)
group <- group[group!="loc"]

det$qest <- round(as.array(as.list(sd,"Est",report = TRUE), drop = FALSE)$pred_p,3)
det$sd_for_logit <- round(as.array(as.list(sd,"Std. Error",report = TRUE), drop = FALSE)$pred_p,3)
det$est <- round(plogis(det$qest),3)
det$lwr <- round(plogis(det$qest - 1.96 * det$sd_for_logit),3)
det$upr <- round(plogis(det$qest + 1.96 * det$sd_for_logit),3)

ymax <- max(det$ymax[det$loc =="WEN"])
p <- as.data.frame(det) %>%
  # filter(loc == "WEN") %>%
  # group_by(loc) %>%
  mutate(ymax_lim = max(ymax)) #
if(length(group)>0){
  p <- p  %>%
    filter(loc == "WEN") %>%
    ggplot2::ggplot(aes(y = est, x = loc , group = !!sym(group), color = !!sym(group))) +
    theme_bw() +
    # facet_wrap(~sp) +
    ylab('Detection probability') +
    geom_point(size = 5, position = position_dodge(width = 1)) +
    geom_errorbar(aes(ymin = lwr, ymax = upr),
                  position = position_dodge(width = 1), size = 1.5, width = 0.15) +
    theme_classic() +
    xlab("Detection location") +
    theme(text = element_text(size = 16)) +
    ylim(0,0.65) +
    labs(color='Fish grouping')
}else{
  p <- p  %>%
    filter(loc == "WEN") %>%
    ggplot2::ggplot(aes(y = est, x = loc )) +
    theme_bw() +
    # facet_wrap(~sp) +
    ylab('Detection probability') +
    geom_point(size = 5, position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = lwr, ymax = upr),
                  position = position_dodge(width = 0.5), size = 1.5, width = 0.15) +
    theme_classic() +
    xlab("Detection location") +
    theme(text = element_text(size = 16)) +
    ylim(0,1) +
    labs(color='Fish grouping')
}
print(p)
print(det)
