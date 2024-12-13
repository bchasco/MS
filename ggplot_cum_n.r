plot_data <- function(my_fit){
  df <- cbind(tmb.data$data,
              cum_phi = my_fit@TMB$rep$data_est,
              omega = my_fit@TMB$rep$data_omega,
              cum_n = my_fit@TMB$rep$data_est_n,
              pred_time = my_fit@TMB$rep$data_tau) %>%
    mutate(obs_n = n*!is.na(time)) %>%
    mutate(pred_n = n*omega*cum_phi) %>%
    select(loc,time,ReleaseSite,n,id,wk,pred_time,pred_n,obs_n) %>%
    filter(loc == "RRJ") %>%
    mutate(arrival_day = wk + pred_time/7) %>%
    mutate(obs_day = wk + time/7) %>%
    group_by(ReleaseSite,arrival_day) %>%
    summarise(obs_n = sum(obs_n),
              pred_n = sum(pred_n)) %>%
    pivot_longer(cols = c(obs_n,pred_n), names_to = "est") #%>%

  df <- df %>%
    mutate(
      est = factor(est, levels = c("obs_n", "pred_n"))
    )
  return(df)
}

