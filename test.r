library(dplyr)
library(tidyr)
rm(list = ls())
data("WEN")

# Define MR_settings
MR_settings <- list(state = 'state',
                            frms = list(phi = ' ~ loc',
                                        p = ' ~ loc',
                                        lam = ' ~ 1',
                                        eta = ' ~ loc_eta'),
                            dm = list(phi = NA,
                                                 eta = NA,
                                                 p = NA),
                            mod = "MS")

phi_m <- list()
for(i in c("sub","yr")){
  tmp_WEN <- WEN %>%
    mutate(state = ifelse(is.na(stage), "unk", stage)) %>%
    mutate(state = factor(state, levels = c("sub","yr","unk"))) %>%
    mutate(loc = factor(loc, levels = c("Trib","First_Trap","Wen","MCJ","JDJ","BON","TWX_EST"))) %>%
    mutate(loc_phi = case_when(
      i == "sub" & loc %in% c("Wen","MCJ", "JDJ", "BON") ~ "Wen",
      i == "sub" ~ as.character(loc),
      i == "yr" & loc %in% c("Trib", "First_Trap", 'Wen') ~ "Wen",
      i == "yr" ~ as.character(loc)
    )) %>%
    mutate(loc_phi = factor(loc_phi, levels = c("Trib","First_Trap","Wen","MCJ","JDJ","BON","TWX_EST"))) %>%
    droplevels() %>%
    mutate(state = !!i) %>%
    mutate(tag_site = "Trib",
           last_site = "TWX_EST") %>%
    filter(loc != tag_site,
           loc != last_site) %>%
    droplevels()

  MR_settings$dm$phi[[i]] <-  tryCatch(Matrix::sparse.model.matrix(formula(MR_settings_example$frm$phi), tmp_WEN),
                     error = function(e) e,
                     warning = function(w) w)
}

eta_m <- list()
for(i in c("sub")){
  tmp_WEN <- WEN %>%
    mutate(state = ifelse(is.na(stageID), "unk", stageID)) %>%
    mutate(state = factor(state, levels = c("sub","yr","unk"))) %>%
    mutate(loc = factor(loc, levels = c("Trib","First_Trap","Wen","MCJ","JDJ","BON","TWX_EST"))) %>%
    mutate(loc_eta = case_when(
      i == "sub" & loc %in% c("Wen","MCJ","JDJ", "BON") ~ "Wen",
      i == "sub" ~ as.character(loc)
    )) %>%
    mutate(loc_eta = factor(loc_eta, levels = c("Trib","First_Trap","Wen","MCJ","JDJ","BON","TWX_EST"))) %>%
    droplevels() %>%
    mutate(state = !!i) %>%
    mutate(tag_site = "Trib",
           last_site = "TWX_EST")

  tmp_WEN <- tmp_WEN[tmp_WEN$loc!=tmp_WEN$tag_site &
                       tmp_WEN$loc!=tmp_WEN$last_site,] %>%
    droplevels()

  MR_setting$dm$eta[[i]] <-  tryCatch(Matrix::sparse.model.matrix(formula('state ~ loc_eta'), tmp_WEN),
                          error = function(e) e,
                          warning = function(w) w)
}



# Create an instance of the tmb_list class
input <- new("tmb_list",
                   MR_settings = MR_settings_example,
                   data = WEN)


#take the input and put it into a tmb model
fit <- MStmb(input)
sd <- RTMB::sdreport(fit@TMB$obj)

print(sapply(fit@TMB$rep$phi_s, function(x){round(plogis(x),2)}))
print(sapply(fit@TMB$rep$eta, function(x){round(plogis(x),2)}))

print(fit@TMB$rep$eta_s)
print(fit@TMB$rep$phi_s)

