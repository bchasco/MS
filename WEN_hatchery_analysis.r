library(dplyr)
library(tidyr)
rm(list = ls())
data("WEN_hatchery")
WENh <- WENh %>%
  filter(loc %in% c("Trib","First_Trap","WEN","MCJ","BON","TWX_EST"))

# WENh <- WENh %>%
#   mutate(state = ifelse(is.na(stage),"unk",stage))
# Define MR_settings
MR_settings <- list(state = 'LH',
                            frms = list(phi = ' ~ loc', #survival
                                        p = ' ~ loc', #detection
                                        lam = ' ~ 1', #nuisance parameter
                                        eta = ' ~ loc_eta'), #transition probability
                            dm = list(),
                            mod = "MS")


phi_m <- list()
for(i in c("yr")){
  tmp_WEN <- WENh %>%
    mutate(state = ifelse(is.na(stage),"unk","yr")) %>%
    mutate(state = factor(state, levels = c("yr","unk"))) %>%
    mutate(loc = factor(loc, levels = c("Trib","First_Trap","WEN","MCJ","JDJ","BON","TWX_EST"))) %>%
    mutate(loc_phi = factor(loc, levels = c("Trib","First_Trap","WEN","MCJ","JDJ","BON","TWX_EST"))) %>%
    droplevels() %>%
    # mutate(state = !!i) %>%
    mutate(tag_site = "Trib",
           last_site = "TWX_EST") %>%
    filter(loc != tag_site,
           loc != last_site) %>%
    droplevels()


  #For each state
  phi_m[[i]] <-  tryCatch(Matrix::sparse.model.matrix(formula(MR_settings$frm$phi), tmp_WEN),
                     error = function(e) e,
                     warning = function(w) w)
  if(class(phi_m[[i]])!="dgCMatrix"){
    phi_m[[i]] <- matrix(rep(1, nrow(tmp_WEN)))
  }
}
MR_settings$dm$phi <- phi_m

eta_m <- list()
for(i in c("sub")){
  tmp_WEN <- WENh %>%
    mutate(state = ifelse(is.na(stage), "unk", stage)) %>%
    mutate(state = factor(state, levels = c("yr","unk"))) %>%
    mutate(loc = factor(loc, levels = c("Trib","First_Trap","WEN","MCJ","JDJ","BON","TWX_EST"))) %>%
    mutate(loc_eta = case_when(
      i == "sub" & loc %in% c("Wen","MCJ","JDJ", "BON") ~ "Wen",
      i == "sub" ~ as.character(loc)
    )) %>%
    mutate(loc_eta = factor(loc_eta, levels = c("Trib","First_Trap","WEN","MCJ","JDJ","BON","TWX_EST"))) %>%
    droplevels() %>%
    mutate(state = !!i) %>%
    mutate(tag_site = "Trib",
           last_site = "TWX_EST")

  tmp_WEN <- tmp_WEN[tmp_WEN$loc!=tmp_WEN$tag_site &
                       tmp_WEN$loc!=tmp_WEN$last_site,] %>%
    droplevels()

  eta_m[[i]] <-  tryCatch(Matrix::sparse.model.matrix(formula('state ~ loc_eta'), tmp_WEN),
                          error = function(e) e,
                          warning = function(w) w)
}
MR_settings$dm$eta <- eta_m

WEN_data <- WENh %>%
  mutate(state = ifelse(is.na(stage), "unk", "yr")) %>%
  mutate(state = factor(state, levels = c("yr","unk"))) %>%
  mutate(loc = factor(loc, levels = c("Trib","First_Trap","WEN","MCJ","JDJ","BON","TWX_EST"))) %>%
  mutate(tag_site = "Trib",
         last_site = "MCJ") %>%
  group_by(state) %>%
  droplevels()

WEN_data



# Create an instance of the tmb_list class
input <- new("tmb_list",
                   MR_settings = MR_settings,
                   data = WEN_data)


#take the input and put it into a tmb model
fit <- MStmb(input)
# plot(fit)
sd <- RTMB::sdreport(fit@TMB$obj)
round((fit@TMB$rep$gamma),4)
round((fit@TMB$rep$omega),2)
