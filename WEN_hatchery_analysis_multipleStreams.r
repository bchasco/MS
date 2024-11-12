library(dplyr)
library(tidyr)
rm(list = ls())
data("WEN_hatchery_multipleStreams")
WENh <- WENh %>%
  filter(loc %in% c("Trib","First_Trap","WEN","MCJ","BON","TWX_EST"))
names(WENh)[3] <- "ReleaseSite"
# WENh <- WENh %>%
#   filter(ReleaseSite == "NASONR")

# Define MR_settings
# Get unique values for each of the first three columns
unique_locs <- unique(WENh$loc)
unique_sites <- unique(WENh$ReleaseSite)

# Use expand.grid to create all combinations
expanded_grid <- data.frame(expand.grid(
  loc = unique_locs,
  ReleaseSite = unique_sites
))


MR_settings <- list(state = 'LH',
                    frms = list(phi = ' ~ loc : ReleaseSite', #survival
                                p = ' ~ loc : ReleaseSite', #detection
                                lam = ' ~ 1', #nuisance parameter
                                eta = ' ~ loc_eta'), #transition probability
                    dm = list(),
                    pdm = list(),
                    pred_data = expanded_grid,
                    mod = "MS")


phi_m <- list()
pred_m <- list()
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
  phi_m[[i]] <-  unlist(tryCatch(Matrix::sparse.model.matrix(formula(MR_settings$frm$phi), tmp_WEN),
                     error = function(e) e,
                     warning = function(w) w))

  names_phi <- colnames(phi_m[[i]])

  pred_data <- MR_settings$pred_data %>%
    mutate(tag_site = "Trib", last_site = "TWX_EST") %>%
    filter(loc != tag_site,loc != last_site)

  pred_m[[i]] <-  unlist(tryCatch(Matrix::sparse.model.matrix(formula(MR_settings$frm$phi), pred_data),
                                 error = function(e) e,
                                 warning = function(w) w))
  names_pred <- colnames(pred_m[[i]])

  pred_m <- pred_m[[i]][,names_pred %in% names_phi]
  # phi_m[[i]] <- as.data.frame(phi_m[[i]])
  # if(class(phi_m[[i]])!="dgCMatrix"){
  #   phi_m[[i]] <- matrix(rep(1, nrow(tmp_WEN)))
  # }
}
MR_settings$pred_data <- pred_data
MR_settings$dm$phi <- phi_m
MR_settings$dm$pred_phi <- pred_m

p_m <- list()
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
  p_m[[i]] <-  unlist(tryCatch(Matrix::sparse.model.matrix(formula(MR_settings$frm$p), tmp_WEN),
                                 error = function(e) e,
                                 warning = function(w) w))

}
MR_settings$dm$p <- p_m

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
         last_site = "TWX_EST") %>%
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
round((fit@TMB$rep$gamma_pred),3)
