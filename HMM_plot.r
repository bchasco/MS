library(ggplot2)
library(ggforce)
library(dplyr)
library(reshape2)

# Simulate some example data (use your actual data)
# x <- tmb_model@TMB$rep$gamma
x <- reshape2::melt(fit@TMB$rep$gamma)
x <- x[x$i%in%c(1),]
x$ReleaseSite <- fit@TMB$obj$tmb.data$data$ReleaseSite[match(x$i,fit@TMB$obj$tmb.data$data$id)]
# x <- x[x$ReleaseSite=="Chinook_Chiwawa_W",]
x$next2 <- x$current_state
x$init2 <- x$next_state

# Process the data for plotting
g <- x %>%
  filter(value >=0, value <1) %>%
  mutate(loc_id = as.numeric(factor(j, labels = locs))) %>%
  # filter(paste(next2,init2)!="unk unk",
  #        init2!="unk") %>%
  mutate(x1 = loc_id -1,
         x2 = ifelse(next2 == "unk",loc_id-1,loc_id  - 0.1),
         y1 = ifelse(init2 == "sub",1.0,2),
         y2 = ifelse(next2 == "unk",ifelse(init2=="yr",2.5,0.5),
                     ifelse(init2=="sub",ifelse(next2=="yr",1.9,1),2))) %>%
  mutate(lx1 = loc_id -1,
         lx2 = ifelse(next2 == "unk",loc_id-1.1,ifelse(init2==next2,loc_id,loc_id - 0.3)),
         ly1 = ifelse(init2 == "sub",1.0,2),
         ly2 = ifelse(next2 == "unk",ifelse(init2=="yr",2.65,0.75),
                      ifelse(init2=="sub",ifelse(next2=="yr",2.15,1.05),2))) %>%
  mutate(rot = ifelse(next2 == "unk",90,ifelse(next2 == init2,0,60)),
         angle = atan2(abs(y2 - y1), abs(x2 - x1)) * 180 / pi) %>%
  filter(loc_id>1,
         value !=0) %>%
  mutate(label = paste0(round(value,2)))# %>%
  # filter(j%in%c("Release","First_Trap","WEN")) #%>%
  # filter(init2 == "sub")

# Plot
p <- g %>%
  ggplot(aes(x = x1, y = y1)) +
  facet_wrap(~ReleaseSite, ncol = 1) +
  xlim(0, 4) +
  scale_x_continuous(
    breaks = seq_along(locs),
    labels = locs,
    expand = expansion(mult = c(0.05, 0.2))
    # breaks = seq_along(locs),
    # labels = locs
  ) +
  ylab("State") +
  xlab("Detection location") +
  ylim(0,2.5) +
  scale_y_continuous(
    breaks = c(0.5, 1, 2, 2.5),
    labels = c("Dead", "Subyearling", "Yearling", "Dead"),
    expand = expansion(mult = c(0.05, 0.2))

  ) +

  # Add arrows for transitions
  geom_segment(aes(xend = x2, yend = y2),
               size = 0.1,
               linejoin = "mitre",
               color = "grey",
               arrow = arrow(type = "closed")) +

  # Add dynamic text labels along the arrows
  geom_text(aes(x = (lx1 + lx2) / 2,
                y = (ly1 + ly2) / 2,
                label = label,
                angle   = angle),  # Use dynamic angle for rotation
            size = 6, hjust = 0.5, vjust = -0.5) +

  theme_classic() +
  theme(
    strip.text = element_text(size = 22),
    text = element_text(size = 22)# Adjust font size here
  )

# Print and save the plot
print(p)
# ggsave(p, filename = "output/HMM_plot_wo_barge.png", width = 10, height = 5)
