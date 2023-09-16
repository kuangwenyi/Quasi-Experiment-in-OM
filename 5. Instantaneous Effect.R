
# This is Section 5 Instantaneous Effect Replication of Avramidis et al. 2022 #

# Note 1: ggplot code is primarily based on 
#         1) Callaway's DID quick guide: https://bcallaway11.github.io/did/articles/TWFE.html
#         2) code from Baker et al. 2022. 
#         We thank Callaway and Baker et al. for making the plotting process easier. 
# Note 2: Replication studies of Avramidis et al. 2022
#         Data can be downloaded from Mangement Science online supplements.

# data used by Avramidis et al. is called master.dta
avra = haven::read_dta("master.dta")
avra = avra[avra$sample_flag ==1,] # only keep sample_flag ==1 based on the source code of Avramidis et al.2022

# set up variables based on Avramidis et al.2022
Y = "wdif_num_branches_bm"
G = "zip3code"
T = "year"
D = "bankmerger"

# create control variables
control_1 = c("demand_shift","i2009","i2010","i2011", "i2012", 
              "i2013", "i2014", "i2015", "i2016") # column 2 in Table 4 (Avramidis et al.2022)
control_2 = c("demand_shift","pdiflag_wnum_req_all","i2009","i2010","i2011", "i2012", 
              "i2013", "i2014", "i2015", "i2016") # column 3 in Table 4 (Avramidis et al.2022)

# decompose weights
twowayfeweights(avra, Y, G, T, D, cmd_type = "feTR") # without controls
twowayfeweights(avra, Y, G, T, D, cmd_type = "feTR", controls = control_1) # with controls

# estimate instantaneous effect
# Following are the initial Stata codes for models from Column 2 to Column 4 in Avramidis et al.2022
# xtreg wdif_num_branches_bm bankmerger demand_shift i.year if sample_flag==1,fe robust
# xtreg wdif_num_branches_bm bankmerger demand_shift pdiflag_wnum_req_all i.year if sample_flag==1,fe robust

# Column 2 in Table 4 (Avramidis et al.2022)
set.seed(1234)
es_col2 = did_multiplegt(avra, Y, G, T, D, placebo = 4, dynamic = 0, controls = control_1, 
                         brep = 2, cluster = "zip3code") # standard errors clustered as was in Avramidis et al.2022

# Column 3 in Table 4 (Avramidis et al.2022)
es_col3 = did_multiplegt(avra, Y, G, T, D, placebo = 4, dynamic =0, controls = control_2,
                         brep = 2, cluster = "zip3code")# standard errors clustered as was in Avramidis et al.2022

# unlist and combine model outputs for plotting
total = list(es_col2, es_col3)
test = as.data.frame(do.call(cbind, total))
test$Names = rownames(test)

# extract estimates and se for plot
estimates =rbind(test[grepl('^placebo_', test$Names), ], test[grepl('^effect', test$Names), ])
se = rbind(test[grepl('^se_placebo_', test$Names), ], test[grepl('^se_effect', test$Names), ])

estimates$Names = NULL
se$Names = NULL

# conver to numeric
for (i in 1:ncol(estimates)){
  estimates[[i]] <- as.numeric(estimates[[i]])}

for (j in 1:ncol(se)){
  se[[j]] <- as.numeric(se[[j]])}


# Column 2 Plot
df <- data.frame(x =-4:0,
                 F =estimates$V1,
                 L =estimates$V1 - se$V1*1.96,
                 U =estimates$V1 + se$V1*1.96)

plot_col2 =  df %>% 
  mutate(group = as.factor(case_when(
    x < 0~1,
    x >= 0~2 ))) %>% 
  ggplot(aes(x = x, y = F)) + 
  geom_point(aes(fill= factor(group)), shape = 21) + geom_line() + 
  scale_fill_manual(values = c("#993441", "#0029a5")) + 
  ggtitle("Bank Branch Changes") + 
  geom_errorbar(aes(ymin = L, ymax = U, 
                    color = factor(group)), 
                linetype = "longdash", show.legend = FALSE) + 
  scale_color_manual(values = c("#993441", "#0029a5"))+
  geom_hline(yintercept = 0,  linetype = "longdash", color = "gray") + 
  geom_vline(xintercept = 0,  linetype = "longdash", color = "gray") + 
  labs(y = "Treatment Effect", x = "Time Since Treatment",
       subtitle = "Column 2 of Table 4 in Avramidis et al.2022") + 
  scale_x_continuous(breaks = seq(-4, 0, by = 1)) + 
  scale_y_continuous(breaks = seq(-4, 1.5, by = 1)) + 
  theme(axis.title.y = element_text(hjust = 0.5, vjust = 0.5, angle = 90),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")

# Column 3 Plot
df <- data.frame(x =-4:0,
                 F =estimates$V2,
                 L =estimates$V2 - se$V2*1.96,
                 U =estimates$V2 + se$V2*1.96)

plot_col3 =  df %>% 
  mutate(group = as.factor(case_when(
    x < 0~1,
    x >= 0~2 ))) %>% 
  ggplot(aes(x = x, y = F)) + 
  geom_point(aes(fill= factor(group)), shape = 21) + geom_line() + 
  scale_fill_manual(values = c("#993441", "#0029a5")) + 
  ggtitle("Bank Branch Changes") + 
  geom_errorbar(aes(ymin = L, ymax = U, 
                    color = factor(group)), 
                linetype = "longdash", show.legend = FALSE) + 
  scale_color_manual(values = c("#993441", "#0029a5"))+
  geom_hline(yintercept = 0,  linetype = "longdash", color = "gray") + 
  geom_vline(xintercept = 0,  linetype = "longdash", color = "gray") + 
  labs(y = "Treatment Effect", x = "Time Since Treatment",
       subtitle = "Column 3 of Table 4 in Avramidis et al.2022") + 
  scale_x_continuous(breaks = seq(-4, 0, by = 1)) + 
  scale_y_continuous(breaks = seq(-1.5,1.5, by = 0.5)) + 
  theme(axis.title.y = element_text(hjust = 0.5, vjust = 0.5, angle = 90),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")


# save all the plots (This is Figure 9 in the manuscript)
library(ggpubr)
totoalplot = ggarrange(plot_col2, plot_col3,
                       labels = c("A", "B"),
                       ncol = 2, nrow = 1)

ggsave(totoalplot, filename = "TotalPlot_Avra_Col23_8*3.5.png", 
       dpi = 500, width = 8, height = 3.5)

