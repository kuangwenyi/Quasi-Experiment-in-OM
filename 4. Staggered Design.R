
# This is Section 4.1.4 and Section 4.2.4 Replication of Wang and Overby 2022 #

# Note 1: ggplot code is primarily based on 
#         1) Callaway's DID quick guide: https://bcallaway11.github.io/did/articles/TWFE.html
#         2) code from Baker et al. 2022. 
#         We thank Callaway and Baker et al. for making the plotting process easier. 
# Note 2: Replication studies of Wang and Overby 2022
#         Data can be downloaded from Mangement Science online supplements.



#### 1. Data and Variable Preparation -------------
# read in original files from Wang and Overyby 2022 
allstate = read.csv("all_state.csv", stringsAsFactors = F)
ninestate = read.csv("nine_state.csv", stringsAsFactors = F)

# create First Treatment indicator for analysis
indicator= c(1,1+which(diff(allstate$lending_club_available)==1))
test = allstate[indicator,]
test = test[test$lending_club_available!=0,]

firstT = test[,c("county_code", "new_timeperiod")]
names(firstT) = c("county_code", "FirstTreat")
unique(firstT$FirstTreat)

# bind to original data
allstate_n = allstate %>% left_join(firstT, by = "county_code")
ninestate_n = ninestate%>% left_join(firstT, by = "county_code")

# deal with control group. There is no treatment hence assigned 0 baased on Callways and Sant'Anna 2021
allstate_n$FirstTreat[is.na(allstate_n$FirstTreat)] = 0
ninestate_n$FirstTreat[is.na(ninestate_n$FirstTreat)] = 0

# check treatment occasions
unique(allstate_n$FirstTreat)
unique(ninestate_n$FirstTreat)

# weight variables based on Wang and Overby 2022 Table 4 and Table 5 Stata Code
allstate_n$bankruptcy_per_capita_aw = allstate_n$bankruptcy_per_capita*sqrt(allstate_n$cem_weights)
allstate_n$log_bankruptcy_aw = allstate_n$log_bankruptcy*sqrt(allstate_n$cem_weights)
allstate_n$population_estimation_aw = allstate_n$population_estimation*sqrt(allstate_n$cem_weights)
allstate_n$employed_individuals_aw = allstate_n$employed_individuals*sqrt(allstate_n$cem_weights)
allstate_n$monthly_earnings_aw = allstate_n$monthly_earnings*sqrt(allstate_n$cem_weights)
allstate_n$medianhouseholdincome_aw = allstate_n$medianhouseholdincome*sqrt(allstate_n$cem_weights)
allstate_n$labor_force_aw = allstate_n$labor_force*sqrt(allstate_n$cem_weights)

ninestate_n$bankruptcy_per_capita_aw = ninestate_n$bankruptcy_per_capita*sqrt(ninestate_n$cem_weights)
ninestate_n$log_bankruptcy_aw = ninestate_n$log_bankruptcy*sqrt(ninestate_n$cem_weights)
ninestate_n$population_estimation_aw = ninestate_n$population_estimation*sqrt(ninestate_n$cem_weights)
ninestate_n$employed_individuals_aw = ninestate_n$employed_individuals*sqrt(ninestate_n$cem_weights)
ninestate_n$monthly_earnings_aw = ninestate_n$monthly_earnings*sqrt(ninestate_n$cem_weights)
ninestate_n$medianhouseholdincome_aw = ninestate_n$medianhouseholdincome*sqrt(ninestate_n$cem_weights)
ninestate_n$labor_force_aw = ninestate_n$labor_force*sqrt(ninestate_n$cem_weights)


#### 2. Estimate using Callaway and Sant'Anna 2021 Estimator -----------

##### 2.1 All States Estimation --------
# DV is per capital bankruptcy filing
cs_allstate <- att_gt(yname = "bankruptcy_per_capita_aw",
                      data = allstate_n,
                      tname = "new_timeperiod",
                      idname = "county_code",
                      gname = "FirstTreat",
                      xformla = ~population_estimation_aw +
                        employed_individuals_aw + monthly_earnings_aw + medianhouseholdincome_aw +
                        labor_force_aw, # controls
                      control_group = c("notyettreated","nevertreated"),
                      est_method = "reg",
                      print_details = FALSE,
                      bstrap = TRUE,
                      cband = TRUE,
                      base_period = "varyiing",
                      clustervars = "county_code",
                      panel = TRUE,
                      anticipation = 0,
                      allow_unbalanced_panel = TRUE)

cs_es_allstate <- aggte(cs_allstate, type="dynamic", min_e = -8,  
                        max_e = 8, na.rm = TRUE)
summary(cs_es_allstate)

# DV is logged bankruptcy filing (not reported in the manuscript)
cs_allstate_log <- att_gt(yname = "log_bankruptcy_aw",
                          data = allstate_n,
                          tname = "new_timeperiod",
                          idname = "county_code",
                          gname = "FirstTreat",
                          xformla = ~population_estimation_aw +
                            employed_individuals_aw + monthly_earnings_aw + medianhouseholdincome_aw +
                            labor_force_aw, # controls
                          control_group = c("notyettreated","nevertreated"),
                          est_method = "reg",
                          print_details = FALSE,
                          bstrap = TRUE,
                          cband = TRUE,
                          base_period = "varyiing",
                          clustervars = "county_code",
                          panel = TRUE,
                          anticipation = 0,
                          allow_unbalanced_panel = TRUE)

cs_es_allstate_log <- aggte(cs_allstate_log, type="dynamic", min_e = -8,  
                            max_e = 8, na.rm = TRUE)
summary(cs_es_allstate_log)


##### 2.2 Nine States Estimation -----------
# DV is per capital bankruptcy filing
cs_9state <- att_gt(yname = "bankruptcy_per_capita_aw",
                    data = ninestate_n,
                    tname = "new_timeperiod",
                    idname = "county_code",
                    gname = "FirstTreat",
                    xformla = ~population_estimation_aw +
                      employed_individuals_aw + monthly_earnings_aw + medianhouseholdincome_aw +
                      labor_force_aw,
                    control_group = c("notyettreated","nevertreated"),
                    est_method = "reg",
                    print_details = FALSE,
                    bstrap = TRUE,
                    cband = TRUE,
                    base_period = "varyiing",
                    clustervars = "county_code",
                    panel = TRUE,
                    anticipation = 0,
                    allow_unbalanced_panel = TRUE)

cs_es_9state <- aggte(cs_9state, type="dynamic", min_e = -8,  
                      max_e = 8, na.rm = TRUE)
summary(cs_es_9state)

# DV is logged bankruptcy filing (not reported in the manuscript)
cs_9state_log <- att_gt(yname = "log_bankruptcy_aw",
                        data = ninestate_n,
                        tname = "new_timeperiod",
                        idname = "county_code",
                        gname = "FirstTreat",
                        xformla = ~population_estimation_aw +
                          employed_individuals_aw + monthly_earnings_aw + medianhouseholdincome_aw +
                          labor_force_aw,
                        control_group = c("notyettreated","nevertreated"),
                        est_method = "reg",
                        print_details = FALSE,
                        bstrap = TRUE,
                        cband = TRUE,
                        base_period = "varyiing",
                        clustervars = "county_code",
                        panel = TRUE,
                        anticipation = 0,
                        allow_unbalanced_panel = TRUE)

cs_es_9state_log <- aggte(cs_9state_log, type="dynamic", min_e = -8,  
                          max_e = 8, na.rm = TRUE)
summary(cs_es_9state_log)


#### 3. Event Study Type of Plot using Callaway and Sant'Anna Estimator ----------
# set plot theme
theme_set(theme_clean() + theme(plot.background = element_blank(),
                                legend.background = element_blank()))

##### 3.1 All state plot ------------------
# event study type aggregation
cs_es_allstate_all <- aggte(cs_allstate, type="dynamic", min_e = -8,  
                            max_e = 30, na.rm = TRUE)
# get the coefficients
summary(cs_es_allstate_all)

# plot
allstate_plot_all <- tidy(cs_es_allstate_all) %>% 
  as_tibble() %>% 
  mutate(group = as.factor(case_when(
    event.time < 0 ~ 1,
    event.time >= 0 ~ 2
  ))) %>% 
  # plot
  ggplot(aes(x = event.time, y = estimate)) + 
  geom_point(aes(fill= factor(group)), shape = 21) + geom_line() + 
  scale_fill_manual(values = c("#993441", "#0029a5")) + 
  ggtitle("The Dynamic Effect of Lending Club on Bankruptcy Filing - All States") + 
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, 
                    color = factor(group)), 
                linetype = "longdash", show.legend = FALSE) + 
  scale_color_manual(values = c("#993441", "#0029a5"))+
  geom_hline(yintercept = 0,  linetype = "longdash", color = "gray") + 
  geom_vline(xintercept = 0,  linetype = "longdash", color = "gray") + 
  labs(y = "Treatment Effect", x = "Quarters Relative to Treatment",
       subtitle = "Callaway and Sant’Anna (2021) Estimator") + 
  scale_x_continuous(breaks = seq(-8, 23, by = 1)) + 
  scale_y_continuous(breaks = seq(-4, 4, by = 0.05)) + 
  theme(axis.title.y = element_text(hjust = 0.5, vjust = 0.5, angle = 90),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")

# save (This is Figure 13)
ggsave(allstate_plot_all, filename = 'Allstate_DynamicPlot_all.png', 
       dpi = 500, width = 7.2, height = 4.5)


##### 3.2 Nine State Plot --------------
# event study type aggregation
cs_es_9state_all <- aggte(cs_9state, type="dynamic", min_e = -8,  
                          max_e = 30, na.rm = TRUE)

# get the coefficients
summary(cs_es_9state_all)

# plot
ninestate_plot_all <- tidy(cs_es_9state_all) %>% 
  as_tibble() %>% 
  mutate(group = as.factor(case_when(
    event.time < 0 ~ 1,
    event.time >= 0 ~ 2
  ))) %>% 
  # plot
  ggplot(aes(x = event.time, y = estimate)) + 
  geom_point(aes(fill= factor(group)), shape = 21) + geom_line() + 
  scale_fill_manual(values = c("#993441", "#0029a5")) + 
  ggtitle("The Dynamic Effect of Lending Club on Bankruptcy Filing - Nine States") + 
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, 
                    color = factor(group)), 
                linetype = "longdash", show.legend = FALSE) + 
  scale_color_manual(values = c("#993441", "#0029a5"))+
  geom_hline(yintercept = 0,  linetype = "longdash", color = "gray") + 
  geom_vline(xintercept = 0,  linetype = "longdash", color = "gray") + 
  labs(y = "Treatment Effect", x = "Quarters Relative to Treatment",
       subtitle = "Callaway and Sant’Anna (2021) Estimator") + 
  scale_x_continuous(breaks = seq(-8, 20, by = 1)) + 
  scale_y_continuous(breaks = seq(-4, 4, by = 0.05)) + 
  theme(axis.title.y = element_text(hjust = 0.5, vjust = 0.5, angle = 90),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")

# save (This is Figure 14)
ggsave(ninestate_plot_20, filename = 'Ninestate_DynamicPlot_all.png', 
       dpi = 500, width = 7.2, height = 4.5)


#### 4. Bacon Decomposition-------------
##### 4.1 Calculate Decomposition -----------------------
# All States (takes a while to run)
bacon_out <- bacon(bankruptcy_per_capita~ lending_club_available ,
                   data = allstate_n,
                   id_var = "county_code",
                   time_var = "new_timeperiod") 

save(bacon_out, file = "bacon_out_WangOveryAllState.rda")

# Nine States (takes a while to run)
bacon_out_9 <- bacon( bankruptcy_per_capita~ lending_club_available ,
                      data = ninestate_n,
                      id_var = "county_code",
                      time_var = "new_timeperiod") 

save(bacon_out_9, file = "bacon_out_WangOverby_9.rda")

##### 4.2 Plot of Bacon Decomposition ------------------

###### 4.2.1 All States Bacon Decomposition Plot -------
# get the total weight for each group
total_weights <- bacon_out %>% 
  group_by(type) %>% 
  summarize(weight = sum(weight))
# get the weighted average within group
group_avg <- bacon_out %>% 
  group_by(type) %>% 
  summarize(avg = weighted.mean(estimate, weight),
            weights = sum(weight))


#### Make Bacon Decomposition Plots
# Group 1: early v late 
EvL <- bacon_out %>% 
  filter(type == "Earlier vs Later Treated") %>% 
  ggplot(aes(x = weight, y = estimate)) + 
  geom_point(size = 3, alpha = 1/2) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = group_avg$avg[1], color = "#014182", size = 1.5) + 
  labs(x = "Weight", y = expression(widehat(delta^'DD'))) + 
  ggtitle(paste0("Early vs Later Treated \n Total Weight = ", scales::percent(total_weights$weight[1]))) + 
  scale_y_continuous(limits = c(-.12, 0.12)) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.y = element_text(angle = 360, hjust = 0.5, vjust = 0.5))

# Group 2: late v early
LvE <- bacon_out %>% 
  filter(type == "Later vs Earlier Treated") %>% 
  ggplot(aes(x = weight, y = estimate)) + 
  geom_point(size = 3, alpha = 1/2) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = group_avg$avg[2], color = "#014182", size = 1.5) + 
  labs(x = "Weight", y = expression(widehat(delta^'DD'))) + 
  scale_y_continuous(limits = c(-.12, 0.12)) + 
  ggtitle(paste0("Later vs Early Treated \n Total Weight = ", scales::percent(total_weights$weight[2]))) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.y = element_text(angle = 360, hjust = 0.5, vjust = 0.5))

# Group 3: Treated VS Untreated
TvU <- bacon_out %>% 
  filter(type == "Treated vs Untreated") %>% 
  ggplot(aes(x = weight, y = estimate)) + 
  geom_point(size = 3, alpha = 1/2) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = group_avg$avg[[3]], color = "#014182", size = 1.5) + 
  labs(x = "Weight", y = expression(widehat(delta^'DD'))) + 
  scale_y_continuous(limits = c(-.12, 0.12)) + 
  ggtitle(paste0("Treated vs Untreated \n Total Weight = ", scales::percent(total_weights$weight[3]))) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.y = element_text(angle = 360, hjust = 0.5, vjust = 0.5))

# combine the figures and save
BLL_decomp_plot <- EvL + LvE + TvU
plot(BLL_decomp_plot) # check plot before saving

# save (This is Figure 11)
ggsave(BLL_decomp_plot, filename = "AllState_decomp_plot.png", 
       dpi = 500, width = 10, height = 4)


###### 4.2.2 Nine States Bacon Decomposition Plot ----------
# get the total weight for each group. 
total_weights_9 <- bacon_out_9 %>% 
  group_by(type) %>% 
  summarize(weight = sum(weight))
# get the weighted average within group
group_avg_9 <- bacon_out_9 %>% 
  group_by(type) %>% 
  summarize(avg = weighted.mean(estimate, weight),
            weights = sum(weight))

# make the table
BLL_decomp_9 <- group_avg_9 %>% 
  kable("latex", digits = 3, align = 'lcc',
        booktabs = T,
        col.names = c("Type", "Weighted \n Average", "Total \n Weight"),
        label = "BLL_decomp") %>% 
  kable_styling(position = "center", font_size = 8,
                latex_options = c("HOLD_position", "scale_down"))

# save
write_lines(BLL_decomp_9, file =  "BLL_decomp_NineState.tex")


#### Make Bacon Decomp Figures 
# Group 1: first early v late
bacon_out_9
EvL_9 <- bacon_out_9 %>% 
  filter(type == "Earlier vs Later Treated") %>% 
  ggplot(aes(x = weight, y = estimate)) + 
  geom_point(size = 3, alpha = 1/2) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = group_avg_9$avg[1], color = "#014182", size = 1.5) + 
  labs(x = "Weight", y = expression(widehat(delta^'DD'))) + 
  ggtitle(paste0("Early vs Later Treated \n Total Weight = ", scales::percent(total_weights_9$weight[1]))) + 
  scale_y_continuous(limits = c(-.12, 0.15)) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.y = element_text(angle = 360, hjust = 0.5, vjust = 0.5))

# Group2: late v early
LvE_9 <- bacon_out_9 %>% 
  filter(type == "Later vs Earlier Treated") %>% 
  ggplot(aes(x = weight, y = estimate)) + 
  geom_point(size = 3, alpha = 1/2) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = group_avg_9$avg[2], color = "#014182", size = 1.5) + 
  labs(x = "Weight", y = expression(widehat(delta^'DD'))) + 
  scale_y_continuous(limits = c(-.12, 0.12)) + 
  ggtitle(paste0("Later vs Early Treated \n Total Weight = ", scales::percent(total_weights_9$weight[2]))) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.y = element_text(angle = 360, hjust = 0.5, vjust = 0.5))

# Group 3: Treated VS Untreated
TvU_9 <- bacon_out_9 %>% 
  filter(type == "Treated vs Untreated") %>% 
  ggplot(aes(x = weight, y = estimate)) + 
  geom_point(size = 3, alpha = 1/2) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = group_avg_9$avg[[3]], color = "#014182", size = 1.5) + 
  labs(x = "Weight", y = expression(widehat(delta^'DD'))) + 
  scale_y_continuous(limits = c(-.12, 0.12)) + 
  ggtitle(paste0("Treated vs Untreated \n Total Weight = ", scales::percent(total_weights_9$weight[3]))) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.y = element_text(angle = 360, hjust = 0.5, vjust = 0.5))

# combine the figures and save
BLL_decomp_plot_9 <- EvL_9 + LvE_9 + TvU_9
plot(BLL_decomp_plot_9) # check plots before saving

# save (This is Figure 12)
ggsave(BLL_decomp_plot_9, filename = "NineState_decomp_plot.png", 
       dpi = 500, width = 10, height = 4)


#### 5. Treatment Timing Plot --------------- 
##### 5.1 All states ---------
allstate_timing <- allstate_n %>% 
  select(county_code, new_timeperiod, FirstTreat) %>% 
  mutate(county_code_n = as.factor(county_code)) %>%
  mutate(County = fct_reorder(county_code_n, rank(desc(county_code_n)))) %>% 
  mutate(post = if_else(new_timeperiod < FirstTreat, "Pre", "Post")) %>% 
  mutate(post = factor(post, levels = c("Pre", "Post"))) %>% 
  ggplot(aes(x = new_timeperiod, y = County)) + 
  geom_tile(aes(fill = as.factor(post)), alpha = 3/4) + 
  scale_fill_manual(values = c("#5e001f", "#030E4F")) + 
  ggtitle("All States Lending Club Adoption")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom',
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.background = element_rect(color = "white")) 

plot(allstate_timing) # check plot

##### 5.2 Nine States --------
ninestate_timing <- ninestate_n %>% 
  select(county_code, new_timeperiod, FirstTreat) %>% 
  mutate(county_code_n = as.factor(county_code)) %>%
  mutate(County = fct_reorder(county_code_n, rank(desc(county_code_n)))) %>% 
  mutate(post = if_else(new_timeperiod < FirstTreat, "Pre", "Post")) %>% 
  mutate(post = factor(post, levels = c("Pre", "Post"))) %>% 
  ggplot(aes(x = new_timeperiod, y = County)) + 
  geom_tile(aes(fill = as.factor(post)), alpha = 3/4) + 
  scale_fill_manual(values = c("#5e001f", "#030E4F")) +
  ggtitle("Nine States Lending Club Adoption")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom',
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.background = element_rect(color = "white")) 

plot(ninestate_timing) # check plot

# Combine the two timing plots and save (This is Figure 10)
FTotalTiming = allstate_timing + ninestate_timing

ggsave(FTotalTiming, filename = "AllTiming_All+9.png", 
       dpi = 500, width = 12, height = 7)


