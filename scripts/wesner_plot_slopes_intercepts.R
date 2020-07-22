# libraries
library(tidyverse)
library(brms)
library(viridis)
library(ggridges)
library(gridExtra)
library(janitor)
library(tidyselect)
library(tidybayes)
library(nlme)
library(ggdist)
library(cowplot)

# load all NEON data
dat_all <- readRDS(file = "data/dat_all.rds")
data_conditioned_on <- dat_all %>% filter(year != 2019)
data_new <- dat_all %>% filter(year == 2019)
dat <- dat_all %>% filter(year != 2019)

# load model
mod <- readRDS("data/dd_as_r_slope_6-30-20.RDS") # r slopes after talking with JSW 6/30

# mod <- brm(data = dat,
#     log_count_corrected ~
#       log_mids_center * log10(degree_days) +
#       (1+ log10(degree_days)*log_mids_center|siteID) + (1|year), 
#     family = gaussian(),
#     prior = c(prior(normal(-1, 1),
#                     class = "b", coef = "log_mids_center"),
#               prior(normal(0,1), class = "b"),
#               prior(normal(4.5, 1.5), class = "Intercept"),
#               prior(cauchy(0, 1), class = "sigma"),
#               prior(cauchy(0,1), class = "sd")),
#     chains = 4, 
#     iter = 4000,
#     cores = 4)
# 
# 
# 
# 

# SLOPES and INTERCEPTS ------------------------------------------------------------------

# site and degree day slopes
sitedd_a <- posterior_samples(mod) %>% clean_names() %>% mutate(iter = 1:nrow(.)) %>% 
  select(!contains("cor_site")) %>% 
  pivot_longer(contains("r_site"), names_to = "ranef", values_to = "offset") %>%
  mutate(site = str_sub(ranef, 11,14),
         ranef = paste0("r_",str_sub(ranef, 16))) %>%
  pivot_wider(names_from = ranef, values_from = offset) %>%
  mutate(site = toupper(site)) %>% glimpse() %>% 
  left_join(data_conditioned_on %>% select(siteID, degree_days, year) %>% rename(site = siteID) %>% distinct(site, degree_days, year)) %>%
  mutate(Slope = b_log_mids_center + r_log_mids_center + (b_log_mids_center_log10degree_days + 
                                                            r_log10degree_days_log_mids_center)*log10(degree_days),
         Intercept = b_intercept + r_intercept + 
           (b_log_mids_center + r_log_mids_center)*mean(dat$log_mids_center) +
           (b_log10degree_days + r_log10degree_days)*log10(degree_days) + 
           (b_log_mids_center_log10degree_days + r_log10degree_days_log_mids_center)*log10(degree_days)*mean(dat$log_mids_center)) %>%
  mutate(Intercept = case_when(year == "2017" ~ Intercept + r_year_2017_intercept, 
                               TRUE ~ Intercept + r_year_2018_intercept),
         group = paste0(site,round(degree_days,0)),
         prediction_level = "Sites in original model") %>% 
  glimpse()

# posterior slopes for sites with data only in 2019
sitedd_b <- sitedd_a %>% filter(site %in% c("WLOU", "WALK", "REDB")) %>% #just chose 3 existing sites as placeholders
  mutate(site = case_when(site == "WLOU" ~ "SYCA",
                          site == "WALK" ~ "BLUE",
                          site == "REDB" ~ "BLDE"),
         siteID = site) %>%  
  left_join(data_new %>% select(siteID, degree_days) %>% rename(site = siteID) %>% distinct(site, degree_days)) %>%
  mutate(Slope= b_log_mids_center + rnorm(nrow(.), 0, sd_site_id_log_mids_center) +
           (b_log_mids_center_log10degree_days + rnorm(nrow(.), 0, sd_site_id_log10degree_days_log_mids_center))*log10(degree_days),
         Intercept = b_intercept + rnorm(nrow(.), 0, sd_site_id_intercept) + rnorm(nrow(.), 0, sd_year_intercept) +
           (b_log10degree_days + rnorm(nrow(.), 0, sd_site_id_log10degree_days))*log10(degree_days) + 
           (b_log_mids_center + rnorm(nrow(.), 0, sd_site_id_log_mids_center))*mean(dat$log_mids_center) + 
           (b_log_mids_center_log10degree_days + rnorm(nrow(.), 0, sd_site_id_log10degree_days_log_mids_center))*mean(dat$log_mids_center)*log10(degree_days)) %>%
  mutate(group = paste0(site,round(degree_days,0)),
         prediction_level = "Sites not in original model") %>% 
  glimpse()

# combine predictions
sitedd.all <- bind_rows(sitedd_a, sitedd_b) %>% 
  mutate(order = case_when(prediction_level == "Sites not in original model" ~ -100,
                           TRUE ~ Slope),
         degree_days = round(degree_days, 0),
         siteID = site) %>% 
  select(siteID, Slope, Intercept, degree_days, prediction_level) %>%
  group_by(siteID) %>% 
  mutate(median = median(Slope),
         median = case_when(siteID %in% c("SYCA", "BLDE", "BLUE") ~ 0,  #place 2019-only sites at the bottom
                            TRUE ~ median),
         color = case_when(siteID %in% c("SYCA", "BLDE", "BLUE") ~ "same",
                           TRUE ~ siteID)) %>% 
  #filter(siteID == "ARIK" | siteID == "WLOU"| siteID == "LECO") %>%
  gather(Parameter, value, c("Slope", "Intercept")) %>% 
  mutate(Parameter = fct_relevel(Parameter, "Slope")) %>% 
  ungroup() %>% 
  glimpse()

saveRDS(sitedd.all, file = "data/sitedd.all.rds")

# Plot
plot_slopes_wesner <- sitedd.all %>% 
  filter(Parameter == "Slope") %>% 
  ggplot() +
  geom_density_ridges_gradient(aes(x = value, y = reorder(siteID, -median), 
                                   fill = degree_days, group = interaction(degree_days, siteID)),
                               scale = 2,
                               rel_min_height = 0.01,
                               quantile_lines = TRUE, quantiles = 2,
                               alpha = 0.5) +
  # facet_wrap(~Parameter, scales = "free_x") +
  xlim(c(-1.8, -0.9)) +
  scale_fill_viridis(alpha = 0.5) +
  # geom_vline(aes(xintercept = -1.56),
  #            color = "black",
  #            linetype = "dashed",
  #            size = 1) +
  # geom_vline(aes(xintercept = -1.07),
  #            color = "black",
  #            linetype = "dashed",
  #            size = 1) +
  theme_bw() +
  # guides(fill = F) +
  labs(y = "Site",
       x = "Slope",
       subtitle = "") +
  NULL

plot_intercepts_wesner <- sitedd.all %>% 
  filter(Parameter != "Slope") %>% 
  ggplot() +
  geom_density_ridges_gradient(aes(x = value, y = reorder(siteID, -median), 
                                   fill = degree_days, group = interaction(degree_days, siteID)),
                               scale = 2,
                               rel_min_height = 0.01,
                               quantile_lines = TRUE, quantiles = 2,
                               alpha = 0.5) +
  # facet_wrap(~Parameter, scales = "free_x") +
  xlim(c(2.5, 6)) +
  scale_fill_viridis(alpha = 0.5) +
  # geom_vline(aes(xintercept = -1.56),
  #            color = "black",
  #            linetype = "dashed",
  #            size = 1) +
  # geom_vline(aes(xintercept = -1.07),
  #            color = "black",
  #            linetype = "dashed",
  #            size = 1) +
  theme_bw() +
  # guides(fill = F) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  labs(y = "Site",
       x = "Intercept (Community Height, log10N')",
       subtitle = "") +
  NULL


# combine plots 
p <- plot_grid(plot_slopes_wesner+ theme(legend.position = "none") ,
                                             plot_intercepts_wesner + theme(legend.position = "none"),
               align = "h", rel_widths = c(1, 0.8),
               labels = "auto", hjust = c(-5.9, -0.5))

# get legend
legend <- get_legend(
  # create some space to the left of the legend
  plot_intercepts_wesner + theme(legend.box.margin = margin(0, 0, 0, 12)))

# add legend
plot_slopeandintercept <- plot_grid(p, legend,
                                    rel_widths = c(3, 0.7))

plot_slopeandintercept 

ggsave(plot_slopeandintercept , file = "plots/plot_slopeandintercept.png", dpi = 500, width = 6, height = 8)


# summary stats
sitedd.all %>% 
  group_by(siteID, Parameter, degree_days) %>% 
  median_qi(value) %>% 
  pivot_wider(names_from = Parameter, values_from = c(value, .lower, .upper)) %>% 
  arrange(value_Slope) 



# Plot degree days on x and slope and comm height on y --------------------

# extract fitted values of log_count_corrected along the range of degree days - assumes log_mids_center is at it's mean
c_eff_int <- conditional_effects(mod, effects = c("degree_days"), re_formula = NA, method = "fitted")

# convert c_eff to list
c_eff_list <- plot(c_eff, plot = FALSE)

# convert c_eff_list to tibble for plotting
c_eff_toplot <- c_eff_list$degree_days$data %>% bind_rows() %>% clean_names() %>% 
  select(degree_days, estimate, lower, upper) %>% 
  mutate(Parameter = "Community Height")



# do same as above, but for slopes
min_dd <- min(dat$degree_days)
max_dd <- max(dat$degree_days)

c_eff_slopes <- posterior_samples(mod) %>% as_tibble() %>% clean_names() %>% mutate(iter = 1:nrow(.)) %>% 
  expand_grid(degree_days = seq(min_dd, max_dd, length.out = 100)) %>% 
  mutate(Slope = b_log_mids_center + b_log_mids_center_log10degree_days*log10(degree_days)) %>% 
  select(degree_days, Slope, iter) %>% glimpse() %>% 
  group_by(degree_days) %>% 
  median_qi(Slope) %>% 
  rename(estimate = Slope, lower = .lower, upper = .upper) %>% 
  mutate(Parameter = "Slope")

#combine 
slope_int = c_eff_toplot %>% bind_rows(c_eff_slopes) %>% glimpse()

# extract mean and CrI for each site*dd combo
sitedd_posts <- sitedd.all %>% 
  # filter(prediction_level == "Sites in original model") %>% 
  group_by(siteID, Parameter, degree_days, prediction_level) %>% 
  median_qi(value) %>% 
  pivot_wider(names_from = Parameter, values_from = c(value, .lower, .upper)) %>% 
  arrange(value_Slope) 

# plot slope versus degree_days
plot_slope_dd <- slope_int %>% filter(Parameter == "Slope") %>% clean_names() %>%
  ggplot() + 
  geom_line(aes(x = degree_days, y = estimate)) + 
  # geom_label(data = sitedd_posts, aes(x = degree_days, y = value_Intercept, label = siteID)) +
  geom_pointrange(data = sitedd_posts %>% filter(prediction_level == "Sites in original model"),
               aes(x = degree_days, y = value_Slope, ymin = .lower_Slope, ymax = .upper_Slope,
                   group = degree_days), size = 0.2) +
  geom_ribbon(alpha = 0.2, aes(x = degree_days, ymin = lower, ymax = upper)) +
  scale_x_log10() +
  theme_bw() +
  theme_bw() +
  theme(plot.margin = margin(0, 30, 0, 0)) +
  ylim(-1.75, -0.75) +
  labs(y = "Slope",
       x = "Degree Days") +
  # coord_cartesian(ylim = c(3, 5.5)) +
  NULL

# plot change of community height as a function of degree_days
plot_int_dd <- slope_int %>% filter(Parameter == "Community Height") %>% clean_names() %>%
  ggplot() + 
  geom_line(aes(x = degree_days, y = estimate)) + 
  # geom_label(data = sitedd_posts, aes(x = degree_days, y = value_Intercept, label = siteID)) +
  geom_pointrange(data = sitedd_posts %>% filter(prediction_level == "Sites in original model"),
                  aes(x = degree_days, y = value_Intercept, ymin = .lower_Intercept, ymax = .upper_Intercept,
                      group = degree_days), size = 0.2) +
  geom_ribbon(alpha = 0.2, aes(x = degree_days, ymin = lower, ymax = upper)) +
  scale_x_log10() +
  theme_bw() +
  theme(plot.margin = margin(0, 30, 0, 0)) +
  ylim(2.5, 5.5) +
  labs(y = "Intercept (Community Height)\nlog10N'",
       x = "Degree Days") +
  # coord_cartesian(ylim = c(3, 5.5)) +
  NULL


plot_slopeint_dd <- plot_grid(plot_slope_dd, plot_int_dd, labels = "auto", hjust = c(-5.7, -5.))

plot_slopeint_dd
ggsave(plot_slopeint_dd, file = "plots/plot_slope_dd.png", dpi = 600, 
       width = 7, height = 2.5)

