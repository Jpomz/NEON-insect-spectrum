# 2a modeling filtered binned data

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
source("scripts/custom_functions.R")

# create a new folder to save outputs for MS
if(!dir.exists("ms")){
  dir.create("ms")  
}

# load data
# this data set has been filtered to remove the individuals which had damage that affected their measurements
dat <- readRDS("data/binned_no_damage_MN_degree_days_2017_2019.RDS")

# combine date columns into one collect date factor
dat <- dat %>%
  unite(collectDate, year:day, remove = FALSE)


# Subset data which has degree_days values and more than 200 daily observations per year
dat <- dat %>%
  filter(!is.na(degree_days), 
         yearly.n >= 200)

# Make dataframe with all years
dat_all <- dat

# filter out 2017-2018 data for bayesian model
dat <- dat %>%
  filter(year != 2019)

dat.2019 <- dat_all %>%
  filter(year == 2019)

# number of populated bins per sample?
dat %>%
  group_by(siteID, collectDate) %>%
  count() %>%
  pull(n) %>%
  range()


# plot data by siteID
# colorcoded by degree days 
ggplot(dat,
       aes(y = log_count_corrected, 
           x = log_mids_center)) +
  geom_point(aes(fill = log10(degree_days)),
             alpha = 0.7,
             size = 2,
             pch = 21, 
             color = "black") +
  facet_wrap(siteID~.,
             ncol = 3) +
  theme_bw() +
  scale_fill_viridis_c(option = "plasma") +
  labs(y = "Log10 Normalized Abundance (N)",
       x = "Log10 M (mg, dry weight) centered",
       title = "NEON Macroinvertebrate 2017-2018 M-N observations"
       ) 


# Bayesian model ----------------------------------------------------------

get_prior(data = dat,
    log_count_corrected ~
      log_mids_center * log10(degree_days) +
      (1 + log10(degree_days)|siteID) + (1|year), 
    family = gaussian())

# define my priors
my_priors <- c(prior(normal(-2, 2), class = "b", coef = "log_mids_center"),
               prior(normal(0,1), class = "b"),
               prior(normal(4.25, 2), class = "Intercept"),
               prior(exponential(1), class = "sigma"),
               prior(exponential(1), class = "sd"))

# Bayesian Model
mod <- brm(data = dat,
           log_count_corrected ~
             log_mids_center * log10(degree_days) +
             (1 + log_mids_center*log10(degree_days)|siteID) + (1|year), 
           family = gaussian(),
           prior = my_priors,
           chains = 4, 
           iter = 4000,
           cores = 4)
print(mod)

# prior vs posterior plot ####
# re do this with posterior draws instead of rnorm. 
data.frame(value = c(
  # Int_prior = 
  rnorm(1000, 4.25, 2),
  # Int_post = 
  rnorm(1000, fixef(mod)[1,1], fixef(mod)[1,2]),
  # M_prior = 
  rnorm(1000, -2, 2),
  # M_post =
  rnorm(1000, fixef(mod)[2,1], fixef(mod)[2,2]),
  # b_prior =
  rnorm(1000, 0, 1),
  # DD_post = 
  rnorm(1000, fixef(mod)[3,1], fixef(mod)[3,2]),
  # b_prior =
  rnorm(1000, 0, 1),
  # M_DD_post =
  rnorm(1000, fixef(mod)[4,1], fixef(mod)[4,2])),
  b_name = rep(c("Intercept", "Intercept", "M", "M",
                 "Degree_days", "Degree_days", "M:DD", "M:DD"),
               each = 1000),
  sample = rep(c("prior", "post",
                 "prior", "post",
                 "prior", "post",
                 "prior", "post"), each = 1000)) %>%
  ggplot() +
  geom_density(aes(x = value,
                   fill = sample),
               alpha = 0.5) +
  facet_wrap(.~b_name,
             #cols = 2,
             scales = "free") +
  theme_bw() +
  labs(title =
         "Prior and posterior distributions of parameter coefficients")
ggsave("ms/SI_fig1.png")

pp_check(mod, type = "boxplot")
ggsave("ms/SI_fig2.png")

plot(conditional_effects(mod), points = TRUE)

# saveRDS(mod, "data/dd_as_r_slope_6-30-20.RDS") # r slopes after talking with JSW 6/30

# mod <- readRDS("data/dd_as_r_slope_6-30-20.RDS") # r slopes after talking with JSW 6/30


# slope beta distributions ------------------------------------------------

# plotting slopes across site:degree_days
mod.slopes <- as_tibble(
  posterior_samples(
    mod, pars = c("log_mids_center", "log10degree_days")))
site_dd <- dat %>%
  select(siteID, degree_days) %>%
  unique()
# 
# # pivot longer to get site_M and site_M_DD offsets
site.slopes <- mod.slopes %>%
  pivot_longer(10:28, names_to = "M_key", values_to = "site_M") %>%
  mutate(siteID = str_sub(M_key, 10, 13)) %>%
  pivot_longer(10:28, names_to = "M_DD_key", values_to = "site_M_DD") %>%
  left_join(site_dd) %>%
  # not 100% certain that the below calculation is correct
  mutate(site_slope = b_log_mids_center + site_M,
    site_slope_dd = b_log_mids_center + site_M +
           (`b_log_mids_center:log10degree_days` + site_M_DD) *
           log10(degree_days)) %>%
  group_by(siteID, degree_days) %>%
  mutate(median.slope = median(site_slope)) %>%
  ungroup()

#### NEW CODE FROM WESNER - REPLACES the two code chunks above #####
# SITE LEVEL PREDICTIONS
site_slopes_a <- posterior_samples(mod) %>%
  clean_names() %>%
  mutate(iter = 1:nrow(.)) %>% 
  select(!contains("cor_site")) %>% 
  pivot_longer(contains("r_site"), names_to = "ranef", values_to = "offset") %>%
  mutate(site = str_sub(ranef, 11,14),
         ranef = paste0("r_",str_sub(ranef, 16))) %>% 
  pivot_wider(names_from = ranef, values_from = offset) %>%
  mutate(site = toupper(site)) %>% 
  # left_join(lms_combine) %>% 
  mutate(site_slope = b_log_mids_center + r_log_mids_center) %>%
  mutate(group = site,
         prediction_level = "Sites in original model") %>% 
  glimpse()

# posterior slopes for sites with data only in 2019
site_slopes_b <- site_slopes_a %>% filter(site %in% c("WLOU", "WALK", "REDB")) %>% #just chose 3 existing sites as placeholders
  # left_join(lms_combine) %>%
  mutate(site = case_when(site == "WLOU" ~ "SYCA",
                          site == "WALK" ~ "BLUE",
                          site == "REDB" ~ "BLDE"),
         site_slope= b_log_mids_center + rnorm(nrow(.), 0, sd_site_id_log_mids_center)) %>% 
mutate(group = site,
       model = "siterandom_year",
       prediction_level = "Sites not in original model") %>% 
  glimpse()

# combine predictions
site.slopes <- bind_rows(site_slopes_a, site_slopes_b) %>% 
  mutate(order =
           case_when(prediction_level ==
                       "Sites not in original model" ~ -100,
                           TRUE ~ site_slope),
         siteID = site)


#SITE X DEGREE DAY LEVEL PREDICTIONS
sitedd_slopes_a <- posterior_samples(mod) %>%
  clean_names() %>%
  mutate(iter = 1:nrow(.)) %>% 
  select(!contains("cor_site")) %>% 
  pivot_longer(contains("r_site"),
               names_to = "ranef",
               values_to = "offset") %>%
  mutate(site = str_sub(ranef, 11,14),
         ranef = paste0("r_",str_sub(ranef, 16))) %>% 
  pivot_wider(names_from = ranef, values_from = offset) %>%
  mutate(site = toupper(site)) %>% 
  # left_join(lms_combine) %>%
  # glimpse() %>% 
  mutate(site_slope_dd20172018 =
           b_log_mids_center + r_log_mids_center +
           (b_log_mids_center_log10degree_days +
              r_log10degree_days_log_mids_center)*
           log10(degree_days)) %>%
  mutate(group = paste0(site,round(degree_days,0)),
         prediction_level = "Sites in original model") %>% 
  glimpse()

# posterior slopes for sites with data only in 2019
sitedd_slopes_b <- site_slopes_a %>%
  filter(site %in% c("WLOU", "WALK", "REDB")) %>%
  #just chose 3 existing sites as placeholders
  mutate(site = case_when(site == "WLOU" ~ "SYCA",
                          site == "WALK" ~ "BLUE",
                          site == "REDB" ~ "BLDE"),
         siteID = site) %>%  
  left_join(lms_combine) %>%
  mutate(site_slope_dd= b_log_mids_center + rnorm(nrow(.), 0, sd_site_id_log_mids_center) +
           (b_log_mids_center_log10degree_days + rnorm(nrow(.), 0, sd_site_id_log10degree_days_log_mids_center))*log10(degree_days)) %>%
  mutate(group = paste0(site,round(degree_days,0)),
         prediction_level = "Sites not in original model") %>% 
  glimpse()

# combine predictions
site.slopes <- bind_rows(site_slopes_a, site_slopes_b) %>% 
  mutate(order =
           case_when(prediction_level ==
                       "Sites not in original model" ~ -100,
                     TRUE ~ site_slope),
         siteID = site)

sitedd.slopes <- bind_rows(sitedd_slopes_a, sitedd_slopes_b) %>% 
  mutate(order = case_when(prediction_level == "Sites not in original model" ~ -100,
                           TRUE ~ site_slope_dd),
         siteID = site)




# spatial variation ####
# 95% CrI for site-specific slopes
site.slopes %>%
  group_by(siteID) %>%
  summarize(l95 = quantile(site_slope, probs = 0.025),
            q50 = quantile(site_slope, probs = 0.5),
            u95 = quantile(site_slope, probs = 0.975)) %>%
  arrange(q50)

# range of site-specific slopes, ignoring DD
site.slopes %>%
  group_by(siteID) %>%
  summarize(l95 = quantile(site_slope, probs = 0.025),
            q50 = quantile(site_slope, probs = 0.5),
            u95 = quantile(site_slope, probs = 0.975)) %>%
  pull(q50) %>%
  range()

site.slopes %>%
  group_by(siteID) %>%
  summarize(l95 = quantile(site_slope, probs = 0.025),
            q50 = quantile(site_slope, probs = 0.5),
            u95 = quantile(site_slope, probs = 0.975)) %>%
  mutate(delta = u95 - l95) %>%
  ungroup() %>% 
  arrange(delta)

site.slopes %>%
  group_by(siteID, degree_days) %>%
  summarize(l95 = quantile(site_slope_dd, probs = 0.025),
            q50 = quantile(site_slope_dd, probs = 0.5),
            u95 = quantile(site_slope_dd, probs = 0.975)) %>%
  pull(q50) %>%
  range()

# difference in median slope parameters #####NOT SURE WHAT THIS IS? - WESNER ####
site.slopes %>%
  group_by(siteID) %>%
  summarize(med = median(site_slope)) %>%
  ungroup() %>%
  summarize(delta = max(med) -  min(med))
  

# 95% CrI for site-specific slopes, with degree days
# sitedd.slopes object above is not correct
# sitedd.slopes %>%
#   group_by(siteID, degree_days) %>%
#   summarize(l95 = quantile(site_slope_dd, probs = 0.025),
#             q50 = quantile(site_slope_dd, probs = 0.5),
#             u95 = quantile(site_slope_dd, probs = 0.975)) %>%
#   write.csv("ms/SI_site_dd_slopes.csv", row.names = FALSE)
# 
# # range of seasonal variation
# sitedd.slopes %>%
#   group_by(siteID, degree_days) %>%
#   summarize(l95 = quantile(site_slope_dd, probs = 0.025),
#             q50 = quantile(site_slope_dd, probs = 0.5),
#             u95 = quantile(site_slope_dd, probs = 0.975)) %>%
#   ungroup() %>%
#   mutate(delta = u95 - l95) %>%
#   arrange(delta)

# difference in seasonal median slope parameters within a site
site.slopes %>%
  group_by(siteID, degree_days) %>%
  summarize(med = median(site_slope_dd)) %>%
  ungroup() %>%
  group_by(siteID) %>%
  summarize(delta = max(med) -  min(med)) %>%
  arrange(delta)


# 95% CrI for median slopes
site.slopes %>%
  summarize(l95 = quantile(site_slope, probs = 0.025),
            q50 = quantile(site_slope, probs = 0.5),
            u95 = quantile(site_slope, probs = 0.975)) %>%
  arrange(q50)

#median site slope across latitudes?
site.info <- readr::read_csv("data/aquatic-field-sites.csv")

site.slopes %>%
  group_by(siteID) %>%
  summarize(q50 = quantile(site_slope, probs = 0.5)) %>%
  left_join(site.info[, c("Site ID", "Latitude")], by = c("siteID" = "Site ID")) %>%
  ggplot(aes(x = Latitude, y = q50)) +
  geom_point() +
  stat_smooth(method = "lm")
site.slopes %>%
  group_by(siteID) %>%
  summarize(q50 = quantile(site_slope, probs = 0.5)) %>%
  left_join(site.info[, c("Site ID", "Latitude")], by = c("siteID" = "Site ID")) %>%
  select(q50, Latitude) %>%
  cor()

  
# median slope across mean annual temperature?
site.slopes %>%
  group_by(siteID) %>%
  summarize(q50 = quantile(site_slope, probs = 0.5)) %>%
  left_join(site.info[, c("Site ID", "mat.c")], by = c("siteID" = "Site ID")) %>%
  ggplot(aes(x = mat.c, y = q50)) +
  geom_point() +
  stat_smooth(method = "lm")

site.slopes %>%
  group_by(siteID) %>%
  summarize(q50 = quantile(site_slope, probs = 0.5)) %>%
  left_join(site.info[, c("Site ID", "mat.c")], by = c("siteID" = "Site ID")) %>%
  select(q50, mat.c) %>%
  cor()

# plots for slope distributions across sites ####
# still a work in progress
site.slopes %>%
  ggplot(aes(y = fct_reorder(siteID, median.slope),
             x = site_slope)) +
  tidybayes::stat_pointinterval() +
  geom_vline(xintercept = -1.33,
             color = "red",
             linetype = "longdash") +
  geom_vline(xintercept = -1.61, 
             color = "red",
             linetype = "dashed")+
  geom_vline(xintercept = -1.05, 
             color = "red",
             linetype = "dashed") +
  theme_bw() +
  labs(y = "Site",
       x = "slope coefficient",
       title = "Site-specific slope coefficient distribution")
ggsave("ms/SI_fig6.png")

# Plot showing effects of degree_days within a site - WESNER NEW
# sitedd.slopes %>%
#   select(siteID, site_slope_dd, degree_days) %>%
#   #filter(siteID == "ARIK" | siteID == "WLOU"| siteID == "LECO") %>%
#   ggplot() +
#   geom_density_ridges_gradient(
#     aes(x = site_slope_dd,
#         y = interaction(degree_days, siteID),
#         group = interaction(siteID, degree_days),
#         fill = siteID),
#     scale = 2,
#     rel_min_height = 0.01,
#     quantile_lines = TRUE, quantiles = 2,
#     alpha = 0.5) +
#   xlim(c(-1.8, -0.9)) +
#   #scale_fill_viridis_d(alpha = 0.5) +
#   # geom_vline(aes(xintercept = -1.56),
#   #            color = "black",
#   #            linetype = "dashed",
#   #            size = 1) +
#   # geom_vline(aes(xintercept = -1.07),
#   #            color = "black",
#   #            linetype = "dashed",
#   #            size = 1) +
#   theme_bw() +
#   # facet_wrap(.~siteID,
#   #            scales = "free_y")+
#   labs(y = "Site",
#        x = "Slope coefficient",
#        title = "Effect of degree days on Site-specific \nslope coefficient distributions",
#        subtitle = "") +
#   NULL

# WESNER ADD 2019 DATA POSTERIOR PLOTS
# compute least squares regressions for each site/year/degree_day combo
dat_formod <- dat_all %>%
  mutate(group = paste0(siteID,round(degree_days,0), "_",year))

lms <- lmList(log_count_corrected ~ log_mids_center | group, data = dat_formod)

#get coefficients
lms_coefs <- coef(lms) %>% rownames_to_column() %>% clean_names() %>% 
  rename(group = rowname,
         int= intercept,
         slope = log_mids_center) %>% 
  mutate(year = str_sub(group,-4,-1),
         site = str_sub(group, 1, 4),
         degree_days = parse_number(str_sub(group, -50, -6))) %>% 
  glimpse()

#select columns
lms_combine <- lms_coefs %>%
  select(site, year, degree_days, slope) %>% 
  mutate(siteID = site)

# add posterior medians for plotting order
lms_toplot <- sitedd.slopes %>%
           select(siteID, site_slope_dd, degree_days) %>% 
           group_by(siteID) %>% 
           summarize(median = median(site_slope_dd)) %>% 
  left_join(lms_combine) %>% 
  mutate(fit_pred = case_when(year == "2019" ~ "Prediction",
                                  TRUE ~ "Fitted"))


plot_all_slopes_with2019data <- sitedd.slopes %>%
  select(siteID, site_slope_dd, degree_days) %>%
  group_by(siteID) %>% 
  left_join(lms_toplot %>% select(degree_days, year)) %>% 
  mutate(median = median(site_slope_dd),
         median = case_when(siteID %in% c("SYCA", "BLDE", "BLUE") ~ 0,  #place 2019-only sites at the bottom
                            TRUE ~ median),
         color = case_when(siteID %in% c("SYCA", "BLDE", "BLUE") ~ "same",
                           TRUE ~ siteID),
         fit_pred = case_when(year == "2019" ~ "Prediction",
                              TRUE ~ "Fitted")) %>% 
  #filter(siteID == "ARIK" | siteID == "WLOU"| siteID == "LECO") %>%
  ggplot() +
  geom_density_ridges_gradient(aes(x = site_slope_dd, y = reorder(as.character(degree_days), degree_days), 
                                   fill = reorder(color, median), group = interaction(degree_days, siteID)),
                               scale = 2,
                               rel_min_height = 0.01,
                               quantile_lines = TRUE, quantiles = 2,
                               alpha = 0.5) +
  geom_point(data = lms_toplot ,
             aes(x = slope, y = reorder(as.character(degree_days), degree_days),
                 shape = year,
                 size = year)) +
  facet_grid(reorder(siteID, median) ~ fit_pred, scales = "free_y", switch = "y") +
  xlim(c(-1.8, -0.9)) +
  scale_fill_viridis_d(alpha = 0.5) +
  # scale_color_manual(values = c("black", "black", "grey90")) +
  scale_shape_manual(values = c(16, 16, 4)) +
  scale_size_manual(values = c(1, 1, 2)) +
  # geom_vline(aes(xintercept = -1.56),
  #            color = "black",
  #            linetype = "dashed",
  #            size = 1) +
  # geom_vline(aes(xintercept = -1.07),
  #            color = "black",
  #            linetype = "dashed",
  #            size = 1) +
  theme_bw() +
  theme(panel.spacing.x=unit(0, "lines"),panel.spacing.y=unit(0, "lines"),
        axis.text.y = element_text(size = 5)) +
  guides(fill = F) +
  labs(y = "Site",
       x = "Slope coefficient",
       # title = "Effect of degree days on Site-specific \nslope coefficient distributions",
       subtitle = "") +
  NULL


plot_all_slopes_with2019data


ggsave(plot_all_slopes_with2019data,file = "plots/plot_all_slopeswith2019data.png",
       width = 6, height = 10 )


# model intercepts --------------------------------------------------------

mod.int <- as_tibble(
  posterior_samples(
    mod, pars = c("Intercept", "log10degree_days")))

# pivot longer to get site_int and site_int_DD offsets
site.int <- mod.int %>%
  pivot_longer(7:25, names_to = "int_key",
               values_to = "site_int_offset") %>%
  mutate(siteID = str_sub(int_key, 10, 13)) %>%
  pivot_longer(16:34, names_to = "int_DD_key",
               values_to = "site_int_DD_offset") %>%
  left_join(site_dd) %>%
  # not 100% certain that the below calculation is correct
  mutate(site_int = b_Intercept + site_int_offset +
           b_log10degree_days * log10(degree_days),
         site_int_dd = b_Intercept + site_int_offset +
           (site_int_DD_offset + b_log10degree_days) *
           log10(degree_days)) %>%
  group_by(siteID, degree_days) %>%
  mutate(median.int = median(site_int)) %>%
  ungroup()

#site specific intercepts, averaged over degree days
site.int %>%
  group_by(siteID) %>%
  summarize(l95 = quantile(site_int, probs = 0.025),
            q50 = quantile(site_int, probs = 0.5),
            u95 = quantile(site_int, probs = 0.975)) %>%
  arrange(q50)

#site specific intercepts including effects of degree days
site.int %>%
  group_by(siteID, degree_days) %>%
  summarize(l95 = quantile(site_int_dd, probs = 0.025),
            q50 = quantile(site_int_dd, probs = 0.5),
            u95 = quantile(site_int_dd, probs = 0.975)) %>%
  arrange(q50)

site.int %>%
  ggplot(aes(y = fct_reorder(siteID, median.int),
             x = site_int)) +
  ggdist::stat_pointinterval() +
  # geom_vline(xintercept = -1.33,
  #            color = "red",
  #            linetype = "longdash") +
  # geom_vline(xintercept = -1.61, 
  #            color = "red",
  #            linetype = "dashed")+
  # geom_vline(xintercept = -1.05, 
  #            color = "red",
  #            linetype = "dashed") +
  theme_bw() +
  labs(y = "Site",
       x = "Intercept coefficient",
       title = "Site-specific Intercept coefficient distribution")
ggsave("ms/site-intercepts.png")

# Plot showing effects of degree_days within a site
site.int %>%
  ggplot() +
  geom_density_ridges_gradient(
    aes(x = site_int_dd,
        y = interaction(degree_days, siteID),
        group = interaction(siteID, degree_days),
        fill = siteID),
    scale = 2,
    rel_min_height = 0.01,
    quantile_lines = TRUE, quantiles = 2,
    alpha = 0.5) +
  xlim(c(2.5, 6)) +
  #scale_fill_viridis_d(alpha = 0.5) +
  # geom_vline(aes(xintercept = -1.56),
  #            color = "black",
  #            linetype = "dashed",
  #            size = 1) +
  # geom_vline(aes(xintercept = -1.07),
  #            color = "black",
  #            linetype = "dashed",
  #            size = 1) +
  theme_bw() +
  # facet_wrap(.~siteID,
  #            scales = "free_y")+
  labs(y = "Site",
       x = "Intercept coefficient",
       title = "Effect of degree days on Site-specific \nintercept coefficient distributions",
       subtitle = "") +
  NULL
ggsave("ms/site-dd-intercepts.png")

# intercepts across latitudes?
site.info <- readr::read_csv("data/aquatic-field-sites.csv")
site.int %>%
  group_by(siteID) %>%
  summarize(q50 = quantile(site_int, probs = 0.5)) %>%
  left_join(site.info[, c("Site ID", "Latitude")], by = c("siteID" = "Site ID")) %>%
  ggplot(aes(x = Latitude, y = q50)) +
  geom_point() +
  stat_smooth(method = "lm")

site.int %>%
  group_by(siteID) %>%
  summarize(q50 = quantile(site_int, probs = 0.5)) %>%
  left_join(site.info[, c("Site ID", "Latitude")], by = c("siteID" = "Site ID")) %>%
  select(q50, Latitude) %>%
  cor()

# mean annual temperature?
site.int %>%
  group_by(siteID) %>%
  summarize(q50 = quantile(site_int, probs = 0.5)) %>%
  left_join(site.info[, c("Site ID", "mat.c")], by = c("siteID" = "Site ID")) %>%
  ggplot(aes(x = mat.c, y = q50)) +
  geom_point() +
  stat_smooth(method = "lm")

site.int %>%
  group_by(siteID) %>%
  summarize(q50 = quantile(site_int, probs = 0.5)) %>%
  left_join(site.info[, c("Site ID", "mat.c")], by = c("siteID" = "Site ID")) %>%
  select(q50, mat.c) %>%
  cor()
# fitted model ------------------------------------------------------------


# fitted model
mod.fit <- fitted(mod)
mod.fit <- bind_cols(as_tibble(mod.fit),
          dat[,c("log_mids_center", "degree_days", "siteID",
                 "year", "log_count_corrected")])

# RMSE of fitted model

mod.fit %>%
  mutate(resi = (log_count_corrected - Estimate)**2) %>%
  group_by(siteID, year, degree_days) %>%
  summarize(rmse = sqrt(mean(resi))) %>%
  ungroup() %>%
  summarize(avg.rmse = mean(rmse),
            min.rmse = min(rmse),
            max.rmse = max(rmse))


# plot fitted intervals by site
ggplot(mod.fit,
       aes(ymin = Q2.5, ymax = Q97.5,
           x = log_mids_center,
           group = interaction(siteID, year, degree_days),
           fill = log10(degree_days))) +
  geom_ribbon(alpha = 0.2) +
  scale_fill_viridis_c (option = "plasma") +
  facet_wrap(.~siteID)

# plot all fitted lines
ggplot(mod.fit,
       aes(y = Estimate,
           x = log_mids_center,
           group = interaction(siteID, year, degree_days),
           color = log10(degree_days))) +
  geom_line() +
  scale_color_viridis_c (option = "plasma") +
  theme_bw()


# # make new data for prediction at new sites
# new.site <- dat %>% 
#   select(log_mids_center, degree_days) %>%
#   unique() %>%
#   mutate(siteID = "new", 
#          year = "new")
# 
# pred <- predict(mod, 
#                 type = "response",
#                 newdata = new.site,
#                 allow_new_levels = TRUE,
#                 summary = TRUE,
#                 sample_new_levels = "gaussian")
# pred <- bind_cols(as_tibble(pred), new.site)

# plot fitted intervals with observed points ####

# first, add a new column combining siteID:degree-days for facetting
mod.fit <- mod.fit %>%
  mutate(ldd = round(log10(degree_days),3)) %>%
  unite(facetID, siteID, ldd, remove = FALSE)
  
# plot by facetID
ggplot(mod.fit,
       aes(y = Estimate,
           ymin = Q2.5,
           ymax = Q97.5,
           x = log_mids_center,
           group = interaction(siteID, degree_days),
           fill = log10(degree_days))) +
  theme_bw() +
  geom_ribbon(alpha = 0.6) +
  geom_point(aes(y = log_count_corrected,
                 x = log_mids_center, 
                 fill = log10(degree_days)),
             pch = 21, 
             color = "black",
             size = 1,
             inherit.aes = FALSE,
             alpha = 0.7) +
  facet_wrap(.~facetID,
             ncol = 10) +
  scale_fill_viridis(option = "plasma") +
  scale_color_viridis(option = "plasma") +
  labs(y = "Log10 Normalized Abundance (N)",
       x = "Log10 M (mg, dry weight) centered",
       title = "NEON Macroinvertebrate M-N observations",
       subtitle = "95% fitted interval")

# by site
ggplot(mod.fit,
       aes(y = Estimate,
           ymin = Q2.5,
           ymax = Q97.5,
           x = log_mids_center,
           group = interaction(siteID, degree_days),
           fill = degree_days)) +
  theme_bw() +
  geom_ribbon(alpha = 0.4) +
  geom_point(aes(y = log_count_corrected,
                 x = log_mids_center, 
                 color = degree_days),
             # pch = 21, 
             # color = "black",
             size = 0.9,
             inherit.aes = FALSE,
             alpha = 0.4) +
  facet_wrap(.~siteID,
             ncol = 3) +
  scale_fill_viridis(option = "plasma") +
  scale_color_viridis(option = "plasma") +
  labs(y = "Log10 Normalized Abundance (N)",
       x = "Log10 M (mg, dry weight) centered"#,
       # title = "2017-2018 NEON Macroinvertebrate M-N observations",
       # subtitle = "95% fitted interval"
       )
ggsave("ms/fitted_fig.png",
       width = 6.5,
       height = 7,
       units = "in")


# Predict 2019 data -------------------------------------------------------

# predict new observations for sites that already have been seen

old.site.2019.dat <- dat.2019 %>%
  filter(siteID %in% dat$siteID) %>%
  select(siteID, log_mids_center, degree_days, year, log_count_corrected)

pred.2019 <- predict(mod, 
                     type = "response", 
                     newdata = old.site.2019.dat,
                     allow_new_levels = TRUE)

pred.2019 <- bind_cols(as_tibble(pred.2019),
                       old.site.2019.dat)


# calculate RMSE for pred.2019
pred.2019 %>%
  mutate(resi = (log_count_corrected - Estimate)**2) %>%
  group_by(siteID, year, degree_days) %>%
  summarize(rmse = sqrt(mean(resi))) %>%
  ungroup() %>%
  summarize(avg.rmse = mean(rmse),
            min.rmse = min(rmse),
            max.rmse = max(rmse))

# % observations within prediction interval?
pred.2019 %>%
  mutate(correct = log_count_corrected > Q2.5 &
           log_count_corrected < Q97.5) %>%
  group_by(siteID, degree_days) %>%
  summarize(p.correct = mean(correct)) %>%
  ungroup() %>%
  summarize(avg.percent = mean(p.correct),
            min.percent = min(p.correct),
            max.percent = max(p.correct))

pred.2019 <- pred.2019 %>%
  mutate(ldd = round(log10(degree_days),3)) %>%
  unite(facetID, siteID, ldd, remove = FALSE)

# plot prediction interval and observed points for 2019
ggplot(pred.2019,
       aes(y = Estimate,
           ymin = Q2.5,
           ymax = Q97.5,
           x = log_mids_center,
           group = interaction(siteID, degree_days),
           fill = degree_days)) +
  # geom_line(aes(y = Estimate, x = log_mids_center),
  #           color = "black") +
  theme_bw() +
  geom_ribbon(alpha = 0.3) +
  geom_point(aes(y = log_count_corrected,
                 x = log_mids_center, 
                 color = degree_days),
             size = .9,
             inherit.aes = FALSE,
             alpha = 1) +
  facet_wrap(.~siteID,
             ncol = 3) +
  scale_fill_viridis(option = "plasma") +
  scale_color_viridis(option = "plasma") +
  labs(y = "Log10 Normalized Abundance (N)",
       x = "Log10 M (mg, dry weight) centered",
       title = "NEON Macroinvertebrate M-N 2019 observations",
       subtitle = "95% Prediction interval")
ggsave("ms/predict_2019_old.png",
       width = 6.5,
       height = 7,
       units = "in")

# predict for new sites in 2019
new.site.2019 <- dat.2019 %>%
  filter(!siteID %in% dat$siteID) %>%
  select(siteID, log_mids_center, degree_days, year, log_count_corrected)

new.pred.2019 <- predict(mod, 
                     type = "response", 
                     newdata = new.site.2019,
                     allow_new_levels = TRUE)

new.pred.2019 <- bind_cols(as_tibble(new.pred.2019),
                       new.site.2019)
new.pred.2019 %>%
  mutate(resi = (log_count_corrected - Estimate)**2) %>%
  group_by(siteID, year, degree_days) %>%
  summarize(rmse = sqrt(mean(resi))) %>%
  ungroup() %>%
  summarize(avg.rmse = mean(rmse),
            min.rmse = min(rmse),
            max.rmse = max(rmse))

# % observations within prediction interval?
new.pred.2019 %>%
  mutate(correct = log_count_corrected > Q2.5 &
           log_count_corrected < Q97.5)  %>%
  group_by(siteID, degree_days) %>%
  summarize(p.correct = mean(correct)) %>%
  ungroup() %>%
  summarize(avg.percent = mean(p.correct),
            min.percent = min(p.correct),
            max.percent = max(p.correct))

# plot prediction interval and observed points for 2019
ggplot(new.pred.2019,
       aes(y = Estimate,
           ymin = Q2.5,
           ymax = Q97.5,
           x = log_mids_center,
           group = interaction(siteID, degree_days),
           fill = degree_days)) +
  geom_line(aes(y = Estimate, x = log_mids_center),
            color = "black") +
  theme_bw() +
  geom_ribbon(alpha = 0.4) +
  geom_point(aes(y = log_count_corrected,
                 x = log_mids_center, 
                 fill = degree_days),
             shape = 21, 
             color = "black",
             size = 2,
             inherit.aes = FALSE,
             alpha = 1) +
  facet_wrap(.~siteID) +
  scale_fill_viridis(option = "plasma") +
  scale_color_viridis(option = "plasma") +
  labs(y = "Log10 Normalized Abundance (N)",
       x = "Log10 M (mg, dry weight) centered",
       title = "M-N 2019 new site observations",
       subtitle = "95% Prediction interval")
ggsave("ms/predict_2019_new.png")



# Old figure 2, don't think we're using this anymore. 

# figure 2
# multipanel figure of one site
# top row: 2017 and 2018 fit
# second row: 2019 prediction, 2019 prediction for new site
fig.2a <- mod.fit %>%
    filter(siteID == "ARIK", year == 2017) %>%
    ggplot(aes(y = Estimate,
               ymin = Q2.5,
               ymax = Q97.5,
               x = log_mids_center,
               group = interaction(siteID, degree_days),
               fill = log10(degree_days))) +
    theme_bw() +
    geom_ribbon(alpha = 0.3) +
    geom_point(aes(y = log_count_corrected,
                   x = log_mids_center, 
                   fill = log10(degree_days)),
               pch = 21, 
               color = "black",
               size = 2,
               inherit.aes = FALSE,
               alpha = 0.5) +
    scale_fill_viridis(option = "plasma") +
    scale_color_viridis(option = "plasma") +
    labs(y = "Log10 N",
         x = "Log10 M",
         title = "ARIK 2017 fit")


fig.2b <- mod.fit %>%
    filter(siteID == "ARIK", year == 2018) %>%
    ggplot(aes(y = Estimate,
               ymin = Q2.5,
               ymax = Q97.5,
               x = log_mids_center,
               group = interaction(siteID, degree_days),
               fill = log10(degree_days))) +
    theme_bw() +
    geom_ribbon(alpha = 0.3) +
    geom_point(aes(y = log_count_corrected,
                   x = log_mids_center, 
                   fill = log10(degree_days)),
               pch = 21, 
               color = "black",
               size = 2,
               inherit.aes = FALSE,
               alpha = 0.5) +
    scale_fill_viridis(option = "plasma") +
    scale_color_viridis(option = "plasma") +
    labs(y = "Log10 N",
         x = "Log10 M",
         title = "ARIK 2018 fit")


fig.2c <- pred.2019 %>%
    filter(siteID == "ARIK", year == 2019) %>%
    ggplot(aes(y = Estimate,
               ymin = Q2.5,
               ymax = Q97.5,
               x = log_mids_center,
               group = interaction(siteID, degree_days),
               fill = log10(degree_days))) +
    theme_bw() +
    geom_ribbon(alpha = 0.3) +
    geom_point(aes(y = log_count_corrected,
                   x = log_mids_center, 
                   fill = log10(degree_days)),
               pch = 21, 
               color = "black",
               size = 2,
               inherit.aes = FALSE,
               alpha = 0.5) +
    scale_fill_viridis(option = "plasma") +
    scale_color_viridis(option = "plasma") +
    labs(y = "Log10 N",
         x = "Log10 M",
         title = "ARIK 2019 prediction")


fig.2d <- new.pred.2019 %>%
    filter(siteID == "BLDE", year == 2019) %>%
    ggplot(aes(y = Estimate,
               ymin = Q2.5,
               ymax = Q97.5,
               x = log_mids_center,
               group = interaction(siteID, degree_days),
               fill = log10(degree_days))) +
    theme_bw() +
    geom_ribbon(alpha = 0.3) +
    geom_point(aes(y = log_count_corrected,
                   x = log_mids_center, 
                   fill = log10(degree_days)),
               pch = 21, 
               color = "black",
               size = 2,
               inherit.aes = FALSE,
               alpha = 0.5) +
    scale_fill_viridis(option = "plasma") +
    scale_color_viridis(option = "plasma") +
    labs(y = "Log10 N",
         x = "Log10 M",
         title = "New Site (BLDE) 2019 prediction")


grid.arrange(fig.2a, fig.2b, fig.2c, fig.2d)
g <- arrangeGrob(fig.2a, fig.2b, fig.2c, fig.2d, nrow = 2)
ggsave("ms/fig2.png", g)




# Wesner_Predict slopes-------------------------------------------------

# Get least squares estimates
# compute least squares regressions for each site/year/degree_day combo
dat_formod <- dat_all %>% mutate(group = paste0(siteID,round(degree_days,0), "_",year))

lms <- lmList(log_count_corrected ~ log_mids_center | group, data = dat_formod)

#get coefficients
lms_coefs <- coef(lms) %>% rownames_to_column() %>% clean_names() %>% 
  rename(group = rowname,
         int= intercept,
         slope = log_mids_center) %>% 
  mutate(year = str_sub(group,-4,-1),
         site = str_sub(group, 1, 4),
         degree_days = parse_number(str_sub(group, -50, -6))) %>% 
  glimpse()

#select columns
lms_combine <- lms_coefs %>% select(site, year, degree_days, slope)


# posterior slopes for sites with past data
postsnew4 <- posterior_samples(mod) %>% clean_names() %>% mutate(iter = 1:nrow(.)) %>% 
  select(!contains("cor_site")) %>% 
  pivot_longer(contains("r_site"), names_to = "ranef", values_to = "offset") %>%
  mutate(site = str_sub(ranef, 11,14),
         ranef = str_sub(ranef, 16)) %>% 
  pivot_wider(names_from = ranef, values_from = offset) %>%
  mutate(site = toupper(site)) %>% 
  left_join(lms_combine) %>%
  # mutate(slope_mod = b_log_mids_center + log_mids_center + (b_log_mids_center_log10degree_days + log_mids_center_log10degree_days)*log10(degree_days)) %>%
  # mutate(slope_mod = b_log_mids_center + rnorm(nrow(.), 0, sd_site_id_log_mids_center) + (b_log_mids_center_log10degree_days + rnorm(nrow(.), 0, sd_site_id_log_mids_center_log10degree_days))*log10(degree_days)) %>%
  mutate(slope_mod = b_log_mids_center + log_mids_center) %>%
  mutate(group = paste0(site, round(degree_days,0)),
         model = "siterandom_year",
         prediction_level = "Sites in original model") %>% 
  glimpse()

# posterior slopes for sites with data only in 2019
postsnew4_2019 <- postsnew4 %>% filter(site %in% c("WLOU", "WALK", "REDB")) %>% #just chose 3 existing sites as placeholders
  left_join(lms_combine) %>%
  mutate(site = case_when(site == "WLOU" ~ "SYCA",
                          site == "WALK" ~ "BLUE",
                          site == "REDB" ~ "BLDE"),
         slope_mod = b_log_mids_center + rnorm(nrow(.), 0, sd_site_id_log_mids_center)) %>% 
  # slope_mod = b_log_mids_center + rnorm(nrow(.), 0, sd_site_id_log_mids_center) + (b_log_mids_center_log10degree_days + rnorm(nrow(.), 0, sd_site_id_log_mids_center_log10degree_days))*log10(degree_days)) %>%
  mutate(group = paste0(site, round(degree_days,0)),
         model = "siterandom_year",
         prediction_level = "Sites not in original model") %>% 
  glimpse()

# combine predictions
postsnew4_all <- postsnew4 %>% bind_rows(postsnew4_2019) %>% 
  mutate(order = case_when(prediction_level == "Sites not in original model" ~ -100,
                           TRUE ~ slope_mod))

#make medians for sorting later         
postsnew4_medians <- postsnew4_all %>%
  group_by(site, order) %>% 
  summarize(slope_mod = median(slope_mod))

#add bayes medians to lms coefs
lms_coefs_plot <- left_join(lms_coefs, postsnew4_medians) %>%
  glimpse

# combine posteriors and plot against least squares slopes from all samples. Highlight 2019 samples, since that's what we're trying to predict.
plot_predict_slopes2019 <-  postsnew4_all %>% 
  # drop_na(slope_mod) %>% 
  ggplot(aes(x = slope_mod, y = reorder(site, order))) + 
  geom_halfeyeh(aes(alpha = prediction_level)) +
  geom_point(data = lms_coefs_plot %>% filter(year == 2019), aes(x = slope), size = 2.5, shape = 21) +
  # scale_fill_manual(values = "green") +
  scale_alpha_manual(values = c(1, 0.4)) + 
  # scale_fill_brewer(type = "seq", palette = 1) +
  coord_cartesian(xlim = c(-1.8, -1)) +
  labs(y = "Stream",
       x = "N ~ M slope") +
  theme_bw() +
  theme(legend.title = element_blank())

plot_predict_slopes2019  
ggsave(plot_predict_slopes2019, file = "plots/plot_predict_slopes2019.png", dpi = 500, width = 6, height = 8)



# Wesner_Compare predictions to an overfitted model ------------------------
# simple Bayesian Model - no random effects. Overfitted
my_newpriors <- c(prior(normal(-2, 2), class = "b", coef = "log_mids_center"),
                  prior(normal(0,1), class = "b"),
                  prior(normal(4.25, 2), class = "Intercept"),
                  prior(cauchy(0, 0.1), class = "sigma"))

mod.new2 <- brm(data = dat,
                log_count_corrected ~
                  log_mids_center*log10(degree_days)*siteID, 
                family = gaussian(),
                prior = my_newpriors,
                chains = 1, 
                iter = 2000,
                cores = 4)

postsnew2 <- posterior_samples(mod.new2) %>% as_tibble() %>% mutate(iter = 1:nrow(.)) %>% 
  pivot_longer(contains("siteID"), names_to = "group", values_to = "value") %>% 
  mutate(site = str_sub(group, -4),
         group = paste0("site",str_sub(group, -50, -12))) %>% 
  pivot_wider(names_from = group, values_from = value) %>% 
  left_join(lms_combine) %>% 
  clean_names() %>% 
  # mutate(slope_mod = b_log_mids_center + siteb_log_mids_center + (b_log_mids_center_log10degree_days + siteb_log_mids_center_log10degree_days)*log10(degree_days)) %>%
  mutate(slope_mod = b_log_mids_center + siteb_log_mids_center) %>%
  mutate(group = paste0(site, round(degree_days,0)),
         model = "Overfit",
         prediction_level = "Sites in original model") %>% 
  glimpse()

postsnew2_2019 <- postsnew2 %>% filter(site %in% c("WLOU", "WALK", "REDB")) %>% #just chose 3 existing sites as placeholders
  left_join(lms_combine) %>%
  mutate(site = case_when(site == "WLOU" ~ "SYCA",
                          site == "WALK" ~ "BLUE",
                          site == "REDB" ~ "BLDE"),
         slope_mod = b_log_mids_center) %>% 
  # slope_mod = b_log_mids_center + rnorm(nrow(.), 0, sd_site_id_log_mids_center) + (b_log_mids_center_log10degree_days + rnorm(nrow(.), 0, sd_site_id_log_mids_center_log10degree_days))*log10(degree_days)) %>%
  mutate(group = paste0(site, round(degree_days,0)),
         model = "Overfit",
         prediction_level = "Sites not in original model") %>% 
  glimpse()

postsnew2_all <- postsnew2 %>% bind_rows(postsnew2_2019) %>% 
  mutate(order = case_when(prediction_level == "Site not in original model" ~ -100,
                           TRUE ~ slope_mod))

#plot to show partial pooling
plot_predict_compare2019 <-  postsnew2_all %>% 
  bind_rows(postsnew4_all %>% mutate(model = "Full model")) %>% 
  # drop_na(slope_mod) %>% 
  ggplot(aes(x = slope_mod, y = reorder(site, order))) + 
  geom_halfeyeh(aes(alpha = prediction_level, fill = model)) +
  geom_point(data = lms_coefs_plot %>% filter(year == 2019), aes(x = slope), size = 2.5, shape = 21) +
  scale_fill_brewer(type = "qual")+
  scale_alpha_manual(values = c(1, 0.4)) + 
  # scale_fill_brewer(type = "seq", palette = 1) +
  # coord_cartesian(xlim = c(-1.8, -1)) +
  labs(y = "Stream",
       x = "N ~ M slope") +
  theme_bw() +
  theme(legend.title = element_blank())

plot_predict_compare2019  
ggsave(plot_predict_compare2019, file = "plots/plot_predict_compare2019.png", dpi = 500, width = 6, height = 8)



# estimate prediction error - how far are median predictions from each model from the actual data in 2019
prediction_error <- bind_rows(postsnew2, postsnew4) %>% 
  group_by(model, group) %>% 
  summarize(median_pred = median(slope_mod)) %>% 
  mutate(site = str_sub(group, 1,4),
         degree_days = parse_number(str_sub(group, 5))) %>% 
  left_join(lms_combine) %>% filter(year == 2019) %>% 
  mutate(prediction_error = median_pred - slope)

prediction_error %>% 
  group_by(model) %>% 
  # median_qi(absolute_error = abs(prediction_error), .width = 1) %>% 
  ggplot(aes(x = model, y = prediction_error)) +
  geom_point(position = position_jitter(width = 0.05))
# geom_pointrange(aes(ymin = .lower, ymax = .upper))

