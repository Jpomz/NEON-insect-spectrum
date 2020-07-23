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

temp <- posterior_samples(mod) %>% clean_names() %>% mutate(iter = 1:nrow(.)) %>% 
  select(!contains("cor_site")) %>% 
  pivot_longer(contains("r_site"), names_to = "ranef", values_to = "offset") %>%
  mutate(site = str_sub(ranef, 11,14),
         ranef = paste0("r_",str_sub(ranef, 16))) %>% 
  pivot_wider(names_from = ranef, values_from = offset) %>%
  mutate(site = toupper(site)) %>% 
  # pivot_longer(cols = c("r_year_2017_intercept", "r_year_2018_intercept"),names_to = "mod_year", values_to = "r_intercept_year") %>% 
  left_join(data_conditioned_on %>% select(siteID, degree_days, year) %>% 
              rename(site = siteID) %>% distinct(site, degree_days, year)) %>%
  # mutate(year = parse_number(mod_year)) %>% 
  mutate(Slope = (b_log_mids_center + r_log_mids_center) + (b_log_mids_center_log10degree_days +
                                                              r_log10degree_days_log_mids_center)*log10(degree_days),
         Intercept = b_intercept + r_intercept + (b_log10degree_days + r_log10degree_days)*log10(degree_days)) %>%
  mutate(Intercept = case_when(year == "2017" ~ Intercept + r_year_2017_intercept, 
                           TRUE ~ Intercept + r_year_2018_intercept), 
         group = paste0(site,round(degree_days,0)),
         prediction_level = "Sites in original model") %>% 
  glimpse()


newdata <- dat %>% distinct(siteID, degree_days, year, log_mids_center)

mod_predict <- fitted(mod_widecauchy, newdata = newdata) %>% as_tibble() %>% 
  clean_names() %>% mutate(log_mids_center = newdata$log_mids_center,
                           siteID = newdata$siteID,
                           degree_days = newdata$degree_days,
                           year = newdata$year) %>%
  rename(N = estimate, 
         .lower = q2_5,
         .upper = q97_5) %>% 
  mutate(method = "predict()")



fit_checks <- mod_predict %>% 
  ggplot(aes(x = log_mids_center)) + 
  geom_ribbon(aes(ymin = .lower, ymax = .upper, group = degree_days), alpha = 0.2) +
  geom_line(aes(y = N, color = degree_days, group = degree_days)) +
  facet_wrap(siteID ~ as.integer(degree_days)) +
  geom_point(data = dat, aes(y = log_count_corrected), shape = 21, size = 0.5) +
  geom_smooth(se = F, data = dat, aes(y = log_count_corrected), method = "lm") +
  theme_void() +
  labs(title = "Posterior predictions (grey) versus unpooled estimates at each site (blue lines)")

fit_checks

ggsave(fit_checks, file = "plots/fit_checks.png")




# extract predictions by hand
dat_to_add <- dat %>% select(siteID, degree_days) %>% distinct(siteID, degree_days) %>% 
  rename(site = siteID)

by_hand_predict <- temp %>% 
  filter(iter <=1000) %>% 
  expand_grid(log_mids_center = seq(-2, 3, by = 0.15))  %>%
  mutate(N = Intercept + Slope*log_mids_center) %>% 
  group_by(site, degree_days, log_mids_center) %>% 
  mean_qi(N) %>% 
  rename(siteID = site) %>% 
  mutate(method = "by_hand")


mod_predict %>% 
  bind_rows(by_hand_predict) %>%
  filter(siteID == "ARIK") %>% 
  ggplot(aes(x = log_mids_center, y = N, fill = method, color = method)) + 
  geom_line() +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.2) +
  facet_wrap(~degree_days)




# Compare pooling ---------------------------------------------------------

dat_arik <- dat %>% filter(siteID == "ARIK")

# Model with no pooling (separate regressions at each site)
dat_lists <- split(dat, f = dat$siteID)

# mod_nopool <- brm_multiple(log_count_corrected ~ log_mids_center*log10(degree_days),
#                   family = gaussian(),
#                   data = dat_lists,
#                   prior = c(prior(normal(-1, 1),
#                                   class = "b", coef = "log_mids_center"),
#                             prior(normal(0,1), class = "b"),
#                             # prior(normal(4.5, 1.5), class = "Intercept"),
#                             prior(cauchy(0,1), class = "sigma")),
#                   chains = 2,
#                   iter = 2000)
# 
# saveRDS(mod_nopool, file = "data/mod_nopool.rds")
mod_nopool <- readRDS(file = "data/mod_nopool.rds")

# 
# mod_arik <- brm(log_count_corrected ~ log_mids_center*log10(degree_days)*siteID,
#     family = gaussian(),
#     data = dat,
#     prior = c(prior(normal(-1, 1),
#                     class = "b", coef = "log_mids_center"),
#               prior(normal(0,1), class = "b"),
#               # prior(normal(4.5, 1.5), class = "Intercept"),
#               prior(cauchy(0,1), class = "sigma")),
#     chains = 1,
#     iter = 1000)

sites_rep <- dat %>% select(siteID) %>% distinct(siteID) %>% slice(rep(1:n(), each=2000))


post_het <- posterior_samples(mod_nopool, summary = F) %>% as_tibble() %>% clean_names() %>% mutate(iter = 1:nrow(.),
                                                                                                  siteID = sites_rep$siteID) %>% 
  select(-sigma, -lp) %>%
  left_join(data_conditioned_on %>% select(siteID, degree_days) %>% distinct(siteID, degree_days)) %>% 
  mutate(value = b_log_mids_center + b_log_mids_center_log10degree_days*log10(degree_days)) %>% 
  mutate(prediction_level = "Sites in original model",
         Parameter = "Slope",
         method = "No pooling") %>% 
  glimpse()


# sitedd.all %>% mutate(method = "Hierarchical") %>% 
#   bind_rows(post_het) %>% 
#   filter(Parameter == "Slope") %>% 
#   ggplot() +
#   geom_density_ridges_gradient(aes(x = value, y = reorder(siteID, -median), 
#                                    fill = degree_days, group = interaction(degree_days, siteID)),
#                                scale = 2,
#                                rel_min_height = 0.01,
#                                quantile_lines = TRUE, quantiles = 2,
#                                alpha = 0.5) +
#   facet_wrap(~method, scales = "free_x") +
#   xlim(c(-1.8, -0.9)) +
#   scale_fill_viridis(alpha = 0.5) +
#   # geom_vline(aes(xintercept = -1.56),
#   #            color = "black",
#   #            linetype = "dashed",
#   #            size = 1) +
#   # geom_vline(aes(xintercept = -1.07),
#   #            color = "black",
#   #            linetype = "dashed",
#   #            size = 1) +
#   theme_bw() +
#   # guides(fill = F) +
#   labs(y = "Site",
#        x = "Slope",
#        subtitle = "") +
#   NULL



pool_nopool <- sitedd.all %>% mutate(method = "Hierarchical") %>% 
  bind_rows(post_het) %>% 
  filter(Parameter == "Slope") %>% 
  group_by(method, degree_days) %>% 
  median_qi(value) %>% 
  ggplot(aes(x = degree_days, y = value, ymin = .lower, ymax = .upper)) +
  geom_pointrange(aes(color = method), position = position_dodge(width = 0.2)) +
  facet_wrap(~method) +
  labs(y = "Slope") +
  guides(color = F) +
  theme_bw()

ggsave(pool_nopool, file = "plots/pool_nopool.png", dpi = 600, width = 7, height = 5)
