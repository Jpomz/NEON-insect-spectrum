temp <- posterior_samples(mod_widecauchy) %>% clean_names() %>% mutate(iter = 1:nrow(.)) %>% 
  select(!contains("cor_site")) %>% 
  pivot_longer(contains("r_site"), names_to = "ranef", values_to = "offset") %>%
  mutate(site = str_sub(ranef, 11,14),
         ranef = paste0("r_",str_sub(ranef, 16))) %>% 
  pivot_wider(names_from = ranef, values_from = offset) %>%
  mutate(site = toupper(site)) %>% 
  left_join(data_conditioned_on %>% select(siteID, degree_days) %>% 
              rename(site = siteID) %>% distinct(site, degree_days)) %>% 
  mutate(Slope = (b_log_mids_center + r_log_mids_center) + (b_log_mids_center_log10degree_days +
                                                              r_log_mids_center_log10degree_days)*log10(degree_days),
         Intercept = b_intercept + r_intercept) %>%
  mutate(group = paste0(site,round(degree_days,0)),
         prediction_level = "Sites in original model") %>% 
  glimpse()


temp %>% 
  filter(site %in% c("WLOU", "REDB")) %>% 
  group_by(site, degree_days) %>% 
  median_qi(Slope)


wlou_dat <- dat %>% filter(siteID == "REDB") %>% 
  mutate(degree_days = as.integer(degree_days)) %>%
  filter(degree_days == 1522)


WLOU403 <- temp %>% 
  filter(site == "REDB") %>% 
  mutate(degree_days = as.integer(degree_days)) %>% 
  filter(degree_days == 1522) %>% 
  glimpse()

WLOU403_outcome <- WLOU403 %>% 
  mutate(int = b_intercept + r_intercept,
         lmc = (b_log_mids_center + r_log_mids_center)+ 
           (b_log10degree_days + r_log10degree_days)*log10(degree_days) +
           (b_log_mids_center_log10degree_days + r_log_mids_center_log10degree_days)*log10(degree_days)*log_mids_center) %>% 
  group_by(log_mids_center) %>% 
  median_qi(N) %>% 
  mutate(method = "by_hand")



newdata <- dat %>% distinct(siteID, degree_days, year, log_mids_center)

wlou_fitted <- fitted(mod_widecauchy, newdata = newdata) %>% as_tibble() %>% 
  clean_names() %>% mutate(log_mids_center = newdata$log_mids_center,
                           siteID = newdata$siteID,
                           degree_days = newdata$degree_days,
                           year = newdata$year) %>%
  rename(N = estimate, 
         .lower = q2_5,
         .upper = q97_5) %>% 
  mutate(method = "predict()")



fit_checks <- wlou_fitted %>% 
  ggplot(aes(x = log_mids_center)) + 
  geom_ribbon(aes(ymin = .lower, ymax = .upper, group = degree_days), alpha = 0.2) +
  geom_line(aes(y = N, color = degree_days, group = degree_days)) +
  facet_wrap(siteID ~ as.integer(degree_days)) +
  geom_point(data = dat, aes(y = log_count_corrected), shape = 21, size = 0.5) +
  geom_smooth(se = F, data = dat, aes(y = log_count_corrected), method = "lm") +
  theme_void() +
  labs(title = "Posterior predictions (grey) versus unpooled estimates at each site (blue lines)")

ggsave(fit_checks, file = "plots/fit_checks.png")
