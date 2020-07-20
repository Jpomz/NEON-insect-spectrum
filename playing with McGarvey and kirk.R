# Predicting Other Studies

# "correct" mcgarvey intercepts so they match NEON, see if slopes are the same?

library(tidyverse)
library(brms)

mcgarvey <- read.csv("data/mcgarvey_kirk_SI.csv")
model <- readRDS("data/dd_as_r_slope_6-30-20.RDS")

macros <- mcgarvey %>%
  as_tibble() %>%
  filter(TaxaGroup == 1) %>%
  mutate(Month = as.factor(Month))

# center bin midpoint to match with NEON model
macros <- macros %>%
  group_by(Site, Month) %>%
  mutate(log_mids_center = scale(Log10_BinMidpoint, scale = FALSE)) %>%
    ungroup()

# rename columns so they match with NEON model, and add "year" column
macros <- macros %>%
  rename(siteID = Site) %>%
  mutate(year = "2014",
         degree_days = 955) # median value in NEON data set4

# plot each month within site
ggplot(macros, 
       aes(y = Log10_NormalizedDensity,
           x = log_mids_center,
           color = Month)) +
  geom_point() +
  facet_wrap(.~siteID) +
  theme_bw()

# plot each site within month (match fig 1 in McGarvey and Kirk)
ggplot(macros, 
       aes(y = Log10_NormalizedDensity,
           x = log_mids_center,
           color = siteID)) +
  geom_point() +
  facet_wrap(.~Month) +
  theme_bw()

# predict McGarvey and Kirk data using NEON model
MK.forecast <- predict(model,
                       type = "response",
                       newdata = macros,
                       allow_new_levels = TRUE) %>%
  as_tibble() %>%
  mutate(siteID = macros$siteID,
         Month = macros$Month,
         log_mids_center = macros$log_mids_center,
         year = macros$year,
         obs.y = macros$Log10_NormalizedDensity,
         # "correct" y- values for plots below
         y_plus_3 = obs.y + 3,
         y_plus_4 = obs.y + 4,
         y_plus_5 = obs.y + 5)

MK.rmse <- MK.forecast %>%
  mutate(resi.orginal = (obs.y - Estimate)**2,
         resi.3 = (y_plus_3 - Estimate)**2,
         resi.4 = (y_plus_4 - Estimate)**2,
         resi.5 = (y_plus_5 - Estimate)**2,) %>%
  group_by(siteID, Month) %>%
  summarize(rmse.original = sqrt(sum(resi.orginal) / n()),
            rmse.3 = sqrt(sum(resi.3) / n()),
            rmse.4 = sqrt(sum(resi.4) / n()),
            rmse.5 = sqrt(sum(resi.5) / n()),
            n = n())%>%
  ungroup() %>%
  summarize(avg.rmse.o = mean(rmse.original),
            min.rmse.o = min(rmse.original),
            max.rmse.o = max(rmse.original),
            avg.rmse.3 = mean(rmse.3),
            min.rmse.3 = min(rmse.3),
            max.rmse.3 = max(rmse.3),
            avg.rmse.4 = mean(rmse.4),
            min.rmse.4 = min(rmse.4),
            max.rmse.4 = max(rmse.4),
            avg.rmse.5 = mean(rmse.5),
            min.rmse.5 = min(rmse.5),
            max.rmse.5 = max(rmse.5))

points.in.pred.int <- MK.forecast %>%
  mutate(correct.original = obs.y > Q2.5 & obs.y < Q97.5,
         correct.3 = y_plus_3 > Q2.5 & y_plus_3 < Q97.5,
         correct.4 = y_plus_4 > Q2.5 & y_plus_4 < Q97.5,
         correct.5 = y_plus_5 > Q2.5 & y_plus_5 < Q97.5) %>%
  group_by(siteID, Month) %>%
  summarize(n.correct.o = sum(correct.original) / n(),
            n.correct.3 = sum(correct.3) / n(),
            n.correct.4 = sum(correct.4) / n(),
            n.correct.5 = sum(correct.5) / n())


# plot prediction interval with observed points from McGarvey and Kirk
ggplot(MK.forecast,
       aes(y = Estimate,
           ymin = Q2.5,
           ymax = Q97.5,
           x = log_mids_center,
           group = siteID)) +
  geom_line() +
  geom_ribbon(alpha = 0.2) +
  labs(y = "log10N",
       title = "Prediction interval, McGarvey and Kirk observations") +
  theme_bw() +
  geom_point(aes(x = log_mids_center,
                 y = obs.y,
                 color = Month),
             inherit.aes = FALSE) +
  # stat_smooth(aes(x = log_mids_center,
  #                 y = obs.y,
  #                 color = Month),
  #             method = "lm",
  #             se = FALSE,
  #             inherit.aes = FALSE) +
  facet_wrap(.~siteID) +
  theme(legend.position = "none") +
  NULL

# plot "corrected" y-values
# here we are looking for general aggreement in the relationship, after "raising" the MK line up to the NEON community abundance level
ggplot(MK.forecast,
       aes(y = Estimate,
           ymin = Q2.5,
           ymax = Q97.5,
           x = log_mids_center,
           group = siteID)) +
  geom_line() +
  geom_ribbon(alpha = 0.2) +
  labs(y = "Normalized density",
       title = "Prediction interval, McGarvey and Kirk observations",
       subtitle = "Orignal y-values + 3") +
  theme_bw() +
  geom_point(aes(x = log_mids_center,
                 y = y_plus_3,
                 color = Month),
             inherit.aes = FALSE) +
  facet_wrap(.~siteID) +
  theme(legend.position = "none") +
  NULL

ggplot(MK.forecast,
       aes(y = Estimate,
           ymin = Q2.5,
           ymax = Q97.5,
           x = log_mids_center,
           group = siteID)) +
  geom_line() +
  geom_ribbon(alpha = 0.2) +
  labs(y = "Normalized density",
       title = "Prediction interval, McGarvey and Kirk observations",
       subtitle = "Orignal y-values + 4") +
  theme_bw() +
  geom_point(aes(x = log_mids_center,
                 y = y_plus_4,
                 color = Month),
             inherit.aes = FALSE) +
  facet_wrap(.~siteID) +
  theme(legend.position = "none") +
  NULL

ggplot(MK.forecast,
       aes(y = Estimate,
           ymin = Q2.5,
           ymax = Q97.5,
           x = log_mids_center,
           group = siteID)) +
  geom_line() +
  geom_ribbon(alpha = 0.2) +
  labs(y = "Normalized density",
       title = "Prediction interval, McGarvey and Kirk observations",
       subtitle = "Orignal y-values + 5") +
  theme_bw() +
  geom_point(aes(x = log_mids_center,
                 y = y_plus_5,
                 color = Month),
             inherit.aes = FALSE) +
  facet_wrap(.~siteID) +
  theme(legend.position = "none") +
  NULL

