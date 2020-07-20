# predict M-N relationship at 2017-2019 sites without temperature data

library(tidyverse)
library(brms)
library(viridis)

# load full binned data set
dat <- readRDS("data/binned_no_damage_MN_degree_days_2017_2019.RDS")

# load model
#model <- readRDS("data/dd_as_r_slope_6-24-20.RDS") # not actually random slopes
model <- readRDS("data/dd_as_r_slope_6-30-20.RDS") # r slopes after talking with JSW 6/30
# filter data
dat.t.na <- dat %>% 
  filter(is.na(degree_days))

dat.few.t <- dat %>%
  filter(yearly.n < 200)

newdata <- bind_rows(dat.t.na, dat.few.t)

# remove lake and large river sites
lakes <- c("BARC", "CRAM", "LIRO", "PRLA", "PRPO",
           "SUGG", "TOOK", "BLWA", "FLNT", "TOMB")

newdata <- newdata %>%
  filter(!siteID %in% lakes)

# add median degree_days value from modeled data to this data set
newdata$degree_days <- 955

no.t.forecast <- predict(model,
                       type = "response",
                       newdata = newdata,
                       allow_new_levels = TRUE,
                       sample_new_levels = "gaussian") %>%
  as_tibble() %>%
  mutate(siteID = newdata$siteID,
         log_mids_center = newdata$log_mids_center,
         year = newdata$year,
         month = newdata$month,
         obs.y = newdata$log_count_corrected)

no.t.forecast %>%
  mutate(resi = (obs.y - Estimate)**2) %>%
  group_by(siteID, month, year) %>%
  summarize(rmse = sqrt(sum(resi) / n())) %>%
  ungroup() %>%
  summarize(avg.rmse = mean(rmse),
            min.rmse = min(rmse),
            max.rmse = max(rmse))

# % observations within prediction interval?
no.t.forecast %>%
  mutate(correct = obs.y > Q2.5 &
           obs.y < Q97.5)  %>%
  group_by(siteID, month, year) %>%
  summarize(p.correct = mean(correct)) %>%
  ungroup() %>%
  summarize(avg.percent = mean(p.correct),
            min.percent = min(p.correct),
            max.percent = max(p.correct))

# which sites poorly predicted?
no.t.forecast %>%
  mutate(correct = obs.y > Q2.5 &
           obs.y < Q97.5)  %>%
  group_by(siteID, month, year) %>%
  summarize(p.correct = mean(correct)) %>%
  arrange(p.correct)
# TECR in month = 09

# plot predictions
ggplot(no.t.forecast,
       aes(y = Estimate,
           ymin = Q2.5,
           ymax = Q97.5,
           x = log_mids_center,
           group = siteID)) +
  geom_line() +
  geom_ribbon(alpha = 0.4,
              fill = "grey50") +
  labs(y = "log10N",
       x = "Log10 M (mg, dry weight) centered",
       title = "Prediction interval, NEON sites (2017-2019) without Temperature data") +
  theme_bw() +
  geom_point(aes(x = log_mids_center,
                 y = obs.y,
                 color = month,
                 shape = year),
             size = 1.7,
             alpha = 0.6,
             inherit.aes = FALSE) +
  facet_wrap(siteID~.,
             ncol = 3) +
  scale_color_viridis_d(option = "plasma") +
  NULL
ggsave("ms/fig3.png",
       width = 6.5,
       height = 7,
       units = "in")

# poorest prediction was at TECR, plot here for SI
no.t.forecast %>% 
  filter(siteID == "TECR", month == "09") %>%
  ggplot(aes(y = Estimate,
             ymin = Q2.5,
             ymax = Q97.5,
             x = log_mids_center,
             group = siteID)) +
  geom_line() +
  geom_ribbon(alpha = 0.4,
              fill = "grey50") +
  labs(y = "log10N",
       x = "Log10 M (mg, dry weight) centered",
       title = "Prediction interval, TECR, 2019-09") +
  theme_bw() +
  geom_point(aes(x = log_mids_center,
                 y = obs.y),
             color = "red",
             size = 2,
             alpha = 0.6,
             inherit.aes = FALSE) +
  NULL
ggsave("ms/TECR_prediction_SI.png")

# second worst
no.t.forecast %>% 
  filter(siteID == "LECO", month == "10", year == "2018") %>%
  ggplot(aes(y = Estimate,
             ymin = Q2.5,
             ymax = Q97.5,
             x = log_mids_center,
             group = siteID)) +
  geom_line() +
  geom_ribbon(alpha = 0.4,
              fill = "grey50") +
  labs(y = "log10N",
       x = "Log10 M (mg, dry weight) centered",
       title = "Prediction interval, LECO, 2018-10") +
  theme_bw() +
  geom_point(aes(x = log_mids_center,
                 y = obs.y),
             color = "red",
             size = 2,
             alpha = 0.6,
             inherit.aes = FALSE) +
  NULL
ggsave("ms/LECO_prediction_SI.png")



# third worst was OKSR in 08
no.t.forecast %>% 
  filter(siteID == "OKSR", month == "08", year == "2017") %>%
  ggplot(aes(y = Estimate,
             ymin = Q2.5,
             ymax = Q97.5,
             x = log_mids_center,
             group = siteID)) +
  geom_line() +
  geom_ribbon(alpha = 0.4,
              fill = "grey50") +
  labs(y = "log10N",
       x = "Log10 M (mg, dry weight) centered",
       title = "Prediction interval, OKSR, 2017-08") +
  theme_bw() +
  geom_point(aes(x = log_mids_center,
                 y = obs.y),
             color = "red",
             size = 2,
             alpha = 0.6,
             inherit.aes = FALSE) +
  NULL
ggsave("ms/OKSR_prediction_SI.png")

# relationship holds for lakes and large rivers?
non.stream.data <- bind_rows(dat.t.na, dat.few.t)

# keep lake and large river sites
lakes <- c("BARC", "CRAM", "LIRO", "PRLA", "PRPO",
           "SUGG", "TOOK", "BLWA", "FLNT", "TOMB")

non.stream.data <- non.stream.data %>%
  filter(siteID %in% lakes)

# TOMB in 2019-03 only has a single data point, removing
non.stream.data <- non.stream.data %>%
  filter(!((siteID == "TOMB" & month == "03" & year =="2019")))

# add median degree_days value from modeled data to this data set
non.stream.data$degree_days <- 955

non.stream.forecast <- predict(model,
                               type = "response",
                               newdata = non.stream.data,
                               allow_new_levels = TRUE,
                               sample_new_levels = "gaussian") %>%
  as_tibble() %>%
  mutate(siteID = non.stream.data$siteID,
         log_mids_center = non.stream.data$log_mids_center,
         year = non.stream.data$year,
         month = non.stream.data$month,
         obs.y = non.stream.data$log_count_corrected)

non.stream.forecast %>%
  mutate(resi = (obs.y - Estimate)**2) %>%
  group_by(siteID, month, year) %>%
  summarize(rmse = sqrt(sum(resi) / n())) %>%
  ungroup() %>%
  summarize(avg.rmse = mean(rmse),
            min.rmse = min(rmse),
            max.rmse = max(rmse))

# % observations within prediction interval?
non.stream.forecast %>%
  mutate(correct = obs.y > Q2.5 &
           obs.y < Q97.5) %>%
  summarize(p.correct = mean(correct))

ggplot(non.stream.forecast,
       aes(y = Estimate,
           ymin = Q2.5,
           ymax = Q97.5,
           x = log_mids_center,
           group = siteID)) +
  geom_line() +
  geom_ribbon(alpha = 0.2) +
  labs(y = "log10N",
       title = "Prediction interval, NEON non-stream sites (2017-2019)",
       subtitle = "No temperature data") +
  theme_bw() +
  geom_point(aes(x = log_mids_center,
                 y = obs.y,
                 color = month,
                 shape = year),
             size = 2,
             alpha = 0.5,
             inherit.aes = FALSE) +
  facet_wrap(siteID~.) +
  scale_color_viridis_d(option = "plasma") +
  NULL
