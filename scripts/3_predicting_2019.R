# predict 2019 data

# load libraries
library(tidyverse)
library(brms)
library(ggridges)
library(viridis)

# Load data
full.dat <- readRDS("data/binned_no_damage_MN_degree_days_2017_2019.RDS")
# load model1
model1 <- readRDS("data/model1.RDS")

# combine date columns into one collect date factor
full.dat <- full.dat %>%
  unite(collectDate, year:day, remove = FALSE)

# Subset data which has degree_days values and more than 200 daily observations per year
full.dat <- full.dat %>%
  filter(!is.na(degree_days), 
         yearly.n >= 200)

# split data into modeled(2017-2018) and prediction (2019)
model.dat <- full.dat %>%
  filter(year != 2019)

pred.dat <- full.dat %>%
  filter(year == 2019)

# site-specific predictions ####

site.data <- pred.dat %>%
  filter(siteID %in% model.dat$siteID)

site.forecast <- predict(model1,
                         newdata = site.data,
                         allow_new_levels = TRUE) %>%
  as_tibble() %>%
     mutate(siteID = site.data$siteID,
            log_mids_center = site.data$log_mids_center,
            year = site.data$year,
            degree_days = site.data$degree_days,
            obs.y = site.data$log_count_corrected)
 

rmse.model1 <- site.forecast %>%
  filter(!is.na(obs.y)) %>%
  mutate(resi = (obs.y - Estimate)**2) %>%
  summarize(rmse = sqrt(sum(resi) / n()),
            n = n()) 

# points within prediction interval
points.in.pred.int <- site.forecast %>%
  mutate(correct = obs.y > Q2.5 & obs.y < Q97.5) %>%
  group_by(siteID, degree_days) %>%
  summarize(n.correct = sum(correct) / n())
  
sort(points.in.pred.int$n.correct)
# ~ 50% of predicted site:degree points are within prediction interval
# all site:degree > 80% correct predictions
mean(points.in.pred.int$n.correct)
# mean of 93% of points are within prediction interval
plot(density(rmse.model1$rmse))

# plot site specific prediction intervals
ggplot(site.forecast,
       aes(y = Estimate,
           ymin = Q2.5,
           ymax = Q97.5,
           x = log_mids_center,
           group = degree_days)) +
  geom_line() +
  geom_ribbon(alpha = 0.2) +
  labs(y = "log10N",
       title = "Site-specific prediction intervals") +
  theme_bw() +
  facet_wrap(.~siteID)

# plot prediction interval and observed 2019 data
ggplot(site.forecast,
        aes(y = Estimate,
            ymin = Q2.5,
            ymax = Q97.5,
            x = log_mids_center,
            group = degree_days,
            color = as.factor(degree_days),
            fill = as.factor(degree_days))) +
  geom_line() +
  geom_ribbon(alpha = 0.2) +
  labs(y = "log10N",
       title = "Site-specific prediction interval and observed points for 2019") +
  theme_bw() +
  geom_point(aes(x = log_mids_center,
                 y = obs.y,
                 color = as.factor(degree_days)),
             inherit.aes = FALSE) +
  facet_wrap(.~siteID) +
  theme(legend.position = "none") +
  NULL

# new site predictions ####
new.site.data <- pred.dat %>%
  filter(!siteID %in% model.dat$siteID)

new.site.forecast <- predict(model1,
                         newdata = new.site.data,
                         allow_new_levels = TRUE) %>%
  as_tibble() %>%
  mutate(siteID = new.site.data$siteID,
         log_mids_center = new.site.data$log_mids_center,
         year = new.site.data$year,
         degree_days = new.site.data$degree_days,
         obs.y = new.site.data$log_count_corrected)


rmse.new.site <- new.site.forecast %>%
  filter(!is.na(obs.y)) %>%
  mutate(resi = (obs.y - Estimate)**2) %>%
  summarize(rmse = sqrt(sum(resi) / n()),
            n = n()) 

# points within prediction interval
points.in.pred.int <- new.site.forecast %>%
  mutate(correct = obs.y > Q2.5 & obs.y < Q97.5) %>%
  group_by(siteID, degree_days) %>%
  summarize(n.correct = sum(correct) / n())

sort(points.in.pred.int$n.correct)
# all site:degree > 90% correct predictions
mean(points.in.pred.int$n.correct)
# mean of ~95% of points are within prediction interval


# plot site specific prediction intervals
ggplot(new.site.forecast,
       aes(y = Estimate,
           ymin = Q2.5,
           ymax = Q97.5,
           x = log_mids_center,
           group = degree_days)) +
  geom_line() +
  geom_ribbon(alpha = 0.2) +
  labs(y = "log10N",
       title = "Site-specific prediction intervals") +
  theme_bw() +
  facet_wrap(.~siteID)

# plot prediction interval and observed 2019 data
ggplot(new.site.forecast,
       aes(y = Estimate,
           ymin = Q2.5,
           ymax = Q97.5,
           x = log_mids_center,
           group = degree_days,
           color = as.factor(degree_days),
           fill = as.factor(degree_days))) +
  geom_line() +
  geom_ribbon(alpha = 0.2) +
  labs(y = "log10N",
       title = "New Site prediction interval and observed points for 2019") +
  theme_bw() +
  geom_point(aes(x = log_mids_center,
                 y = obs.y,
                 color = as.factor(degree_days)),
             inherit.aes = FALSE) +
  facet_wrap(.~siteID) +
  theme(legend.position = "none") +
  NULL


# slopes across sites

mod1.posts <- as_tibble(posterior_samples(model1))
mod1.slopes <- mod1.posts %>% select(b_log_mids_center, `b_log_mids_center:log10degree_days`)
model.site.dd <- model.dat %>%
  select(siteID, degree_days) %>%
  unique()
mod1.slopes <- expand_grid(mod1.slopes, model.site.dd)

mod1.slopes %>%
  mutate(slope = b_log_mids_center + `b_log_mids_center:log10degree_days`*log10(degree_days)) %>%
  ggplot() +
  geom_density_ridges(
    aes(x = slope,
        y = as.factor(round(log10(degree_days), 2)),
        height = ..density..,
        fill = as.factor(round(log10(degree_days), 2))),
    alpha = 0.6,
    scale = 3) +
  scale_x_continuous(limits= c(-1.5, -1.2)) +
  scale_fill_viridis(discrete = TRUE) 

# no year / degree days data?
test.data <- tibble(siteID = c("ARIK", "NEW1"),
                    degree_days =
                      c(1895,
                        mean(model.dat$degree_days)), #estimateBasedOnMean?
                    year = "2019")

test.data <- expand_grid(test.data, log_mids_center = unique(site.data$log_mids_center))

predict(model1, newdata = test.data, allow_new_levels = TRUE)
