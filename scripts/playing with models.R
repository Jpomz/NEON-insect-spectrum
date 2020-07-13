# modeling filtered binned data

# started out clean, second ~half or so is a mess, started playing around with all sorts of models and got way too carried away. Keeping it here for now in case I need to go back to it, but can probably delete this

# libraries
library(tidyverse)
library(brms)

# load data
# this data set has been filtered to remove the individuals which had damage that affected their measurements
dat <- readRDS("data/binned_no_damage_MN_degree_days_2017_2019.RDS")

# combine date columns into one collect date factor
dat <- dat %>%
  unite(collectDate, year:day, remove = FALSE)

# filter out 2017-2018 data
dat <- dat %>%
  filter(year != 2019)

# Subset data which has degree_days values and more than 200 daily observations per year
dat <- dat %>%
  filter(!is.na(degree_days), 
         yearly.n >= 200)

# make month a factor
dat <- dat %>%
  mutate(month = as.factor(month))

# plot data by siteID
# colorcoded by degree days (as factor)
ggplot(dat,
       aes(y = log_count_corrected, 
           x = log_mids_center,
           color = as.factor(log10(degree_days)))) +
  geom_point() +
  facet_wrap(.~siteID) +
  theme_bw() +
  theme(legend.position = "none") 

# Bayesian model ---------------------------------------------------------

# define my priors
my_priors <- c(prior(normal(-1, 0.5), class = "b"),
               prior(normal(4.5, 1), class = "Intercept"),
               prior(cauchy(0, 0.1), class = "sigma"))
# For models that also include sd priors
my_priors_sd <- c(my_priors, prior(cauchy(0,0.1), class = "sd"))

# random intercept by site and year model
# model log10 N ~ log10 M  * log10(degree_days) + (1|siteID) + (1|year)

model1 <- brm(data = dat,
           log_count_corrected ~
             log_mids_center * log10(degree_days) +
             (1|siteID) + (1|year), 
           family = gaussian(),
           prior = my_priors_sd,
           chains = 4, 
           iter = 2000,
           control = list(adapt_delta = 0.99))

print(model1)
plot(model1)
plot(conditional_effects(model1), points = TRUE)
pp_check(model1)

saveRDS(model1, "data/model1.RDS")

# no year effect
model2 <- update(model1,
                 formula. =  ~ . -(1|year))

# no degree days
model3 <- update(model1,
                 formula. =  ~ . -log10(degree_days) -
                   log_mids_center:log10(degree_days))

saveRDS(model3, "data/no_dd_model.RDS")

model3b <- update(model3,
                  formula. = ~. -log_mids_center + log_mids,
                  newdata = dat)
saveRDS(model3b, "data/no_dd_model_beta.RDS")

# just site
model4 <- update(model1,
                 formula. =  ~ . -log10(degree_days) -
                   log_mids_center:log10(degree_days)
                 -(1|year))


pp_check(model1)
pp_check(model2)
pp_check(model3)
pp_check(model4)

loo.mods <- loo(model1, model2, model3, model4)

# model1 with random slopes
model1.rs <- brm(data = dat,
              log_count_corrected ~
                log_mids_center * log10(degree_days) +
                (1+ log10(degree_days)|siteID) + (1|year), 
              family = gaussian(),
              prior = my_priors_sd,
              chains = 4, 
              iter = 2000,
              control = list(adapt_delta = 0.99))

pp_check(model1.rs, type = "boxplot")
loo.3 <- loo(model1, model1.rs, month.mod3)
# models with month as a factor ####
month.mod <- brm(data = dat,
                 log_count_corrected ~
                   log_mids_center * month +
                   (1|siteID) + (1|year), 
                 family = gaussian(),
                 prior = my_priors_sd,
                 chains = 4, 
                 iter = 2000,
                 control = list(adapt_delta = 0.99))
plot(conditional_effects(month.mod), points = TRUE)

pp_check(month.mod)

month.mod2 <- update(month.mod, 
                     formula. = ~ . +
                       log10(degree_days) +
                       log_mids_center:log10(degree_days),
                     newdata = dat)
plot(conditional_effects(month.mod2))
pp_check(month.mod2)

month.mod3 <-brm(data = dat,
                 log_count_corrected ~
                   log_mids_center * month +
                   (1 + month|siteID) + (1|year), 
                 family = gaussian(),
                 prior = my_priors_sd,
                 chains = 4, 
                 iter = 2000,
                 control = list(adapt_delta = 0.99))
month.mod3
pp_check(month.mod3)
saveRDS(month.mod3, "data/month_fct_rslope_model.RDS")

month.mod.4 <- brm(data = dat,
                   log_count_corrected ~
                     log_mids_center * month *log10(degree_days) +
                     (1 + month*log10(degree_days)|siteID) + (1|year), 
                   family = gaussian(),
                   prior = my_priors_sd,
                   chains = 1, 
                   iter = 1000)

month.mod.5 <- update(month.mod.4,
                      formula. = log_count_corrected ~
                     log_mids_center * month *log10(degree_days) +
                     (1 + month + log10(degree_days)|siteID) + (1|year))

month.mod.6 <- update(month.mod.4,
                      formula. = log_count_corrected ~
                     log_mids_center * month *log10(degree_days) +
                     (1 + log10(degree_days)|siteID) + (1|year))

# loo.month <- loo(month.mod, month.mod2, month.mod3, reloo = TRUE)

month.mod3.fit <- fitted(month.mod3)

month.mod3.fit <- bind_cols(as_tibble(month.mod3.fit),
                            dat[,c("siteID", "log_count_corrected",
                                   "log_mids_center", "month")])

month.mod3.fit %>%
  summarise(rmse =sqrt(
    sum((log_count_corrected - Estimate)**2))
    / n())

# mm3.new <- data.frame(log_mids_center = seq(-3.0, 4.8, length.out = 100),
#                       month = "new",
#                       siteID = "new",
#                       year = "new")

mm3.new <- dat %>% 
  select(log_mids_center, month) %>%
  unique() %>%
  mutate(siteID = "new", 
         year = "new")

month.mod3.pred <- predict(month.mod3, 
                           type = "response",
                           newdata = mm3.new, 
                           re_formula = ~(1 + month | siteID) +
                             (1 | year),
                           allow_new_levels = TRUE,
                           summary = TRUE,
                           sample_new_levels = "gaussian")

month.mod3.pred <- (bind_cols(as_tibble(month.mod3.pred),
                                       mm3.new))

ggplot()+
  geom_ribbon(data = month.mod3.pred,
              aes(x = log_mids_center, 
                  ymin = Q2.5, 
                  ymax = Q97.5),
              alpha = 0.4,
              fill = "grey50") +
  geom_ribbon(data = month.mod3.fit,
              aes(x = log_mids_center, 
                  y = Estimate, 
                  ymin = Q2.5, 
                  ymax = Q97.5,
                  fill = siteID),
              alpha = .4) +
  geom_point(data = month.mod3.fit,
             aes(y = log_count_corrected,
                 x = log_mids_center,
                 color = siteID),
             inherit.aes = FALSE,
             alpha = 0.4) +
  facet_wrap(.~month) +
  theme_bw()

mm3.posts <- posterior_samples(month.mod3)
mm3.posts %>%
  mutate(s2 = b_log_mids_center,
         s3 = b_log_mids_center + `b_log_mids_center:month03`,
         s4 = b_log_mids_center + `b_log_mids_center:month04`,
         s5 = b_log_mids_center + `b_log_mids_center:month05`,
         s6 = b_log_mids_center+ `b_log_mids_center:month06`,
         s7 = b_log_mids_center+ `b_log_mids_center:month07`,
         s8 = b_log_mids_center + `b_log_mids_center:month08`,
         s9 = b_log_mids_center+ `b_log_mids_center:month09`,
         s10 = b_log_mids_center + `b_log_mids_center:month10`) %>%
  select(s2:s10) %>%
  pivot_longer(s2:s10, names_to = "month", values_to = "slope") %>%
  ggplot() +
  geom_density_ridges(
    aes(x = slope,
        y = month,
        #height = ..density..,
        fill = month),
    alpha = 0.6,
    scale = 3) +
  xlim(-2, -0.5)


mod1.fit <- fitted(model1)
mod1.fit <- cbind(obs.y = dat$log_count_corrected,
                  log_mids_center = dat$log_mids_center,
                  degree_days = dat$degree_days,
                  mod1.fit) %>%
  as_tibble()

mod1.fit %>%
  summarise(rmse = sqrt(sum((obs.y - Estimate)**2)) / n())

mod1.new <- data.frame(degree_days = c(min(dat$degree_days),
                                       max(dat$degree_days)),
                       year = "new", 
                       siteID = "new")

# mod1.mids <- data.frame(log_mids_center = 
#                           seq(min(dat$log_mids_center),
#                               max(dat$log_mids_center),length.out = 25))
# mod1.new <- expand_grid(mod1.new, mod1.mids)

mod1.new <- dat %>%
  select(degree_days, log_mids_center) %>%
  unique() %>%
  mutate(year = "new", 
         siteID = "new")

mod1.pred <- predict(model1, 
                     type = "response",
                     newdata = mod1.new, 
                     re_formula = ~(1 | siteID) +
                       (1 | year),
                     allow_new_levels = TRUE,
                     summary = TRUE,
                     sample_new_levels = "gaussian")

mod1.pred <- bind_cols(as_tibble(mod1.pred),
                       mod1.new)

ggplot() +
  geom_ribbon(data = mod1.pred,
              aes(x = log_mids_center, 
                  y = Estimate, 
                  ymin = Q2.5, 
                  ymax = Q97.5,
                  group = degree_days,
                  fill = "predicted"),
              alpha = 0.4) +
  geom_ribbon(data = mod1.fit,
              aes(x = log_mids_center, 
                  y = Estimate, 
                  ymin = Q2.5, 
                  ymax = Q97.5,
                  group = degree_days,
                  fill = "fitted"),
              alpha = 0.4) +
  geom_point(data = mod1.fit,
             aes(x = log_mids_center,
                 y = obs.y)) +
  facet_wrap(.~log10(degree_days)) +
  theme_bw()
