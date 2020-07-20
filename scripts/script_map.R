library(tidyverse)

# get coordinates for polygons - function in ggplot
world <- map_data("world")
states <- map_data("state")

#get NEON coords (uses old dataset I made with lat/longs)
dat_all <- readRDS(file = "data/dat_all.rds")

neon_latlong <- read_csv(file = "data/neon_latlong.csv") %>% 
  select(-X1) %>% 
  right_join(dat_all %>% distinct(siteID) %>% rename(site = siteID))

#make the map
map <- ggplot() + 
  geom_polygon(data = world, aes(x = long, y = lat, group = group), color = "white") +
  geom_polygon(data = states, aes(x = long, y = lat, group =group), color = "white") +
  # geom_polygon(data = alaska, aes(x = long, y = lat, group =group), color = "white") +
  coord_quickmap(ylim = c(18,70),
      xlim = c(-160,-50))+
  geom_point(data = neon_latlong, aes(x = long, y = lat),size = 3, shape = 21,
             fill = "yellow") +
  theme_map()

write.csv(neon_latlong, file = "data/neon_latlong.csv")
ggsave(map, file = "plots/map.png", width = 10, height = 8, dpi = 600)
