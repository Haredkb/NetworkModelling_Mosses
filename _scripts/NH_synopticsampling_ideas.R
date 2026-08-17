##### Moss Proposal pH NH #####################
library(tidyverse)
theme_set(theme_bw())
df <- read.csv("C:/Users/dhare/Dropbox/ResearchProjects/_GeneralData_/NH/GrabSample_PH_Grafton_Carrol_Coos_raw2024.csv")
#filter for stream and station id duplicates

#station <- read.csv("C:/Users/hared/Dropbox/ResearchProjects/_GeneralData_/NH/GrabSample_PH_Grafton_Carrol_Coos_stationfiltered.csv")
station <- read.csv("C:/Users/dhare/Dropbox/ResearchProjects/100_HubbardBrookProjects/1010_Mosses/Proposal_2024/data/final_locations.csv")

station_ids <-  unique(station$Name)

df <- df %>%
  dplyr::filter(STATION.ID %in% station_ids)%>%
  mutate(
    result = as.numeric(QUALIFIER.AND.RESULTS)
  ) %>%
  group_by(STATION.ID)%>%
  reframe(
    x = mean(LATITUDE.DECIMAL.DEGREE, na.rm = TRUE),
    y = mean(LONGITUDE.DECIMAL.DEGREE, na.rm = TRUE),
    huc12 = substr(mean(HUC.12.CODE, na.rm = TRUE), 1, 4),
    pH = median(result, na.rm = TRUE),
    pH_round = round(pH,0)
  )%>%
  left_join(station, by = c("STATION.ID" = "Name"))

a <- ggplot(df)+
  geom_density(aes(x = pH))

b <- ggplot(df)+
  geom_density(aes(x = WETLAND))+
  labs(x = "% Wetland")

c <- ggplot(df)+
  geom_density(aes(x = CONIF))+
  labs(x = "% Conifer")

d <- ggplot(df)+
  geom_density(aes(x = BSLDEM30M))+#mean basin slope
  labs(x = "Mean Basin Slope")

library(cowplot)
plot_grid(a, b, c, d)


obj2 <- station %>%
  left_join(., df)%>%
  drop_na(pH)

write.csv(obj2, "Objective_2_WhiteMountain_LocpH.csv")

a <- ggplot(obj2)+
  geom_density(aes(x = pH))


b <- ggplot(obj2)+
  geom_density(aes(x = WETLAND))+
  labs(x = "% Wetland")

c <- ggplot(obj2)+
  geom_density(aes(x = CONIF))+
  labs(x = "% Conifer")

d <- ggplot(obj2)+
  geom_density(aes(x = BSLDEM30M))+#mean basin slope
  labs(x = "Mean Basin Slope")

plot_grid(a, b, c, d)

sf <- st_as_sf(df, coords = c("y", "x"))

library(maps)

ggplot(sf)+
  geom_sf(aes(color = pH_round))
#0140 Androscoggin
#1060 Saco
#0107 Merrrimack
#0108 Connecticut

ggplot(sf)+
  geom_sf(aes(color = huc12))
