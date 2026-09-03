rm(list = ls())

# PACKAGES AND SOURCE  ----
suppressPackageStartupMessages(library(tidyverse))
source("src/01a-utils.R")

# LOAD DATA AND FILTER ----
## Load data ----
data = readRDS("data/bike-sharing/milan_bikemi_data_for_regression.rds") |> 
  mutate(Weekend = as.numeric(weekdays(commit_at) %in% c("Saturday", "Sunday")),
         Rain = as.numeric(Precipitation > 0))
## Filter 2025 and no bikes ----
nobikes2025 = data |> filter(year(commit_at) == 2025 & unavail_type == "no_bikes")

# Exclude stations close to Pero: 334, 500 and 501. Actually 500 does not have any observation in 2025.
D = nobikes2025 |> filter(ID != 334, ID != 501) |> 
  select(-ID, -commit_at, -Precipitation, -Temperature, -Rain)


## Spatial component ----
X_spat = D |> select(station, longitude, latitude)
X_spat_unique = X_spat |> distinct()
coord_utm = sf::st_as_sf(X_spat_unique, 
                         coords = c("longitude", "latitude"), crs = 4326) |> 
  sf::st_transform(crs = 32632) |> 
  sf::st_coordinates()


# # Option 1: Normalize with observed observed min and max x and y coordinates
min_x = min(coord_utm[ , 1]); max_x = max(coord_utm[ , 1])
min_y = min(coord_utm[ , 2]); max_y = max(coord_utm[ , 2])
coord_scaled = coord_utm
coord_scaled[ , 1] = (coord_utm[ , 1] - min_x) / (max_x - min_x)
coord_scaled[ , 2] = (coord_utm[ , 2] - min_y) / (max_y - min_y)
milano_coords_scaled = milano_coords
milano_coords_scaled[ , 1] = (milano_coords[ , 1] - min_x) / (max_x - min_x)
milano_coords_scaled[ , 2] = (milano_coords[ , 2] - min_y) / (max_y - min_y)

n_bases = 5
out_spatial_DR = construct_DR_Spatial(coord_scaled[ , 1], coord_scaled[ , 2], 
                                      k1 = n_bases, k2 = n_bases)

x1=coord_scaled[,1]
x2=coord_scaled[,2]

ngrid = 100
t1=seq(0, 1, length=ngrid)
t2=seq(0, 1,length=ngrid)

TT=expand.grid(t1,t2)

t1n=TT[,1]
t2n=TT[,2]

BT1=smoothCon(s(t1n,bs="ps",k=n_bases),data=data.frame(t1n))[[1]]$X
BT2=smoothCon(s(t2n,bs="ps",k=n_bases),data=data.frame(t2n))[[1]]$X

BT=tensor.prod.model.matrix(list(BT1,BT2))
XT=BT%*%out_spatial_DR$Trafo

plt_DR = plot_spatial_splines(coord = coord_scaled, grid = TT, bases = XT, type = "DR",
                              idx_to_plot = 1:(n_bases^2-1), normalize = F, 
                              boundaries = milano_coords_scaled,
                              title = "2D DR Splines")

plt_Bspl = plot_spatial_splines(coord = coord_scaled, grid = TT, bases = BT, type = "B-Splines",
                                idx_to_plot = 1:n_bases^2, normalize = F, 
                                boundaries = milano_coords_scaled,
                                title = "2D B-Splines")

ggsave(filename = "spatial_basis_DR.pdf",
       path = "application/bike-sharing/pics/data",
       plot = plt_DR, width = 10, height = 7)

ggsave(filename = "spatial_basis_BSplines.pdf",
       path = "application/bike-sharing/pics/data",
       plot = plt_Bspl, width = 10, height = 7)
