rm(list = ls())
library(mapview)
library(sf)
library(arrow)
library(terra)
library(exactextractr)
suppressPackageStartupMessages(library(tidyverse))

PLOT = TRUE

voronoi_utm = readRDS("data/bike-sharing/milan_voronoi.rds")
poi = readRDS("data/bike-sharing/milan_poi.rds")
stations_utm = readRDS("data/bike-sharing/milan_stations_utm.rds")
bike_history = readRDS("data/bike-sharing/milan_bikemi_history.rds")
# weather = readRDS("data/bike-sharing/milan_weather_ASPRA_15mins.rds")

# Compute population density ----
pop <- rast("data/bike-sharing/ita_pop_2025_CN_100m_R2024B_v1.tif") # download from https://hub.worldpop.org/geodata/summary?id=55090
pop_utm <- project(pop, crs(voronoi_utm))

voronoi_utm$population <- exact_extract(
  pop_utm,
  voronoi_utm,
  "sum"
)
voronoi_utm$population_density = voronoi_utm$population / voronoi_utm$area_km2

# Check that the population density in the Voronoi cells resembles the population density in the raster
if (PLOT){
  area_sf <- st_transform(voronoi_utm, crs(pop_utm))
  area_vect <- vect(area_sf)
  pop_crop <- crop(pop_utm, area_vect)
  pop_mask <- mask(pop_crop, area_vect)
  plot(pop_mask)
  mapview(pop_mask)
  mapview(voronoi_utm, zcol = "population_density")
}

# Plot Voronoi and stations ----
if (PLOT){
  mapview(stations_utm, color = "navy", col.region = "navy", cex = 2) + 
    mapview(voronoi_utm, alpha.regions = 0.5, col.region = "royalblue")
}

# Summarize poi per Voronoi ----
# Add an explicit Voronoi cell id
voronoi_cells <- voronoi_utm %>%
  mutate(cell_id = row_number())

# Make sure POIs are in the same CRS as the Voronoi polygons
poi_utm <- st_transform(poi, st_crs(voronoi_cells))

# Attach each POI to the Voronoi cell containing it
poi_in_cells <- st_join(
  poi_utm,
  voronoi_cells %>% select(cell_id, station, area_km2),
  join = st_within,
  left = FALSE) %>% 
  left_join(as_tibble(stations_utm), by = "station") %>% 
  mutate(dist_from_station = st_distance(geometry.x, geometry.y, by_element = T))

# For each voronoi cell count the number of POIs for each category. Add also density.
summary_poi = poi_in_cells %>% 
  st_drop_geometry() %>% 
  group_by(cell_id, station, category, area_km2) %>% 
  summarise(count_poi = n(),
            closest_dist = min(dist_from_station),
            .groups = "drop") %>% 
  mutate(dens_poi = count_poi / area_km2)

count_poi_wide = summary_poi %>% 
  pivot_wider(id_cols = station, names_from = category, values_from = c(count_poi), values_fill = 0)

dens_poi_wide = summary_poi %>% 
  pivot_wider(id_cols = station, names_from = category, values_from = c(dens_poi), values_fill = 0)

# For each station find the closest POI for each category, even if it is outside the Voronoi cell
categories = poi_utm$category %>% unique() %>% sort
closest_poi_wide = bind_cols(stations_utm %>% select(station) %>% st_drop_geometry(), 
                             map_dfc(categories, function(cat) {
                               poi_cat = poi_utm %>% filter(category == cat)
                               closest = st_nearest_feature(stations_utm, poi_cat)
                               dist = as.numeric(st_distance(stations_utm, poi_cat[closest, ], by_element = T))
                               tibble(!!cat := dist)
                             }))

stations_map <- stations_utm %>%
  st_transform(4326)

counts_map = voronoi_utm %>% 
  right_join(count_poi_wide, by = "station") %>% 
  st_transform(4326)

density_map = voronoi_utm %>% 
  right_join(dens_poi_wide, by = "station") %>% 
  st_transform(4326)

closest_map = voronoi_utm %>% 
  right_join(closest_poi_wide, by = "station") %>% 
  st_transform(4326)

if (PLOT){
  mapview(
    closest_map,
    zcol = "train",
    alpha.regions = 0.7
  ) +
    mapview(
      stations_map,
      col.regions = "black",
      cex = 3
    )
}



# Unavailability times ----
bike_unavail = bike_history %>% 
  filter(
    !(data.table::as.ITime(commit_at) > data.table::as.ITime("02:00:00") & 
        data.table::as.ITime(commit_at) < data.table::as.ITime("06:00:00"))
  ) %>%
  mutate(
    no_bikes = (bikes == 0), 
    no_stands = (stands == 0)
  )

unavail_time = bike_unavail %>% 
  mutate(time = data.table::as.ITime(commit_at)) %>% 
  pivot_longer(cols = c("no_bikes", "no_stands"), names_to = "unavail_type", values_to = "unavail") %>% 
  filter(unavail) %>% 
  mutate(time_pretty = factor(str_sub(time, 1, 5),  
                              levels = c(t(outer(c("06", "07", "08", "09", 10:23, "00", "01"), 
                                                 c("00", "15", "30", "45"), 
                                                 paste, sep = ":")), 
                                         "02:00")))

# Plot distribution of unavailability times ----
if (PLOT){
  # Define layout constraints
  plots_per_page <- 24
  num_rows <- 4
  num_cols <- 6
  
  # Calculate total number of pages needed for all stations
  total_stations <- length(unique(unavail_time$station))
  total_pages <- ceiling(total_stations / plots_per_page)
  
  # 3. Open the PDF device
  pdf("application/bike-sharing/bike_unavailability_all_stations.pdf", width = 11.5, height = 7)
  
  # 4. Loop through each page and print
  for (p in 1:total_pages) {
    
    cat(p, "\t")
    
    print(
      unavail_time %>% 
        ggplot() +
        geom_bar(aes(x = time_pretty, y = after_stat(prop), fill = unavail_type, group = unavail_type), 
                 position = "identity", alpha = 0.6) +
        
        # Use paginate instead of normal facet_wrap
        ggforce::facet_wrap_paginate(~station, nrow = num_rows, ncol = num_cols, page = p, scales = "fixed") +
        scale_y_continuous(limits = c(0, 0.1)) + 
        scale_x_discrete(breaks = c("06:00", "12:00", "18:00", "00:00")) +
        labs(title = paste("Bike Unavailability over Time (Page", p, "of", total_pages, ")"),
             x = "Time of day", 
             y = "Prop. of unavailability") + 
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
              strip.text = element_text(size = 6),
              legend.position = "bottom") 
    )
  }
  
  # 5. Close the PDF device to save the file
  dev.off()
}


# Remove stations with too few unavailability events (less than 50)
count_unavail = unavail_time |> group_by(station, unavail_type) |> summarise(n = n(), .groups = "drop") |> 
  pivot_wider(id_cols = station, names_from = unavail_type, values_from = n, values_fill = 0) |> 
  mutate(total_unavail = no_bikes + no_stands) |> arrange(desc(total_unavail))

quantile(count_unavail$no_bikes, seq(0, 1, by = 0.1))

stations_to_remove = count_unavail |> filter(no_bikes < 50) |> pull(station)

# Modal unavailability times ----
modal_unavail_times = unavail_time %>% 
  group_by(station, unavail_type) %>% 
  summarise(modal_time = names(which.max(table(time_pretty))), .groups = "drop") %>% 
  pivot_wider(id_cols = station, names_from = unavail_type, values_from = modal_time)

# We convert no_bikes into a factor ordered chronologically from 06:00 to 02:00
unavail_data <- modal_unavail_times %>%  
  left_join(stations_map, by = "station") %>%  
  st_as_sf() %>%  
  st_transform(4326) %>% 
  mutate(no_bikes = factor(no_bikes, 
                           levels = c(t(outer(c("06", "07", "08", "09", 10:23, "00", "01"), 
                                              c("00", "15", "30", "45"), paste, sep = ":")), "02:00")),
         no_stands = factor(no_stands, 
                            levels = c(t(outer(c("06", "07", "08", "09", 10:23, "00", "01"), 
                                               c("00", "15", "30", "45"), paste, sep = ":")), "02:00")))


if (PLOT){
  # Palettes
  library(RColorBrewer)
  time_palette = colorRampPalette(brewer.pal(9, "YlOrRd"))(81)
  names(time_palette) = levels(unavail_data$no_bikes)
  
  # Plot the modal unavailability times for no_bikes
  m = mapview(closest_map, zcol = "metro", 
              col.regions = brewer.pal(9, "Blues"),
              map.types = "Stadia",   # force light basemap
              alpha.regions = 0.7) +
    mapview(unavail_data, zcol = "no_bikes", 
            col.regions = time_palette, 
            alpha.regions = 0.7,
            cex = 2,
            legend = FALSE)
  leaflet::addLegend(
    m@map,
    position = "topleft",
    colors = time_palette[c("06:00", "09:00","12:00", "15:00","18:00","21:00","00:00","02:00")],
    labels = c("06:00", "09:00","12:00", "15:00","18:00","21:00","00:00","02:00"),
    title = "No bikes"
  )
  
  
  # Plot the modal unavailability times for no_stands
  m = mapview(closest_map, zcol = "metro", 
              col.regions = brewer.pal(9, "Blues"),
              map.types = "Stadia",   # force light basemap
              alpha.regions = 0.7) +
    mapview(unavail_data, zcol = "no_stands", 
            col.regions = time_palette, 
            alpha.regions = 0.7,
            cex = 2,
            legend = FALSE)
  leaflet::addLegend(
    m@map,
    position = "topleft",
    colors = time_palette[c("06:00", "09:00","12:00", "15:00","18:00","21:00","00:00","02:00")],
    labels = c("06:00", "09:00","12:00", "15:00","18:00","21:00","00:00","02:00"),
    title = "No stands"
  )
  
  
  # Plot vs population density ----
  m = mapview(voronoi_utm, zcol = "population_density", 
              col.regions = brewer.pal(9, "Blues"),
              map.types = "Stadia",   # force light basemap
              alpha.regions = 0.7) +
    mapview(unavail_data, zcol = "no_stands", 
            col.regions = time_palette, 
            alpha.regions = 0.7,
            cex = 2,
            legend = FALSE)
  leaflet::addLegend(
    m@map,
    position = "topleft",
    colors = time_palette[c("06:00", "09:00","12:00", "15:00","18:00","21:00","00:00","02:00")],
    labels = c("06:00", "09:00","12:00", "15:00","18:00","21:00","00:00","02:00"),
    title = "No stands"
  )
}

# Load weather dataset ----
df_weather_15 = readRDS("data/bike-sharing/milan_weather_ASPRA_15mins.rds")


# Link to every station the weather at the time of the bike history commit
# For every bike station I need to link the closest temperature station and the closest precipitation station
bike_stations = unavail_data |> select(station, geometry) |> distinct()

library(data.table)

# convert once (by reference) if not already a data.table
setDT(df_weather_15)

temp_df <- df_weather_15[Index == "Temperature", ] |> rename(Temperature = Value)
precip_df <- df_weather_15[Index == "Precipitation", ] |> rename(Precipitation = Value)

temp_stations <- unique(temp_df[, .(IdSensor, geometry)], by = "IdSensor") |> st_as_sf()
precip_stations <- unique(precip_df[, .(IdSensor, geometry)], by = "IdSensor") |> st_as_sf()


# Check that the geometries in weather coincides with the geometries in the bike stations
st_crs(temp_stations) == st_crs(bike_stations) # FALSE
st_crs(precip_stations) == st_crs(bike_stations) # FALSE

st_crs(bike_stations)[1]; st_crs(temp_stations)[1]; st_crs(precip_stations)[1]

# To compute distances we transform everything to EPSG:32632
bike_stations_32632 <- st_transform(bike_stations, 32632)
temp_stations_32632 <- st_transform(temp_stations, 32632)
precip_stations_32632 <- st_transform(precip_stations, 32632)

# Create conversion from bike station to closest temp/precip station
bike2weather <- bike_stations_32632 |> 
  mutate(
    closest_temp_idx   = st_nearest_feature(bike_stations_32632, temp_stations_32632),
    closest_precip_idx = st_nearest_feature(bike_stations_32632, precip_stations_32632),
    
    # Extract just the ID vector using the indices
    closest_temp_station   = temp_stations_32632$IdSensor[closest_temp_idx],
    closest_precip_station = precip_stations_32632$IdSensor[closest_precip_idx]
  ) |> 
  select(-closest_temp_idx, -closest_precip_idx) |> st_drop_geometry()

# Create the dataset for the final analysis ----
data_for_regression = unavail_time |> 
  # Remove stations without enough unavailability events
  filter(!station %in% stations_to_remove) |>
  # Join weather (bike2weather conversion, temperature, precipitation)
  left_join(bike2weather, by = "station") |> 
  left_join(temp_df |> select(IdSensor, DateTime, Temperature), by = c("closest_temp_station" = "IdSensor", "commit_at" = "DateTime")) |> 
  left_join(precip_df |> select(IdSensor, DateTime, Precipitation), by = c("closest_precip_station" = "IdSensor", "commit_at" = "DateTime")) |>
  # Join population density
  left_join(voronoi_utm |> st_drop_geometry() |> select(station, population_density), by = "station") |>
  # Join spatial covariates: counts, density, closest distance
  left_join(dens_poi_wide |> select(station, food, supermarkets, offices, nightlife, parks, tourism) |> rename_with(~ paste0("density_", .x), -station), 
            by = "station") |>
  left_join(closest_poi_wide |> select(station, education, metro, train, tram, bus) |> rename_with(~ paste0("closest_", .x), -station), 
            by = "station") |>
  select(-c(skipped_updates, bikes, stands, unavail, 
            time_pretty, closest_temp_station, closest_precip_station)) |> 
  mutate(
    # time_factor = factor(str_sub(time, 1, 5),  
    #                           levels = c(t(outer(c("06", "07", "08", "09", 10:23, "00", "01"), c("00", "15", "30", "45"), paste, sep = ":")), "02:00")),
         time_numeric_shifted = (as.numeric(time) - 6*60*60) %% (24*60*60),  # Shift time so that 06:00 is 0, and wrap around at 24 hours
         .after = time)



# Define objects for MRF: adjacency matrix and penalty matrix
stations_names = sort(unique(data_for_regression$station))
B_MRF = sapply(data_for_regression$station, function(x) 1*(stations_names == x)) |> t()
colnames(B_MRF) = stations_names
# Just a couple of checks
all(rowSums(B_MRF) == 1) # OK
all(colnames(B_MRF) == names(table(data_for_regression$station))); all(colSums(B_MRF) == table(data_for_regression$station)) # OK

# I need also to define the penalty for the 

saveRDS(data_for_regression, "data/bike-sharing/milan_bikemi_data_for_regression.rds")

