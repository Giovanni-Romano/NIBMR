# Generic pkgs ----
suppressPackageStartupMessages(library(tidyverse))
library(sf)
library(data.table)

# Bike sharing data ----
## DB setup for bike sharing data ----
library(DBI)
library(duckdb)

con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

dbExecute(con, "INSTALL httpfs;")
dbExecute(con, "LOAD httpfs;")
dbExecute(con, "SET s3_endpoint='storage.googleapis.com';")

## Bikes ----
dbExecute(con, "
  COPY (
    SELECT *
    FROM read_parquet('s3://bike-sharing-history/milan/bikemi/*/*.parquet')
  )
  TO 'data/bike-sharing/milan_bikemi_history_original.parquet'
  (FORMAT PARQUET);
")

## Weather forecast (not used) ----
dbExecute(con, "
  COPY (
    SELECT *
    FROM read_parquet('s3://weather-forecast-history/milan/*/*.parquet')
  )
  TO 'data/bike-sharing/milan_weather_history.parquet'
  (FORMAT PARQUET);
")


# Download weather data from ARPA: https://www.arpalombardia.it/temi-ambientali/meteo-e-clima/form-richiesta-dati/

# Weather data from ARPA (load and clean) ----
files_weather = list.files("data/bike-sharing/weather/ARPA/csv", pattern = "csv$", full.names = TRUE)
list_weather = lapply(files_weather, read_csv, show_col_types = FALSE)
for (i in seq_along(list_weather)){
  colnames(list_weather[[i]])[1:2] = c("IdSensor", "DateTime")
  
  if (colnames(list_weather[[i]])[3] == "Medio"){
    list_weather[[i]]$Index = "Temperature"
  } else {
    list_weather[[i]]$Index = "Precipitation"
  }
  
  colnames(list_weather[[i]])[3] = "Value"
}
df_weather = bind_rows(list_weather)


# Check if there are duplicates
df_weather %>% 
  count(IdSensor, DateTime, Index) %>% 
  filter(n > 1)
# Some duplicates (on 30/12 and 31/12)

# Remove duplicates
df_weather = df_weather %>% distinct(IdSensor, DateTime, Index, .keep_all = TRUE)
df_weather %>% 
  count(IdSensor, DateTime, Index) %>% 
  filter(n > 1)
# No more duplicates

# NAs are denoted with -999, replace them with NA
df_weather$Value[df_weather$Value == -999] = NA
sum(is.na(df_weather$Value)) # 15606
df_weather |> group_by(IdSensor) |> summarize(n_NA = sum(is.na(Value)))

stations_weather = read_csv("data/bike-sharing/weather/ARPA/stazioni.csv", show_col_types = FALSE) %>% 
  filter(Provincia == "MI", Tipologia %in% c("Temperatura", "Precipitazione")) %>% 
  select(IdSensore, Tipologia, NomeStazione, Quota, lng, lat, UTM_Nord, UTM_Est) %>% 
  rename(IdSensor = IdSensore, Index = Tipologia, 
         StationName = NomeStazione, Altitude = Quota, Longitude = lng, Latitude = lat) %>% 
  mutate(Index = ifelse(Index == "Temperatura", "Temperature", "Precipitation"))

# Are there in stations_weather all the stations in df_weather?
setdiff(unique(df_weather$IdSensor), unique(stations_weather$IdSensor)) # Yes, good

df_weather = df_weather |> left_join(stations_weather |> select(IdSensor, Index, UTM_Nord, UTM_Est), by = c("IdSensor", "Index"))

# See temporal distribution of NAs for every station
df_weather |>
  filter(is.na(Value)) |>
  mutate(Date = as.Date(DateTime)) |>
  count(IdSensor, Date) |>
  ggplot(aes(x = Date, y = n)) +
  geom_col() +
  facet_wrap(~IdSensor, scales = "fixed")
# It looks like there are some combinations of stations and day for which there
# are NA for the whole day
# I'll impute temperature/precipitation in a station as a weighted avarage of the 
# other stations, with weights inversely proportional to the distance between stations. 

# Split weather data into temperature and precipitation
temp_df = df_weather |> filter(Index == "Temperature") |> st_as_sf(coords = c("UTM_Est", "UTM_Nord"), crs = 32632, remove = TRUE)
prec_df = df_weather |> filter(Index == "Precipitation") |> st_as_sf(coords = c("UTM_Est", "UTM_Nord"), crs = 32632, remove = TRUE)

temp_df_NA = temp_df |> filter(is.na(Value))
temp_df_notNA = temp_df |> filter(!is.na(Value))
prec_df_NA = prec_df |> filter(is.na(Value))
prec_df_notNA = prec_df |> filter(!is.na(Value))



# Precompute distance matrix between unique sensor locations
temp_coords <- stations_weather |> filter(IdSensor %in% unique(temp_df$IdSensor)) |> st_as_sf(coords = c("UTM_Est", "UTM_Nord"), crs = 32632, remove = TRUE)
temp_dist <- st_distance(temp_coords)
dimnames(temp_dist) <- list(temp_coords$IdSensor, temp_coords$IdSensor)
# Move to data.table
setDT(temp_df_notNA)
setDT(temp_df_NA)
temp_df_NA[, row_id := .I]
# Join on DateTime only. allow.cartesian=TRUE allows for multiple matches (i.e., 
# multiple stations with the same DateTime) so that we match each station in NA 
# to every station in notNA.
temp_joined <- temp_df_notNA[temp_df_NA, on = "DateTime", allow.cartesian = TRUE]
# Drop self-matches, look up distance via matrix indexing (fast, no geometry calc)
temp_joined <- temp_joined[IdSensor != i.IdSensor]
temp_joined[, dist := temp_dist[cbind(as.character(i.IdSensor), as.character(IdSensor))]]
temp_joined[, weight := 1 / dist]
# Aggregate weighted average per original NA row
temp_result <- temp_joined[, .(imputed = sum(Value * weight) / sum(weight)), by = row_id]
# Write back
temp_df_NA[temp_result, on = "row_id", Value := i.imputed]

# Precompute distance matrix between unique sensor locations (once, not per row)
prec_coords <- stations_weather |> filter(IdSensor %in% unique(prec_df$IdSensor)) |> st_as_sf(coords = c("UTM_Est", "UTM_Nord"), crs = 32632, remove = TRUE)
prec_dist <- st_distance(prec_coords)
dimnames(prec_dist) <- list(prec_coords$IdSensor, prec_coords$IdSensor)
# Move to data.table
setDT(prec_df_notNA)
setDT(prec_df_NA)
prec_df_NA[, row_id := .I]
# Cartesian join (see above for temp)
prec_joined <- prec_df_notNA[prec_df_NA, on = "DateTime", allow.cartesian = TRUE]
# Drop self-matches, look up distance via matrix indexing
prec_joined <- prec_joined[IdSensor != i.IdSensor]
prec_joined[, dist := prec_dist[cbind(as.character(i.IdSensor), as.character(IdSensor))]]
prec_joined[, weight := 1 / dist]
# Aggregate weighted average per original NA row
prec_result <- prec_joined[, .(imputed = sum(Value * weight) / sum(weight)), by = row_id]
# Write back
prec_df_NA[prec_result, on = "row_id", Value := i.imputed]

# Re-build weather data.frame
df_weather = bind_rows(
  temp_df_notNA |> st_drop_geometry(),
  temp_df_NA |> st_drop_geometry(),
  prec_df_notNA |> st_drop_geometry(),
  prec_df_NA |> st_drop_geometry()
) |> select(-row_id) |> arrange(IdSensor, DateTime, Index)
# Check no NA
df_weather$Value |> is.na() |> sum() # 0, great
sum(df_weather$Value == -999) # 0, great

# Transform date and time from char
df_weather$DateTime_UTCp1 = as.POSIXct(df_weather$DateTime, format = "%Y/%m/%d %H:%M", tz = "Etc/GMT-1") # Times are in UTC+1, so Etc/GMT-1
sum(is.na(df_weather$DateTime_UTCp1))
df_weather$DateTime = as.POSIXct(df_weather$DateTime_UTCp1, tzone = "Europe/Rome") # Shift to Europe/Rome timezon
sum(is.na(df_weather$DateTime))


# I need weather every 15 mins but dataset has every 10
# :00 and :30 --> keep as it is
df_weather_exact <- df_weather %>%
  filter(minute(DateTime) %in% c(0, 30))

# Use data.table because it's way faster than dplyr for this kind of operation
library(data.table)
dt <- as.data.table(df_weather)

# :15 -> average of :10 and :20
df_weather_at15 <- dt[minute(DateTime) %in% c(10, 20), .(Value = mean(Value, na.rm = TRUE)), 
                      by = .(IdSensor, hour_block = DateTime - minutes(minute(DateTime)))][, DateTime := hour_block + minutes(15)][, hour_block := NULL]

# :45 -> average of :40 and :50
df_weather_at45 <- dt[minute(DateTime) %in% c(40, 50), .(Value = mean(Value, na.rm = TRUE)), 
                      by = .(IdSensor, hour_block = DateTime - minutes(minute(DateTime)))][, DateTime := hour_block + minutes(45)][, hour_block := NULL]



# Transform UTM_Nord and UTM_Est in a geometry column
df_weather_15mins = bind_rows(
  df_weather_exact,
  merge(df_weather_at15, stations_weather |> select(-StationName, -Altitude, -Longitude, -Latitude), by = "IdSensor") |> st_as_sf(coords = c("UTM_Est", "UTM_Nord"), crs = 32632, remove = FALSE),
  merge(df_weather_at45, stations_weather |> select(-StationName, -Altitude, -Longitude, -Latitude), by = "IdSensor") |> st_as_sf(coords = c("UTM_Est", "UTM_Nord"), crs = 32632, remove = FALSE)
)


# Save
saveRDS(df_weather, "data/bike-sharing/milan_weather_ASPRA_original.rds")
saveRDS(df_weather_15mins, "data/bike-sharing/milan_weather_ASPRA_15mins.rds")

# OSM data ----
library(osmdata)
# Set overpass API to perform queries
osmdata::set_overpass_url("https://z.overpass-api.de/api/interpreter")
# Options
#   1. "https://overpass-api.de/api/interpreter"
#   2. "https://overpass.kumi.systems/api/interpreter"
#   3. "https://lz4.overpass-api.de/api/interpreter"
#   4. "https://z.overpass-api.de/api/interpreter"

# Set bounding box for Milan
bb = getbb("Milano, Lombardia, Italy")

# Helper to fix same bounding box
osm_query = function(...) opq(bbox = bb, timeout = 180) %>% add_osm_feature(...)

# Helper to convert polygons and multipolygons to points
osm_as_points = function(osm_result, keep_cols = c("osm_id", "name", "amenity", "geometry")) {
  list(osm_result$osm_points, osm_result$osm_polygons, osm_result$osm_multipolygons) %>%
    # drop NULL or empty layers
    keep(\(x) !is.null(x) && nrow(x) > 0) %>%
    # convert polygon layers to centroids
    map(\(x) if (inherits(st_geometry(x), "sfc_POINT")) x else st_centroid(x)) %>%
    # harmonise columns before binding (mismatched OSM tags would break bind_rows)
    map(\(x) select(x, any_of(keep_cols))) %>%
    bind_rows() %>% 
    distinct(name, .keep_all = T)
}


# Education
# edu = osm_query(key = "amenity",
#                 value = c("school", "university", "college", "kindergarten")) %>%
#   osmdata_sf() %>% osm_as_points() %>% filter(!is.na(name))

# school = osm_query(key = "amenity", value = "school") %>%
#   osmdata_sf() %>% osm_as_points() %>% filter(!is.na(name))

# kindergar = osm_query(key = "amenity", value = "kindergarten") %>%
#   osmdata_sf() %>% osm_as_points() %>% filter(!is.na(name))


uni = osm_query(key = "amenity", value = "university") %>%
  osmdata_sf() %>% osm_as_points() %>% filter(!is.na(name))

college = osm_query(key = "amenity", value = "college") %>%
  osmdata_sf() %>% osm_as_points() %>% filter(!is.na(name))

edu = bind_rows(uni, college) %>% distinct(name, .keep_all = T)

# Food & drink (lunch + evening peaks)
food = osm_query(key = "amenity",
                 value = c("restaurant", "bar", "cafe", "fast_food")) %>%
  osmdata_sf() %>% osm_as_points() %>% filter(!is.na(name))

# Supermarkets
supermarkets = osm_query(key = "shop",
                         value = c("supermarket", "convenience")) %>%
  osmdata_sf() %>% osm_as_points() %>% filter(!is.na(name))

# Public transport
metro = osm_query(key = "railway",
                  value = c("subway_entrance")) %>%
  osmdata_sf() %>% osm_as_points() %>% filter(!is.na(name))

train = osm_query(key = "railway",
                  value = c("station")) %>%
  osmdata_sf() %>% osm_as_points() %>% filter(!is.na(name))

tram = osm_query(key = "railway",
                 value = c("tram_stop")) %>%
  osmdata_sf() %>% osm_as_points() %>% filter(!is.na(name))

bus = osm_query(key = "highway", value = "bus_stop") %>%
  osmdata_sf() %>% osm_as_points() %>% filter(!is.na(name))

# Offices / workplaces
offices = osm_query(key = "office") %>%   # no value = all office types
  osmdata_sf() %>% osm_as_points() %>% filter(!is.na(name))

# Nightlife & entertainment (late-evening peaks)
nightlife = osm_query(key = "amenity",
                      value = c("nightclub", "bar", "pub", "theatre", "cinema")) %>%
  osmdata_sf() %>% osm_as_points() %>% filter(!is.na(name))

# Parks & green spaces (afternoon/weekend leisure)
parks = osm_query(key = "leisure",
                  value = c("park", "garden", "playground", "sports_centre")) %>%
  osmdata_sf() %>% osm_as_points() %>% filter(!is.na(name))

# Healthcare
health = osm_query(key = "amenity",
                   value = c("hospital", "clinic", "pharmacy")) %>%
  osmdata_sf() %>% osm_as_points() %>% filter(!is.na(name))

# Tourism (daytime / weekend peaks)
tourism = osm_query(key = "tourism",
                    value = c("museum", "attraction", "hotel", "gallery")) %>%
  osmdata_sf() %>% osm_as_points() %>% filter(!is.na(name))


# Put everything together in a data frame
poi_df = bind_rows(
  education = edu,
  food = food,
  supermarkets = supermarkets,
  metro = metro,
  train = train,
  tram = tram,
  bus = bus,
  offices = offices,
  nightlife = nightlife,
  parks = parks,
  health = health,
  tourism = tourism,
  .id = "category"
)

saveRDS(poi_df, "data/bike-sharing/milan_poi.rds")


# Make Voronoi of stations ----
library(arrow)
bike_events_original = read_parquet("data/bike-sharing/milan_bikemi_history_original.parquet")

# Commits are made only if there are changes in the number of bikes at a station. I want to create a regular 15 minutes grid.
bike_events <- bike_events_original %>%
  # # First tell R the raw numbers are UTC
  # # Then shift the clock to Europe/Berlin time (accounts for German daylight savings)
  # mutate(commit_at = with_tz(force_tz(commit_at, tzone = "UTC"), tzone = "Europe/Berlin")) %>% 
  arrange(station, latitude, longitude, commit_at) %>%
  group_by(station, latitude, longitude) %>%
  complete(
    commit_at = seq(
      ceiling_date(min(commit_at), "15 minutes"),
      ceiling_date(max(commit_at), "15 minutes"),
      by = "15 min"
    )
  ) %>%
  arrange(commit_at) %>% 
  fill(station, longitude, latitude, bikes, stands, .direction = "down") %>% 
  ungroup() %>% 
  filter((minute(commit_at) %in% c(0, 15, 30, 45) & second(commit_at) == 0 & (hour(commit_at) >= 6 | hour(commit_at) < 02)) | 
           (minute(commit_at) == 0 & second(commit_at) == 0 & hour(commit_at) == 2))


stations_df = bike_events %>% distinct(station, latitude, longitude) %>%
  filter(!is.na(latitude) & !is.na(longitude))

events_per_station = count(bike_events, station, latitude, longitude, name = "n")

stations_sf = st_as_sf(stations_df %>% left_join(events_per_station, by = c("station", "longitude", "latitude")), 
                       coords = c("longitude", "latitude"),
                       crs = 4326, # specify that the coordinates are in lon/lat format
                       remove = FALSE)

# Clean from errors
to_remove = stations_sf %>% 
  filter(str_detect(station, fixed("test")) | str_detect(station, fixed("Test ")) |
           floor(latitude) != 45 | floor(longitude) != 9) %>% # Remove test stations and outliers
  mutate(to_remove = T)

stations_sf = left_join(stations_sf, to_remove %>% select(station, latitude, longitude, to_remove) %>% st_drop_geometry(), by = c("station", "latitude", "longitude")) %>%
  filter(is.na(to_remove)) %>% 
  select(-to_remove)

bike_events = left_join(bike_events, to_remove %>% select(station, latitude, longitude, to_remove) %>% st_drop_geometry(), by = c("station", "latitude", "longitude")) %>%
  filter(is.na(to_remove)) %>% 
  select(-to_remove)


# Check if there are stations with different coordinates but same name
multiple_loc = stations_sf %>% group_by(station) %>% summarise(n_lat = n_distinct(latitude), n_long = n_distinct(longitude)) %>% filter(n_lat > 1 | n_long > 1)
#   station                             n_lat n_long
#   <chr>                               <int>  <int>
# 1 19 - Italia - Santa Sofia            2      2
# 3 61 - Augusto                         2      2
# 4 V1 - Apulejo - Campania              2      2
# --- --- --- --- --- --- --- --- ---

# Check if there are stations with same coordinates but different name
stations_sf %>% group_by(latitude, longitude) %>% summarise(n_stations = n_distinct(station)) %>% filter(n_stations > 1)
# There aren't any
# --- --- --- --- --- --- --- --- ---

# Check multiple stations by ID
multiple_names = stations_sf %>% mutate(ID = str_split_i(station, " - ", 1)) %>% group_by(ID) %>% summarise(n = n()) %>% filter(n > 1) %>% arrange(ID)
#   ID        n
#   <chr>   <int>
# 1 166       2
# 2 19        2
# 3 274       2
# 4 352       2
# 5 61        2
# 6 V1        2

summary_multiple = bike_events %>% 
  mutate(ID = str_split_i(station, " - ", 1), .after = station) %>% 
  filter(ID %in% multiple_names$ID) %>% 
  group_by(station, ID, latitude, longitude) %>%
  summarise(n = n(),
            min_date = min(commit_at),
            max_date = max(commit_at)) %>%
  arrange(station)

#   station                         ID    latitude longitude     n min_date            max_date           
#   <chr>                           <chr>    <dbl>     <dbl> <int> <dttm>              <dttm>             
# 1 166 - Dezza - Caravaggio        166       45.5      9.16  5028 2024-04-10 23:45:34 2025-06-01 01:06:35
# 2 166 - Foppa/Coni Zugna          166       45.5      9.17  4784 2025-10-28 17:08:29 2026-05-31 22:14:18
# 3 19 - Italia - Santa Sofia       19        45.5      9.19 10313 2024-04-10 23:45:34 2025-12-19 07:45:26
# 4 19 - Italia - Santa Sofia       19        45.5      9.19  4050 2025-12-19 08:37:07 2026-06-01 01:09:20
# 5 274 - Belloveso - Monte Rotondo 274       45.5      9.19  2727 2024-04-10 23:45:34 2025-02-20 12:49:12
# 6 274 - Passerini Belloveso       274       45.5      9.19  5746 2025-02-20 13:09:54 2026-06-01 00:08:35
# 7 352 - Cordusio                  352       45.5      9.19  9578 2024-04-10 23:45:34 2025-05-04 23:23:26
# 8 352 - Cordusio - Dante          352       45.5      9.19  3941 2025-12-05 16:07:16 2026-05-31 21:44:39
# 10 61 - Augusto                    61       45.5      9.20 14687 2024-08-22 16:48:42 2026-05-31 23:36:47
# 11 61 - Augusto                    61       45.5      9.20    88 2024-04-10 23:45:34 2025-01-05 08:45:50
# 12 V1 - Apulejo - Campania         V1       45.5      9.22  9243 2024-04-10 23:45:34 2026-05-29 09:54:48
# 13 V1 - Apulejo - Campania         V1       45.5      9.22    18 2026-05-29 15:29:00 2026-05-31 21:44:39

# 166: they probably moved the station -> merge
# 19: they probably moved the station -> merge
# 274: they probably moved the station -> merge
# 352: they probably moved the station -> merge
# 61: temporary move? -> merge
# V1: they propably moved the station -> merge

# Merge stations above
bike_events_merged = bike_events %>% mutate(ID = str_split_i(station, " -", 1), .after = station)
stations_sf_merged = stations_sf %>% mutate(ID = str_split_i(station, " -", 1), .after = station)

for (id in multiple_names$ID){
  stations_to_merge = summary_multiple %>% filter(ID == id) %>% arrange(min_date)
  
  st1 = stations_to_merge[1, ]
  st2 = stations_to_merge[2, ]
  
  if (st1$min_date < st2$min_date){
    old = st1$station
    new = st2$station
  }
  
  geom_tmp = st_centroid(st_combine(st_sfc(st_point(c(stations_to_merge$longitude[1], stations_to_merge$latitude[1])), 
                                           st_point(c(stations_to_merge$longitude[2], stations_to_merge$latitude[2])))))
  
  bike_events_merged$latitude[bike_events$station %in% c(old, new)] = st_coordinates(geom_tmp)[, 2]
  bike_events_merged$longitude[bike_events$station %in% c(old, new)] = st_coordinates(geom_tmp)[, 1]
  bike_events_merged$station[bike_events$station %in% c(old, new)] = new
  
  stations_sf_merged$geometry[stations_sf$station %in% c(old, new)] = geom_tmp
  stations_sf_merged$latitude[stations_sf$station %in% c(old, new)] = st_coordinates(geom_tmp)[, 2]
  stations_sf_merged$longitude[stations_sf$station %in% c(old, new)] = st_coordinates(geom_tmp)[, 1]
  stations_sf_merged$station[stations_sf$station %in% c(old, new)] = new
  stations_sf_merged$n[stations_sf$station %in% c(old, new)] = sum(stations_sf$n[stations_sf$station %in% c(old, new)])
}

stations_sf_merged = stations_sf_merged %>% distinct()

# Check anomalies
stations_sf_merged %>% group_by(station) %>% summarise(n = n()) %>% arrange(desc(n))
# No more duplicates

# There are are locations with multiple stations really close to each other, e.g. Cordusio and Cadorna
# It could be worth to merge them in just one, so that the voronoi cell is more meaningful
stations_dist = st_distance(stations_sf_merged) %>% as.numeric %>% matrix(nrow = nrow(stations_sf_merged), ncol = nrow(stations_sf_merged), 
                                                                          dimnames = list(stations_sf_merged$station, stations_sf_merged$station))

stations_dist_df = stations_dist %>% reshape2::melt(varnames = c("from", "to"), value.name = "distance_m") %>% 
  filter(from != to) %>% 
  mutate(from_id = str_split_i(from, " - ", 1),
         to_id = str_split_i(to, " - ", 1),
         from_to = paste0(pmin(from_id, to_id), "-", pmax(from_id, to_id))) %>% 
  distinct(from_to, .keep_all = T)
stations_dist_df %>% arrange(distance_m)

too_close = stations_dist_df %>% filter(distance_m <= 100) %>% select(from, to)

# Check that no station occurs in multiple pairs to merge (otherwise we would need to merge more than 2 stations together)
table(as.character(c(too_close$from, too_close$to))) %>% sort
# none

stations_sf_merged_twice = stations_sf_merged
bike_events_merged_twice = bike_events_merged

for (i in 1:nrow(too_close)){
  p_from = st_point(stations_sf_merged %>% filter(station == too_close$from[i]) %>% select(geometry) %>% st_coordinates())
  p_to = st_point(stations_sf_merged %>% filter(station == too_close$to[i]) %>% select(geometry) %>% st_coordinates())
  
  points = st_sfc(p_from, p_to, crs = 4326)
  barycenter = st_centroid(st_combine(points))
  
  # Remove the "from" station and update the "to" station with the new geometry and name
  stations_sf_merged_twice[stations_sf_merged_twice$station == too_close$to[i], "geometry"] = st_sf(barycenter)
  stations_sf_merged_twice[stations_sf_merged_twice$station == too_close$to[i], "latitude"] = st_coordinates(barycenter)[, 2]
  stations_sf_merged_twice[stations_sf_merged_twice$station == too_close$to[i], "longitude"] = st_coordinates(barycenter)[, 1]
  stations_sf_merged_twice[stations_sf_merged_twice$station == too_close$to[i], "n"] = sum(stations_sf_merged$n[stations_sf_merged$station %in% c(too_close$from[i], too_close$to[i])])
  stations_sf_merged_twice[stations_sf_merged_twice$station == too_close$to[i], "station"] = paste0(too_close$to[i], " + ", too_close$from[i])
  stations_sf_merged_twice = stations_sf_merged_twice[stations_sf_merged_twice$station != too_close$from[i], ]
  
  # Same for bike events
  bike_events_merged_twice$latitude[bike_events_merged_twice$station %in% c(too_close$from[i], too_close$to[i])] = st_coordinates(barycenter)[, 2]
  bike_events_merged_twice$longitude[bike_events_merged_twice$station %in% c(too_close$from[i], too_close$to[i])] = st_coordinates(barycenter)[, 1]
  bike_events_merged_twice$station[bike_events_merged_twice$station %in% c(too_close$from[i], too_close$to[i])] = paste0(too_close$to[i], " + ", too_close$from[i])
}

# Last check on stations with few data
nmbr_events = bike_events_merged_twice %>% group_by(station) %>% summarise(n = n()) %>% arrange(n)
all((nmbr_events %>% arrange(station) %>% pull(station)) == (stations_sf_merged_twice %>% arrange(station) %>% pull(station)))
all((nmbr_events %>% arrange(station) %>% pull(n)) == (stations_sf_merged_twice %>% arrange(station) %>% pull(n)))
stations_sf_merged_twice %>% arrange(n)

saveRDS(bike_events_merged_twice, "data/bike-sharing/milan_bikemi_history.rds")

# Milan boundaries
q = opq(bbox = bb,  timeout = 60)  %>% 
  # Look for administrative boundaries
  add_osm_feature(key = "boundary", value = "administrative") %>% 
  # admin_level=8 means municipalities
  add_osm_feature(key = "admin_level", value = "8")
candidate_boundary = osmdata_sf(q)
boundary_sf = candidate_boundary %>% 
  # Extract the sf object with the boundaries
  pluck("osm_multipolygons") %>% 
  # Select only city of Milan
  filter(name == "Milano")

utm_crs = 32632 # UTM (Universal Transverse Mercator) zone 32N, which covers Milan
stations_utm = st_transform(stations_sf_merged_twice, utm_crs); saveRDS(stations_utm, "data/bike-sharing/milan_stations_utm.rds")

boundary_utm = st_transform(st_make_valid(boundary_sf), utm_crs)
boundary_union = st_union(st_geometry(boundary_utm))

vor = st_voronoi(st_union(st_geometry(stations_utm)),
                 envelope = st_as_sfc(st_bbox(boundary_utm))) %>%
  st_collection_extract("POLYGON") %>% 
  st_cast("POLYGON") %>%
  st_as_sf() %>%
  st_set_crs(utm_crs) %>% 
  # Clip to your actual envelope
  st_intersection(st_sf(geometry = boundary_union))

buff = st_buffer(stations_utm, dist = 500) # 100 meters buffer
mapview(vor) + mapview(buff)

vor_int_buff = st_intersection(vor, st_union(st_geometry(buff)))

mapview(boundary_utm, col.regions = "red3", alpha.regions = 0.1) + 
  mapview(vor_int_buff, alpha.regions = 0.75) + 
  mapview(stations_utm, cex = 2, col.regions = "black")

# Match each Voronoi polygon to its generating station by point intersection.
point_poly = st_intersects(stations_utm, vor_int_buff)
poly_id = sapply(point_poly, function(z) if (length(z) == 0) NA_integer_ else z[1])
keep = !is.na(poly_id)

vor_geom = vor_int_buff[poly_id[keep], , drop = FALSE]
attrs = st_drop_geometry(stations_utm[keep, ])

vor <- st_sf(
  attrs,
  geometry = st_geometry(vor_geom),
  crs = st_crs(stations_utm)
)

vor$area_km2 <- as.numeric(st_area(vor)) / 1e6


mapview(vor)

saveRDS(vor, "data/bike-sharing/milan_voronoi.rds")
