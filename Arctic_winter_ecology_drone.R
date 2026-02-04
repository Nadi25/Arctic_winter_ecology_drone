

# Arctic winter ecology ---------------------------------------------------


# Drone -------------------------------------------------------------------


# install packages --------------------------------------------------------
# install.packages(c("sf", "terra", "mapview", "RStoolbox", "ggplot2",
#                    "amt", "spatstat", "adehabitatHR", "caret", "randomForest"))



# load packages -----------------------------------------------------------
library(sf)
library(terra)
library(mapview)
library(RStoolbox)
library(ggplot2)
library(amt)
library(spatstat)
library(adehabitatHR)
library(caret)
library(randomForest)



# Adventdalen -------------------------------------------------------------

# import data --------------------------------------------------------------
insitu_adv <- read.csv("Data/insitu_adv/20260202/20260202_adventalen_insitu_sampling_reindeer-traces.csv")


# add an ID column --------------------------------------------------------
insitu_adv <- insitu_adv |> 
  mutate(id = 1:nrow(insitu_adv))

insitu_adv$id 

# now we need emlid data to spatialize in-situ records
fn_emlid_adv <- list.files("Data/insitu_adv/20260202/emlid/", full.names = TRUE)

# loop through
emlid_adv <- lapply(fn_emlid_adv, st_read)

view(emlid_adv[[1]])


emlid_adv <- do.call(rbind, emlid_adv)

# rename the emlid name column
emlid_adv$emlid_id <- emlid_adv$Name

emlid_adv$Name <- NULL


# merge in-situ data with the geospatial emlid data
insitu_adv <- merge(emlid_adv, insitu_adv, by = "emlid_id", sort = FALSE)

# we want to change the order of the rows
insitu_adv <- insitu_adv[order(insitu_adv$id),]


# make map ----------------------------------------------------------------

mapview(vect(insitu_adv))

plot(insitu_adv["emlid_id"])


# export
# st_write(insitu_adv, "Data/insitu_adv/20260202/20260202_adventdalen_emlid_insitu_sampling_reindeer_traces.gpkg")


# differentiate tracks and craters and clusters
craters <- insitu_adv[insitu_adv$feature_type == "crater",]

craters$feature_radius

# define min crater radius for NAs
min_crater_radius <- 30

craters$feature_radius[is.na(craters$feature_radius)] <- min_crater_radius

craters$feature_radius

# now, we can buffer the crater points by the radius
craters <- st_buffer(craters, craters$feature_radius/100)

# now the geometry changed from point to polygon


# clusters, we have to create polygons from the perimeter points we have smapled
clusters <- insitu_adv[insitu_adv$feature_type == "cluster",]



clusters <- do.call(rbind, lapply(unique(clusters$feature_id), function(id) {
  cl <- clusters[clusters$feature_id == id, ]
  cl1 <- cl[1, ]
  
  if (nrow(cl) > 2) {
    st_geometry(cl1) <- st_cast(
      st_combine(cl[1:nrow(cl), ]),
      "POLYGON"
    )
  } else {
    st_geometry(cl1) <- st_union(
      st_buffer(cl, min_crater_radius / 100)
    )
  }
  
  cl1 <- st_zm(cl1)
  return(cl1)
}))

plot(clusters["emlid"])

mapview(clusters)

# combine craters and clusters
craters <- rbind(craters, clusters)

mapview(vect(craters), zcol = "feature_type")

cc <- st_centroid(craters)

mapview(cc)

craters$Description <- NULL
st_write(craters, "insitu_adv/20260202/20260202_adventdalen_emlid_insitu_sampling_reindeer_traces.gpkg")




# Loading UAV data --------------------------------------------------------

uav_adv_dsm <- rast("Data/practical_day/rs_adv/20240720_adv_UAV_dsm_aoi.tif")
plot(uav_adv_dsm)

uav_adv_ms <- rast("Data/practical_day/rs_adv/20240720_adv_UAV_multi-spectral_aoi.tif")
uav_adv_ms

# bands
# 1: panchromatic
# 2: VIS blue
# 3: VIS green
# 4: VIS red
# 5: NIR
# 6: SWIR
# 7: thermal


# convert centikelvin to C
uav_adv_ms[[7]] <- (uav_adv_ms[[7]] /100) - 273.15

plot(uav_adv_ms[[7]], col = map.pal("magma", 100))

uav_adv_sd <- rast("Data/practical_day/rs_adv/20250411_adv_UAV_snowdepth_aoi.tif")

plot(uav_adv_sd)

# we can use mapview
rvis <- as(uav_adv_sd, "Raster")
insitu_adv$Description <- NULL
mapview(rvis) + mapview(vect(insitu_adv), zcol = "feature_type")



# statistics on insitu data distribution ----------------------------------

# test spatial randomness of our insitu data
study_window <- as.owin(st_bbox(st_buffer(cc, dist = 100)))


coords <- st_coordinates(cc)
# create point pattern
cc_ppp <- ppp(coords[, 1], coords[, 2], window = study_window)

# on the point patternm clark evens test
mn_test <- clarkevans.test(cc_ppp)
print(mn_test)
 # R < 1: samples show a spatial non-random clustered distribution


# we could do further test on spatial distribution here
# like Ripley's K from wich we could derive Kernal density metrics



# preparing a spatial selection model -------------------------------------
res(uav_adv_dsm)
res(uav_adv_ms)

# with resampling we could bring data from different sensor sources
# to the same grid and resolution

# this is rescaling the resolution to match
uav_adv_dsm <- aggregate(uav_adv_dsm, fact = 4)
res(uav_adv_dsm)

uav_adv_ms <- resample(uav_adv_ms, uav_adv_dsm, method = "bilinear")
uav_adv_sd <- resample(uav_adv_sd, uav_adv_dsm, method = "bilinear")



# calculating env.covariates ------------------------------------------------
# derive meaningful features for our model

# derive terrain variables: slope
uav_adv_slope <- terrain(uav_adv_dsm, v = "slope")
plot(uav_adv_slope)


# derive NDVI as a very basic resource distribution map
# or resource availability map


uav_adv_ndvi <- (uav_adv_ms[[5]] - uav_adv_ms[[4]]) /
  (uav_adv_ms[[5]] + uav_adv_ms[[4]]) 

plot(uav_adv_ndvi)

# now compose or ev. covariates stack
# just those variables, that we want our model to see

env_adv <- c(
  uav_adv_ms[[7]], uav_adv_ndvi, uav_adv_sd, uav_adv_dsm, uav_adv_slope
)

names(env_adv) <- c("thermal", "NDVI", "snowdepth", "elev", "slope")

# value range is pretty variable dependant
# we want to naturalize and rescale this for modeling
env_adv_nrom <- normImage(env_adv)
env_adv_resc <- rescaleImage(env_adv_nrom, ymin = 0, ymax = 1)
env_adv_resc

# transfrom coodrinate system
st_crs(env_adv_resc)
st_crs(craters)
craters <- st_transform(craters, st_crs(env_adv_resc))




# preparing model data ----------------------------------------------------
set.seed(429891)

# we want to randomly sample our environmental covariate values for the craters
# = used area

min_samples_crater <- min_crater_radius*2/(round(res(env_adv)[1], 2)*100)
min_samples_crater

used_samples <- st_as_sf(st_sample(
  craters, size = nrow(craters)*min_samples_crater, by_polygon = TRUE
))

# extracting the covariate values from the raster for the used sample location
used_val <- extract(env_adv_resc, used_samples, ID = FALSE)
used_val


# now we want to sample the "background" or the available conditions
aoi_ms <- as.polygons(uav_adv_ms[[1]]> 0, values = F, na.rm = T)
aoi <- st_as_sf(aoi_ms)
plot(aoi)

# assumption: we dont sample available area in 10m vicinity of the craters
avail_area <- st_difference(aoi, st_union(st_buffer(craters, 10)))
plot(avail_area)

n_avail <- nrow(used_samples)
# should be 120, not 300

avail_samples <- st_as_sf(st_sample(
  avail_area, size = n_avail, by_polygon = TRUE
))

mapview(avail_samples)


# extract also from avail__samples to represent the env. conditions
# at not used locations
avail_val <- extract(env_adv_resc, avail_samples, ID = FALSE)
avail_val

# add response variable the model should predict from the covariates
used_val$used <- 1 # response variable is 1 if the env.co. were used
avail_val$used <- 0 # response variable is 0 if the env.co were not used

selection_df <- rbind(used_val, avail_val)
selection_df
View(selection_df)


selection_df <- selection_df[!apply(selection_df, 
                                    MARGIN = 1, function(x) any(is.na(x))),]


# selection model #1 ------------------------------------------------------
model <- glm(
  used ~ thermal + NDVI + snowdepth + elev + slope, 
  data = selection_df,
  family = binomial(link = "logit")
)

summary(model)

# we can use this model to make sa prediction of crater selection for the entire area
pred_sel <- predict(env_adv_resc, model, type = "response")
plot(pred_sel, main = "Predicted forage suitability")
points(craters, pch = 16, col = "red")




# ENDALEN -----------------------------------------------------------------

# load the insitu snow depth data from endalen and the emlid data
# and merge then the same way as before

# insitu_env with geo-spatial 

# import data --------------------------------------------------------------
insitu_end <- read.csv("Data/practical_day/insitu_end/20260203/20260203_endalen_insitu_sampling_snow-depth.csv")


# add an ID column --------------------------------------------------------
insitu_end <- insitu_end |> 
  mutate(id = 1:nrow(insitu_end))

# now we need emlid data to spatialize in-situ records
fn_emlid_end <- list.files("Data/practical_day/insitu_end/20260203/emlid/", full.names = TRUE)

# loop through
emlid_end <- lapply(fn_emlid_end, st_read)

View(emlid_end[[1]])


emlid_end <- do.call(rbind, emlid_end)

# rename the emlid name column
emlid_end$emlid_id <- emlid_end$Name

emlid_end$Name <- NULL


# merge in-situ data with the geospatial emlid data
insitu_end <- merge(emlid_end, insitu_end, by = "emlid_id", sort = FALSE)

# we want to change the order of the rows
insitu_end <- insitu_end[order(insitu_end$id),]

# make map ----------------------------------------------------------------

mapview(vect(insitu_end))

plot(insitu_end["emlid_id"])



# betula nana locations ---------------------------------------------------
insitu_end_betula <- st_read("Data/practical_day/insitu_end/20250909/20250909_endalen_betula_nana.gpkg")
insitu_end_betula <- insitu_end_betula[nchar(insitu_end_betula$Name) <= 2,]




# Loading UAV data --------------------------------------------------------

uav_end_dsm <- rast("Data/practical_day/rs_end/20250911_end_UAV_dsm_aoi.tif")

# transform our insitu data to UTM 33N
st_crs(uav_end_dsm)
insitu_end <- st_transform(insitu_end, st_crs(uav_end_dsm))
insitu_end_betula <- st_transform(insitu_end_betula, st_crs(uav_end_dsm))

mapview(vect(insitu_end_betula))

focus <- st_as_sfc(st_bbox(st_buffer(insitu_end_betula[insitu_end_betula$Name == "26", ], 
                                     90)))

mapview(vect(insitu_end_betula)) + mapview(vect(focus))+
  mapview(vect(insitu_end))


# derive terrain variables
uav_end_slope <- terrain(uav_end_dsm, v = "slope")
uav_end_aspect <- terrain(uav_end_dsm, v = "aspect")
uav_end_flowdir <- terrain(uav_end_dsm, v = "flowdir")
uav_end_tpi <- terrain(aggregate(uav_end_dsm, 18), v = "TPI")
uav_end_tri <- terrain(aggregate(uav_end_dsm, 18), v = "TRI")

uav_end_tpi <- resample(uav_end_tpi, uav_end_dsm)
uav_end_tri <- resample(uav_end_tri, uav_end_dsm)

# mow we can stack this
env_end <- c(
  uav_end_dsm, uav_end_slope, uav_end_aspect, uav_end_flowdir,
  uav_end_tpi, uav_end_tri
)
names(env_end) <- c("elev", "slope", "aspect", "flowdir", "tpi", "tri")

# norm and rescaling
env_end_norm <- normImage(env_end)
env_end_resc <- rescaleImage(env_end_norm, ymin = 0, ymax = 0)

# clean our the NA sample
insitu_end <- insitu_end[!is.na(insitu_end$snow_depth),]

# we extract env covariates at the snow depth sampling locations
samples_df <- extract(env_end_resc, insitu_end, ID = FALSE)
samples_df$snow_depth <- insitu_end$snow_depth/max(insitu_end$snow_depth)


# snow depth model --------------------------------------------------------

x <- samples_df[, colnames(samples_df) != "snow_depth"]
y <- samples_df$snow_depth

# train a model: randomForest regression
model <- train(
  x = x, y = y,
  trControl = trainControl(
    p = 0.75,
    method = "cv",
    number = 5,
    verboseIter = TRUE,
    #classProbs = TRUE # active if you want to use radndomforest ML model for presence
    # absence or use/availability modelling
  )
)

model

# predict snow depth for focus area
env_end_resc_crop <- crop(env_end_resc, focus)
pred_snow_depth <- predict(env_end_resc_crop, model, na.rm = T)

plot(pred_snow_depth[[1]] * max(insitu_end$snow_depth),
     main = "Snow depth predicted from topography (cm)")














