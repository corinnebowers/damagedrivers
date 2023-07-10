
#### packages #####################################################################################

suppressPackageStartupMessages({
  ##general/other
  require(censusapi)
  require(lubridate)
  # require(ncdf4)
  require(caret)
  require(pracma)
  require(readxl)
  require(httr)
  require(jsonlite)
  require(curl)
  ##parallel
  require(foreach)
  require(doSNOW)
  require(parallel)
  ##spatial data
  require(sf)
  require(terra)
  require(raster)
  require(tigris) 
  require(exactextractr)
  ##datavis
  require(ggplot2)
  require(GGally)
  require(scales) 
  require(ggspatial)
  require(cowplot)
  require(ggridges)
  # require(leaflet)
  # require(mapboxapi)
  # require(mapview)
  require(scico)
  require(extrafont)
  require(patchwork)
  require(grid)
  require(gridExtra)
  require(ggpubr)
  require(gt)
  require(ggnewscale)
  ##tidyverse
  require(tidyverse)
})

options(tigris_use_cache = TRUE) #tigris
Sys.setenv(CENSUS_KEY = 'f2e090156b02ced027d4ed756f82c9a3a1aa38c9') #censusapi


#### constants & utility functions ################################################################

## set parallel backend
cores <- 5
opts <- list(progress = function(n) setTxtProgressBar(pb, n))

## useful helper functions

add_index <- function(n) (1:n)/(n+1)

toNumber <- function(x) as.numeric(paste(x))
strip <- function(x) unname(unlist(x))

Mean <- function(x) ifelse(all(is.na(x)), NA, mean(x, na.rm = TRUE))
Sum <- function(x) ifelse(all(is.na(x)), NA, sum(x, na.rm = TRUE))
Min <- function(x) ifelse(all(is.na(x)), NA, min(x, na.rm = TRUE))
Max <- function(x) ifelse(all(is.na(x)), NA, max(x, na.rm = TRUE))
mode <- function(x) names(which.max(table(x)))[1]

setNA <- function(x, new) ifelse(is.na(x), new, x)
sum.na <- function(x) sum(is.na(x))

wateryear <- function(d) year(d) + ifelse(month(d) %in% 10:12, 1, 0)

positive <- function(x) {
  x[x < 0] <- 0
  x
}


#### themes & axes ################################################################################

scale_x_origin <- function(...) {
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0,NA), ...) }
scale_y_origin <- function(...) {
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0,NA), ...) }
geom_parity <- function() geom_abline(slope = 1, intercept = 0, linetype = 'dashed')

theme_set(
  theme_classic() + theme(
    text = element_text(family = 'Segoe UI', size = 9),
    axis.line = element_line(size = 0.35),
    axis.ticks = element_line(size = 0.35, color = 'black'),
    legend.key.size = unit(0.35, 'cm')))
theme_bw_custom <- 
  function() {
    theme_bw() + theme(
      text = element_text(family = 'Segoe UI'),
      panel.border = element_rect(fill = NA, color = 'black', size = 0.35),
      axis.ticks = element_line(size = 0.35, color = 'black'),
      legend.key.size = unit(0.35, 'cm'))}


#### colors #######################################################################################

baker <- c()
baker[1] <- rgb(56, 95, 150, maxColorValue = 255)
baker[2] <- rgb(207, 89, 33, maxColorValue = 255)
baker[3] <- rgb(158, 184, 219, maxColorValue = 255)
baker[4] <- rgb(231, 184, 0, maxColorValue = 255)
baker[5] <- rgb(128, 0, 0, maxColorValue = 255)

ggcolor <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

roma.colors <- scico(n = 5, palette = 'roma')


#### geospatial ###################################################################################

## function to plot rasters in ggplot
raster.df <- function(x) x %>% as.data.frame(xy = TRUE) %>% setNames(c('x', 'y', 'value'))

## define coordinate systems
wgs84 <- 4326
nad83 <- 4269
projected <- 6417

## define useful spatial geometries
california <- counties(state = 'CA', year = 2020)
caltracts <- tracts(state = 'CA', year = 2020)
cal <- st_union(california)

## define global raster grid 
xres <- 0.625
yres <- 0.5
grid <- 
  raster(
    xmn = -180, xmx = 180, ymn = -90, ymx = 90, 
    resolution = c(xres,yres), crs = projection(california))

## identify grid cells that cover california
grid_ca <- grid %>% 
  raster::crop(california, snap = 'out') %>% 
  setValues(1:ncell(.))
grid_ca[] <- ifelse(rasterize(california, grid_ca, getCover = TRUE)[]>=0.05, grid_ca[], NA)
index_ca <- grid_ca %>% raster.df %>% filter(!is.na(value)) %>% pull(value)


#### functions: AR catalog ########################################################################

calculate_duration <- function(df, int = 3) {
  ## calculates AR duration
  #' df: dataframe with the columns ar, duration, & count
  #' int: integer measuring the time interval of the data (hours)
  
  counter <- 0
  for (ii in 1:nrow(df)) {
    if (df$ar[ii]==1) {
      if (ii==1) {
        counter <- counter+1
        df$duration[ii] <- int
      } else if (df$ar[ii-1]==0) {
        counter <- counter+1
        df$duration[ii] <- int
      } else {
        df$duration[ii] <- df$duration[ii-1]+int
      }
      df$count[ii] <- counter
    }
  }
  return(df)
}

assign_AR_cat <- function(IVT, duration) {
  ## calculates AR intensity category
  if (length(IVT) != length(duration)) stop('IVT and duration must be the same length.')

  cat <- rep(0, length(IVT))
  ivt48 <- IVT[duration > 48]
  cat[duration > 48] <- 
    case_when(ivt48>=1000 ~ 5, ivt48>=750 ~ 4, ivt48>=500 ~ 3, TRUE ~ 2)
  ivt24 <- IVT[duration > 24 & duration <= 48]
  cat[duration > 24 & duration <= 48] <- 
    case_when(ivt24>=1250 ~ 5, ivt24>=1000 ~ 4, ivt24>=750 ~ 3, ivt24>=500 ~ 2, TRUE ~ 1)
  ivt0 <- IVT[duration <= 24]
  cat[duration <= 24] <- 
    case_when(ivt0>=1250 ~ 4, ivt0>=1000 ~ 3, ivt0>=750 ~ 2, ivt0>=500 ~ 1, TRUE ~ 0)
  return(cat)
}

lag <- function(x, agg, fun, align = 'right') {
  ## performs a running calculation across a certain number of rows
  #' x: vector of numbers 
  #' agg: integer specifying the number of rows to include in the running calculation
  #' fun: string specifying what calculation to perform (accepts sum, mean, min, max, events)
  #' align: should the calculated number be saved to the left, right, or center 
  #'   of the running interval?

  if (align == 'left') {
    offset <- agg
  } else if (align == 'center') {
    offset <- agg/2
  } else if (align == 'right') {
    offset <- 0
  } else stop('Unknown input for "align."')
  
  y <- rep(NA, length(x))
  if (fun == 'sum') {
    for (i in (agg+1):length(x)) y[i-offset] <- sum(x[(i-agg):(i-1)])
  } else if (fun == 'mean') {
    for (i in (agg+1):length(x)) y[i-offset] <- mean(x[(i-agg):(i-1)])
  } else if (fun == 'Mean') {
    for (i in (agg+1):length(x)) y[i-offset] <- Mean(x[(i-agg):(i-1)])
  } else if (fun == 'max') {
    for (i in (agg+1):length(x)) y[i-offset] <- max(x[(i-agg):(i-1)])
  } else if (fun == 'min') { 
    for (i in (agg+1):length(x)) y[i-offset] <- min(x[(i-agg):(i-1)])
  } else if (fun == 'events') {
    for (i in (agg+1):length(x)) y[i-offset] <- sum(unique(x[(i-agg):(i-1)])>0)
  } else stop('Unknown input for "agg."')
  y
}


#### variable names ###############################################################################

dimension <- matrix(c(
  'H', 'ARchar', 'maxivt', 'Maximum IVT (kg/m/s)', T,
  'H', 'ARchar', 'duration', 'AR Duration (hrs)', T,
  'H', 'ARchar', 'cat', 'AR Category', T,
  'H', 'ARchar', 'logprecip_total', 'AR Total Precipitation (mm)', T,
  'H', 'ARchar', 'logprecip_max', 'AR Max Precipitation (mm/hr)', T,
  'H', 'antecedent', 'logprecip_lag03', '3-Day Lagged Total Precipitation (mm)', T,
  'H', 'antecedent', 'logprecip_lag07', '7-Day Lagged Total Precipitation (mm)', T,
  'H', 'antecedent', 'logprecip_lag14', '14-Day Lagged Total Precipitation (mm)', T,
  'H', 'antecedent', 'logprecip_lag30', '30-Day Lagged Total Precipitation (mm)', T,
  'H', 'antecedent', 'sm_lag03', '3-Day Lagged Average Soil Moisture (mm/m)', T,
  'H', 'antecedent', 'sm_lag07', '7-Day Lagged Average Soil Moisture (mm/m)', T, 
  'H', 'antecedent', 'sm_lag14', '14-Day Lagged Average Soil Moisture (mm/m)', T, 
  'H', 'antecedent', 'sm_lag30', '30-Day Lagged Average Soil Moisture (mm/m)', T, 
  'H', 'modes-lsm', 'PDO', 'PDO Climate Index', T, 
  'H', 'modes-lsm', 'ENSO', 'ENSO Climate Index', T, 
  'H', 'modes-lsm', 'pct_imperv', 'Impervious Surface Cover (%)', T, 
  'H', 'modes-lsm', 'pct_developed', 'Developed Land Cover (%)', T, 
  'H', 'modes-lsm', 'pct_wetlands', 'Wetlands Land Cover (%)', T, 
  'E', 'people', 'logpop', 'Total Population', T, 
  'E', 'people', 'logpopdensity', 'Population Density per Sq.Mi.', T, 
  'E', 'housing', 'loghu', 'Total Housing Units', T, 
  'H', 'modes-lsm', 'pct_floodplain', 'Land Area in Floodplain (%)', F, 
  'E', 'people', 'pop_pct_floodplain', 'Population in Floodplain (%)', F, 
  'E', 'housing', 'hu_pct_floodplain', 'Housing Units in Floodplain (%)', F, 
  'E', 'indicators', 'coastal', 'Coastal County Indicator', F,
  'E', 'indicators', 'riverine', 'Riverside Census Tract Indicator', F, 
  'V', 'social', 'cdc_theme1', 'SVI1 - Socioeconomic Status', T, 
  'V', 'social', 'cdc_theme2', 'SVI2 - Household Characteristics', T, 
  'V', 'social', 'cdc_theme3', 'SVI3 - Racial & Ethnic Minority Status', T, 
  'V', 'social', 'cdc_theme4', 'SVI4 - Housing Type & Transportation', T, 
  'V', 'social', 'cal_polburd', 'CalEnviroScreen Pollution Burden', F, 
  'V', 'social', 'cal_popchar', 'CalEnviroScreen Population Characteristics', F, 
  'V', 'social', 'pct_dac', 'Population in Disadvantaged Communities (%)', F, 
  'V', 'social', 'hhincome22', 'Median Household Income', T, 
  'V', 'social', 'pct_white_nonhisp', 'Non-Hispanic White Population (%)', T,
  'V', 'social', 'pct_working', 'Working Age (18-64) Population (%)', T,
  'V', 'infra', 'pct_over40', 'Housing Units over 40 Years Old (%)', F, 
  'V', 'infra', 'med_struct_age', 'Median Housing Unit Age', F, 
  'E', 'housing', 'pct_sfh', 'Single Family Homes (%)', F, 
  'V', 'infra', 'pct_ownocc', 'Owner-Occupied Homes (%)', F,
  'V', 'infra', 'pct_mobile', 'Mobile Homes (%)', F,
  'V', 'memory', 'disasters', 'Federal Disasters in the Past 3 Years', T,
  'V', 'memory', 'CRS', 'CRS Score', T,
  NA, NA, 'noise', 'Noise Variable', NA),
  ncol = 5, byrow = TRUE) %>% 
  as.data.frame %>% setNames(c('dim', 'concept', 'var', 'fullname', 'overtime'))


###################################################################################################
