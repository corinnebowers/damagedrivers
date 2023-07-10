
###################################################################################################
## create_df_functions.R: This script contains functions to support the file create_df.Rmd.
## Created by Corinne Bowers 9/28/22
## Last updated 12/9/22
###################################################################################################


###################################################################################################

add_counter <- function(ar) {
  #' Creates a vector of ordered indices to distinguish sequences of true values.
  #' @param ar vector of true/false values
  #' @return vector of numbers (where ar=true) and NA (where ar=false)
  
  ## identify changepoints when the vector goes from false to true
  change <- rep(0,length(ar))
  change[which(c(0,diff(ar))>0)] <- 1

  ## convert changepoints to event indices
  count <- cumsum(change)
  count[!ar] <- NA
  return(count)
  }


###################################################################################################

create_catalog <- function(df, name, cat = TRUE, interval = 3) {
  #' Creates a catalog of events.
  #' @param df dataframe with column ts (datetime), [name] (logical), and ivt (double)
  #' @param name name of the column to base events on 
  #' @param cat logical indicating whether to calculate AR intensity categories
  #' @param interval temporal resolution of the dataframe, in the units you want "duration" to be in
  #'   (should be hours if cat=TRUE)
  #' @return dataframe of AR events with columns ar (integer), start (datetime), end (datetime), 
  #'   ivt_max (double), duration (double), and cat (integer)
  
  catalog <- df %>% 
    mutate(count = add_counter(get(name))) %>% 
    filter(get(name)) %>% 
    group_by(count) %>% 
    summarize(
      start = min(ts), 
      end = max(ts),
      maxivt = max(ivt), 
      duration = length(ivt)*interval) 
  if (cat) {
    return(catalog %>% mutate(cat = assign_AR_cat(maxivt, duration)))
  } else return(catalog)
}


###################################################################################################

calculate_sequences <- function(data, data_hist, window = 5, startend = 0.5, magcutoff = 250) {
  #' Creates a catalog of sequences.
  #' @param df dataframe with columns date (date) & ivt (double)
  #' @param window integer parameter specifying the width of the rolling average, in days
  #' @param startend double parameter specifying the quantile to use as the start/end threshold 
  #' @param magcutoff double parameter specifying the magnitude cutoff for including sequences
  
  ## calculate start-end threshold
  threshold <- data_hist %>%
    filter(wateryear(date) %in% 1981:2010) %>%
    pull(ivt) %>% 
    quantile(startend, na.rm = TRUE)
  
  ## calculate sequences
  data <- data %>% 
    mutate(
      rolling5 = lag(ivt, window, 'mean', 'center') %>% setNA(.,0),
      seq = rolling5 > threshold,
      seq.count = add_counter(seq)) %>%
    transmute(ts = date, ivt_24hr = ivt, ivt = rolling5, seq, seq.count)
  sequences <- data %>%
    create_catalog(., name = 'seq', interval = 1, cat = FALSE) %>%
    rename(maxrolling = maxivt) %>%
    filter(maxrolling > magcutoff) %>%
    left_join(data, by = c('count' = 'seq.count')) %>% 
    group_by(count) %>% 
    summarize(
      start_date = as.Date(start[1]), 
      maxrolling = maxrolling[1], 
      duration = duration[1], 
      maxivt = max(ivt_24hr))
  return(sequences)
}


###################################################################################################

attach_impacts <- function(daily, i) {
  #' Attaches all impacts data sources to daily dataframe.
  #' @param daily daily dataframe
  #' @param i index of grid cell under consideration
  #' @return daily dataframe with additional columns
  
  ## calculate start-end daily intervals of disaster declarations
  intervals <-
    disasters %>%
    filter(fips == 6e3+grid_county[i]) %>%
    map2(.x = ymd(.$start_day), .y = ymd(.$end_day), .f = ~interval(.x, .y)) %>%
    do.call('c', .)
  
  ## attach information to df_merra
  daily <- daily %>%
    ## disaster declarations
    left_join(
      data.frame(
        date = dates_merra,
        disaster = map_dbl(
          .x = dates_merra,
          .f = ~any(.x %within% intervals)) %>% as.logical),
      by = 'date') %>%
    ## public assistance
    left_join(
      pa_grid %>% filter(cell == i) %>% select(-cell),
      by = 'date') %>%
    mutate(pa_amt = case_when(date >= min(pa$start_day) ~ setNA(pa_amt,0))) %>%
    ## NCEI
    left_join(
      ncei_grid %>% filter(cell == i) %>% select(-cell),
      by = 'date') %>%
    mutate(
      ncei_damage = case_when(year(date) >= 1996 ~ setNA(ncei_damage,0)),
      ncei_event = case_when(year(date) >= 1996 ~ setNA(ncei_event, FALSE))) %>% 
    ## claims
    left_join(
      claims_grid %>% filter(cell == i) %>% select(-cell),
      by = 'date') %>%
    mutate(
      claims_num = setNA(claims_num,0),
      claims_value = setNA(claims_value,0),
      claims_coverage = setNA(claims_coverage,0),
      claims_dr = setNA(claims_dr,0)) %>%
    ## policies
    mutate(wy = wateryear(date)) %>%
    left_join(
      policies_grid %>% filter(cell == i) %>% select(-cell),
      by = c('wy' = 'coverage_wy')) %>%
    mutate(
      policies_num =
        case_when(wy >= min(policies$coverage_wy) ~ setNA(policies_num,0)),
      policies_value_adj =
        case_when(wy >= min(policies$coverage_wy) ~ setNA(policies_value_adj,0))) %>%
    rename(policies_value = policies_value_adj) %>% 
    select(-wy) %>%
    ## population
    mutate(pop = case_when(
      year(date) <= 2000 ~ pop2000_raster[i],
      year(date) <= 2010 ~ pop2010_raster[i],
      year(date) <= 2015 ~ pop2015_raster[i],
      TRUE ~ pop2020_raster[i])) %>%
    ## WWA
    left_join(
      wwa_grid %>% filter(cell == i) %>% select(-cell),
      by = 'date') %>%
    mutate(
      wwa_phenom = case_when(date >= min(wwa_grid$date) ~ setNA(wwa_phenom, 'X')),
      wwa_vtec = case_when(date >= min(wwa_grid$date) ~ setNA(wwa_vtec, 'none'))) %>%
    arrange(date)
  
  ## return dataframe
  return(daily)
}


#### functions: inflation #########################################################################

download_bls <- function() {
  ## downloads CPI data from Bureau of Labor Statistics
  read.table(
    'https://download.bls.gov/pub/time.series/cu/cu.data.2.Summaries',
    sep = '\t', header = TRUE) %>%
    filter(grepl('CUUR0000SA0', series_id)) %>% 
    separate(period, into = c('period', 'month'), sep = 1, remove = TRUE) %>% 
    mutate(month = toNumber(month)) %>%
    filter(period == 'M' & month %in% 1:12) %>% 
    transmute(year, month, cpi = value)
}

calculate_inflation <- function(yr, mo = 0, ref_yr, ref_mo = 0, bls = download_bls()) {
  ## finds inflation rate relative to some reference period 
  #' yr: year of current data
  #' mo: month of current data (0 takes annual average)
  #' ref_yr: year that the data should be converted to 
  #' ref_mo: month that the data should be converted to (0 takes annual average)
  
  if (ref_mo == 0) {
    cpi_ref <- bls %>% 
      group_by(year) %>% 
      summarize(cpi_avg = mean(cpi)) %>%
      filter(year == ref_yr) %>% 
      pull(cpi_avg)
  } else {
    cpi_ref <- bls %>% 
      filter(year == ref_yr & month == ref_mo) %>% 
      pull(cpi)
  }
  
  ## find CPI for period of interest
  if (mo == 0) {
    cpi <- bls %>% 
      group_by(year) %>% 
      summarize(cpi_avg = mean(cpi)) %>%
      filter(year == yr) %>% 
      pull(cpi_avg)
  } else {
    cpi <- bls %>% 
      filter(year == yr & month == mo) %>% 
      pull(cpi)
  }
  
  ## find value of a dollar in period of interest r.t. reference period
  cpi_ref/cpi
}


#### functions: gridded population ################################################################

gridpop <- function(pop, grid, fact = 10, error = TRUE) {
  ## produces gridded population estimates from polygons
  #' pop: sf polygon dataframe with the column pop
  #' grid: raster grid to aggregate to
  #' fact: integer determining how precise calculations should be 
  #' error: logical determining whether to report conversion error
  
  pop <- pop %>% 
    mutate(area = units::drop_units(st_area(.))) %>% #m^2
    mutate(density = pop/area)
  grid <- grid %>% 
    setValues(1) %>% 
    disaggregate(fact)
  dens <- rasterize(pop, grid, field = 'density', fun = 'mean', background = 0)
  cell <- cellSize(rast(grid), unit = 'm') %>% raster
  
  popraster <- dens*cell
  if (error) {
    pct_error <- abs(Sum(pop$pop) - Sum(popraster[])) / Sum(pop$pop)
    cat('\n')
    print(paste0('Percent Error: ', percent(pct_error, accuracy = 0.01)))
  }
  popraster <- popraster * Sum(pop$pop) / Sum(popraster[])
  popraster <- aggregate(popraster, fact) * fact^2
  
  return(popraster)
}

