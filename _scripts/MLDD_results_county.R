
###################################################################################################
## MLDD_results_county.R
## Created 4/24/23 by Corinne Bowers
###################################################################################################

#### user-defined variables #######################################################################

## set working directory
# setwd('/scratch/users/cbowers/6-MLDD/')
setwd('D:/6-MLDD/')

## set random seed
set.seed(2023)

## set number of bootstrapped samples
boot <- 100

## set file save name
filename <- 'seed2023'


#### setup ########################################################################################
cat('setting up...\n')

## load functions & packages
source('_data/setup.R')
source('_data/create_df_functions.R')
source('_data/pca_functions.R')

suppressPackageStartupMessages({
  require(DMwR) #SMOTE
  require(car) #vif
  require(factoextra) #PCA tools
  require(ggforce) #geom_circle
  require(ggrepel) #geom_text_repel
  require(BAMMtools) #getJenksBreaks
  require(randomForest)
  require(import)
  require(e1071)
  require(pdp)
})

## set up parallel backend
cores <- 5 #as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))-1
cat('  number of cores:\n')
print(cores)

## load catalog dataframe (catalog.df)
load('_data/_checkpoints/county_catalog_0426.Rdata')

## load tract geometries
california <- st_read('_data/calcounties/CA_Counties_TIGER2016.shp', quiet = TRUE)
polys <- california %>% 
  arrange(NAME) %>%
  transmute(id = 1:nrow(.), name = NAME, geometry)


#### initalize dataframe ##########################################################################
cat('initializing dataframe...\n')

## create predictor & outcome variables
df <- catalog.df %>% 
  mutate(y = factor(ifelse(claims_num > 0, 'yes', 'no'))) %>% 
  select(-count, -start, -end, -claims_num, -claims_value)

## impute missing precipitation & sm lags by linear regression
f07 <- lm(precip_lag07 ~ precip_lag03 + sm_lag03, data = df)
df$precip_lag07[is.na(df$precip_lag07)] <- 
  predict(f07, df %>% filter(is.na(precip_lag07))) %>% 
  cbind(df$precip_lag03[is.na(df$precip_lag07)]) %>% 
  apply(1,max)
f07 <- lm(sm_lag07 ~ precip_lag03 + sm_lag03, data = df)
df$sm_lag07[is.na(df$sm_lag07)] <- predict(f07, df %>% filter(is.na(sm_lag07)))

f14 <- lm(precip_lag14 ~ precip_lag07 + sm_lag07 + precip_lag03 + sm_lag03, data = df)
df$precip_lag14[is.na(df$precip_lag14)] <- 
  predict(f14, df %>% filter(is.na(precip_lag14))) %>% 
  cbind(df$precip_lag07[is.na(df$precip_lag14)]) %>% 
  apply(1,max)
f14 <- lm(sm_lag14 ~ precip_lag07 + sm_lag07 + precip_lag03 + sm_lag03, data = df)
df$sm_lag14[is.na(df$sm_lag14)] <- predict(f14, df %>% filter(is.na(sm_lag14)))

f30 <- lm(precip_lag30 ~ precip_lag14 + sm_lag14 + precip_lag07 + sm_lag07 + precip_lag03, data = df)
df$precip_lag30[is.na(df$precip_lag30)] <- 
  predict(f30, df %>% filter(is.na(precip_lag30))) %>% 
  cbind(df$precip_lag14[is.na(df$precip_lag30)]) %>% 
  apply(1,max)
f30 <- lm(sm_lag30 ~ precip_lag14 + sm_lag14 + precip_lag07 + sm_lag07 + sm_lag03, data = df)
df$sm_lag30[is.na(df$sm_lag30)] <- predict(f30, df %>% filter(is.na(sm_lag30)))

## split train & test
df <- df %>% mutate(cat = factor(cat, levels = 5:1), CRS = factor(CRS, levels = 10:2), id = factor(id))
train <- map(
  .x = unique(df$id),
  .f = ~sample(x = which(df$id==.x), size = 0.8*sum(df$id==.x))) %>% 
  reduce(c)
df.train <- df[train,]
df.test <- df[-train,]


#### SMOTE ########################################################################################
cat('performing data oversampling...\n')

# ## load policies
# load('_data/policies_save.Rdata')
# 
# ## clean up & aggregate policies by county 
# policies.poly <- policies %>% 
#   transmute(
#     county = str_sub(countyCode, start = 3, end = 5),
#     policy_start = ymd_hms(policyEffectiveDate), 
#     policy_end = ymd_hms(policyTerminationDate)) %>% 
#   filter(complete.cases(.)) %>%
#   filter(ymd('2021-10-01') %within% interval(policy_start, policy_end)) %>% 
#   left_join(california %>% st_drop_geometry %>% select(COUNTYFP,NAME), by = c('county' = 'COUNTYFP')) %>%
#   left_join(polys %>% st_drop_geometry %>% select(id,name), by = c('NAME' = 'name')) %>%
#   filter(!is.na(id))
# policies.count <- policies.poly %>% count(id)
# 
# ## load households
# households <- 
#   read.csv('_data/calcounties/households.csv', header = FALSE) %>%
#   slice(1:58) %>% 
#   mutate(V1 = V1 %>% gsub(' County\\, California', '', .)) %>% 
#   mutate(V2 = gsub('\\,', '', V2) %>% toNumber) %>% 
#   rename(name = V1, households = V2)
# 
# ## calculate insurance takeup rate by tract
# takeup <- 
#   households %>% 
#   left_join(polys %>% st_drop_geometry, by = 'name') %>% 
#   left_join(policies.count %>% rename(policies = n), by = 'id') %>% 
#   select(name, id, policies, households) %>%
#   mutate(takeup = setNA(policies/households,0))
# 
# takeup <- takeup %>% 
#   arrange(takeup) %>% 
#   # group_by(bucket = cut(takeup, 5, labels = paste0('q',1:5))) %>% 
#   group_by(
#     bucket = cut(takeup, getJenksBreaks(takeup, 6), include.lowest = TRUE, labels = paste0('q',1:5))) %>% 
#   mutate(newtakeup = Sum(policies)/sum(households)) %>% 
#   ungroup
# 
# ## checkpoint
# save(takeup, file = '_data/_checkpoints/county_takeup_0424.Rdata')
load('_data/_checkpoints/county_takeup_0424.Rdata')

## generate new SMOTE records
df.smote <- 
  foreach (q = paste0('q',1:5), .combine = 'rbind') %do% {
    data <- df.train %>% 
      filter(id %in% takeup$id[takeup$bucket==q]) %>% 
      select(-id) %>% mutate(bucket = factor(q)) %>% as.data.frame
    SMOTE(
      y ~ ., data = data, 
      perc.over = 25/takeup$newtakeup[takeup$bucket==q][1],  # minority oversampling percentage
      perc.under = 200,       # majority undersampling percentage
      k = ncol(df.train)-2)   # KNN used to generate new examples of minority class
  }

## assign county IDs to new records
vars <- names(df.train)[names(df.train) != 'y' & names(df.train) != 'id']
df.smote <- df.smote %>% 
  filter(complete.cases(.)) %>% 
  left_join(df.train %>% select(-y), by = vars) %>% 
  group_by(bucket) %>% 
  mutate(id = ifelse(is.na(id), sample(id[!is.na(id)], sum.na(id), replace = TRUE), id)) %>% 
  ungroup %>% select(-bucket)

obs.id <- data.frame(i = 1:nrow(df.smote), id = df.smote$id)
df.smote <- df.smote %>% select(-id)
cat('  number of records:\n')
nrow(df.smote)


#### variable selection ###########################################################################
cat('performing variable selection...\n')

## define dataframe
df.pca <- df.smote %>% mutate(cat = toNumber(cat), CRS = toNumber(CRS))

## iteratively remove variables based on PCA and VIF

##
vif(glm(y ~ ., data = df.pca, family = 'binomial')) %>% round(2) %>% sort(decreasing = TRUE)
plot_contribution(df.pca, plot.dim = 10)

pca.var.list <- c('hu', 'pop')
# plot_components(df.pca, pca.var.list)
predictive_power(df.pca, pca.var.list)
df.pca <- df.pca %>% select(-pop)

##
vif(glm(y ~ ., data = df.pca, family = 'binomial')) %>% round(2) %>% sort(decreasing = TRUE)
plot_contribution(df.pca, plot.dim = 10)

pca.var.list <- c('sm_lag03', 'sm_lag07', 'sm_lag14', 'sm_lag30', 'precip_lag03', 'precip_lag07', 'precip_lag14', 'precip_lag30')
plot_components(df.pca, pca.var.list)

pca.var.list <- c('sm_lag03', 'sm_lag07', 'sm_lag14', 'sm_lag30')
predictive_power(df.pca, pca.var.list)
df.pca <- df.pca %>% select(-sm_lag07, -sm_lag14, -sm_lag30)

##
# vif(glm(y ~ ., data = df.pca, family = 'binomial')) %>% round(2) %>% sort(decreasing = TRUE)
# plot_contribution(df.pca, plot.dim = 7)
# 
# pca.var.list <- c('pre1950', 'pre1960', 'pre1970', 'pre1980', 'pre1990', 'pre2000')
# plot_components(df.pca, pca.var.list)
# predictive_power(df.pca, pca.var.list)
df.pca <- df.pca %>% select(-pre1950, -pre1960, -pre1970, -pre1980, -pre2000)

##
# vif(glm(y ~ ., data = df.pca, family = 'binomial')) %>% round(2) %>% sort(decreasing = TRUE)
# plot_contribution(df.pca, plot.dim = 7)
# 
# pca.var.list <- c('pct_floodplain', 'pop_pct_floodplain', 'hu_pct_floodplain', 'popdensity', 'developed', 'pct_sfh', 'pct_mobile')
# plot_components(df.pca, pca.var.list)
# pca.var.list <- c('pct_floodplain', 'pop_pct_floodplain', 'hu_pct_floodplain')
# predictive_power(df.pca, pca.var.list)
df.pca <- df.pca %>% select(-hu_pct_floodplain, -pct_floodplain)

##
# vif(glm(y ~ ., data = df.pca, family = 'binomial')) %>% round(2) %>% sort(decreasing = TRUE)
# plot_contribution(df.pca, plot.dim = 8)
# 
# pca.var.list <- c('popdensity', 'developed', 'pct_mobile', 'pct_sfh', 'pct_dac')
# plot_components(df.pca, pca.var.list)
# 
# pca.var.list <- c('pct_mobile', 'pct_dac')
# predictive_power(df.pca, pca.var.list)
df.pca <- df.pca %>% select(-pct_mobile)

##
# pca.var.list <- c('popdensity', 'developed')
# predictive_power(df.pca, pca.var.list)
df.pca <- df.pca %>% select(-popdensity)

## 
# vif(glm(y ~ ., data = df.pca, family = 'binomial')) %>% round(2) %>% sort(decreasing = TRUE)
# plot_contribution(df.pca, plot.dim = 8)
# 
# pca.var.list <- c('cdc_theme1', 'cdc_theme2', 'cdc_theme3', 'cdc_theme4', 'cal_polburd', 'cal_popchar', 'pct_dac')
# plot_components(df.pca, pca.var.list)
# 
# pca.var.list <- c('cdc_theme4', 'cal_popchar')
# predictive_power(df.pca, pca.var.list)
df.pca <- df.pca %>% select(-cal_popchar)

## report final VIF
vif(glm(y ~ ., data = df.pca, family = 'binomial')) %>% round(2) %>% sort(decreasing = TRUE) 


#### model fit ####################################################################################
cat('fitting random forest model...\n')

## define train control parameters
# ctrl <- trainControl(method = 'repeatedcv', number = 10, repeats = 5, verboseIter = TRUE)
ctrl <- trainControl(method = 'cv', number = 10, verboseIter = TRUE)

## fit random forest model
timer <- Sys.time()
cl <- makeCluster(cores)
registerDoSNOW(cl)
f.rf <- train(
  y ~ ., data = df.pca, 
  method = 'parRF', localImp = TRUE,
  tuneGrid = expand.grid(mtry = c(1:6,8,10)),
  trControl = ctrl,
  ntree = 1000)
stopCluster(cl)
Sys.time() - timer
save(
  df.train, df.test, obs.id, df.pca, f.rf, 
  file = paste0('_results/county_rfmodel_', filename, '.Rdata'))

## report accuracy statistics
confusionMatrix(
  data = predict(f.rf, df.test %>% mutate(cat = as.numeric(cat), CRS = as.numeric(CRS))), 
  reference = df.test$y)


#### variable importance: partial dependence plots ################################################
cat('calculating partial dependence...\n')

importance <- f.rf$finalModel$importance %>% 
  as.data.frame %>% 
  rownames_to_column(var = 'varnames') %>% 
  arrange(desc(MeanDecreaseAccuracy))
vars <- importance$varnames

pb <- txtProgressBar(min = 0, max = 20, style = 3)
timer <- Sys.time()
cl <- makeCluster(min(cores,20))
registerDoSNOW(cl)
partials.boot <- 
  foreach (
    i = 1:20,
    .packages = c('randomForest', 'pdp', 'BAMMtools', 'purrr', 'tidyverse'),
    .options.snow = opts) %dopar% {
      breaks <- df.pca %>% 
        slice(sample(1:nrow(.), 1e4)) %>% 
        pull(get(vars[i])) %>% 
        getJenksBreaks(k = 20)
      map_dfc(
        .x = 1:boot,
        .f = ~pdp::partial(
          object = f.rf$finalModel, 
          train = df.pca %>% slice(sample(1:nrow(.), nrow(.), replace = TRUE)), 
          pred.var = vars[i],
          pred.grid = data.frame(breaks) %>% setNames(vars[i]),
          which.class = 'yes',
          prob = FALSE, progress = FALSE, parallel = FALSE)[,2]) %>% 
        data.frame %>% mutate(breaks) %>% 
        setNames(c(paste0('y',1:boot), vars[i]))
    }
stopCluster(cl)
Sys.time() - timer
save(
  partials.boot, 
  file = paste0('_results/county_partials_', filename, '.Rdata'))


###################################################################################################
