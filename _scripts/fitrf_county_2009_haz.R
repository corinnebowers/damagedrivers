
###################################################################################################
## fitrf_county_2009_haz.R
## Created 7/5/23 by Corinne Bowers
###################################################################################################

#### user-defined variables #######################################################################

## set working directory
setwd('/scratch/users/cbowers/6-MLDD/')
# setwd('D:/04-MLDD')

## set random seed
set.seed(2023)

## define which analysis to run from array number
id <- commandArgs(trailingOnly=TRUE)[1]


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
  require(iml)
})

## set up parallel backend
cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))-1
cat('  number of cores:\n')
print(cores)

## load catalog dataframe (catalog.df)
load('_data/_checkpoints/county2009_catalog_0705.Rdata')

## load tract geometries
california <- st_read('_data/calcounties/CA_Counties_TIGER2016.shp', quiet = TRUE)
polys <- california %>% 
  arrange(NAME) %>%
  transmute(id = 1:nrow(.), name = NAME, geometry)

## load policies
load('_data/policies_save.Rdata')


#### initalize dataframe ##########################################################################
cat('initializing dataframe...\n')

## clean up & aggregate policies by county
policies.poly <- policies %>%
  transmute(
    floodZone,
    county = str_sub(countyCode, start = 3, end = 5),
    tract = str_sub(censusTract, start = 6),
    policy_start = ymd_hms(policyEffectiveDate),
    policy_end = ymd_hms(policyTerminationDate)) %>%
  filter(complete.cases(.)) %>%
  filter(ymd('2021-10-01') %within% interval(policy_start, policy_end)) %>%
  left_join(
    california %>% st_drop_geometry %>% select(COUNTYFP,NAME),
    by = c('county' = 'COUNTYFP')) %>%
  left_join(polys %>% st_drop_geometry %>% select(id,name), by = c('NAME' = 'name')) %>%
  filter(!is.na(id))
policies.median <- policies.poly %>% count(id) %>% pull(n) %>% median

## create predictor & outcome variables
df <- catalog.df %>% 
  left_join(policies.poly %>% count(id, name = 'policies'), by = 'id') %>% 
  mutate(policies = setNA(policies,0)) %>% 
  mutate(y = factor(ifelse(claims_num > policies/policies.median, 'yes', 'no'))) %>% 
  select(-count, -start, -end, -claims_num, -claims_value, -policies)

## keep only hazard variables
df <- df[,(names(df) %in% (dimension %>% filter(dim == 'H') %>% pull(var)) | 
             names(df) %in% c('noise','id','y'))]


#### train & test #################################################################################
cat('splitting training and testing data...\n')

## split train & test
train <- map(
  .x = unique(df$id),
  .f = ~sample(x = which(df$id==.x), size = 0.8*sum(df$id==.x))) %>% 
  reduce(c)
df.train <- df[train,]
df.test <- df[-train,]


#### SMOTE ########################################################################################
cat('performing data oversampling...\n')
 
## generate new SMOTE records
data <- df.train %>% 
  mutate(cat = factor(cat, levels = 5:1)) %>% 
  as.data.frame
ratio <- sum(data$y=='yes')/nrow(data)
df.smote <- SMOTE(
  y ~ ., data = data, 
  perc.over = 200/ratio,  # minority oversampling percentage
  perc.under = 100,       # majority undersampling percentage
  k = 10)   # KNN used to generate new examples of minority class

df.smote <- df.smote %>%
  mutate(cat = toNumber(cat)) %>%
  filter(complete.cases(.))
obs.id <- data.frame(i = 1:nrow(df.smote), id = NA) 

cat('  number of records:\n')
nrow(df.smote)
table(df.smote$y)/nrow(df.smote)


#### variable selection ###########################################################################
cat('performing variable selection...\n')

## iteratively remove variables based on PCA and VIF
df.pca <- df.smote

##
# vif(glm(y ~ ., data = df.pca, family = 'binomial')) %>% round(2) %>% sort(decreasing = TRUE)
# pca.var.list <- c('sm_lag03', 'sm_lag07', 'sm_lag14', 'sm_lag30')
# plot_components(df.pca, pca.var.list)
# predictive_power(df.pca, pca.var.list)
df.pca <- df.pca %>% select(-sm_lag07, -sm_lag14, -sm_lag30)

##
# cor(df.pca %>% select(-y)) %>%
#   as.data.frame %>% rownames_to_column %>% pivot_longer(-rowname) %>%
#   filter(rowname != name) %>% arrange(desc(abs(value))) %>%
#   mutate(pair = rep(c(1,2), nrow(.)/2)) %>% filter(pair == 1) %>% select(-pair) %>%
#   filter(grepl('logprecip',rowname) & grepl('logprecip',name))
# predictive_power(df.pca, c('logprecip_total', 'logprecip_max'))
df.pca <- df.pca %>% select(-logprecip_max)
# pca.var.list <- c('logprecip_lag03', 'logprecip_lag07', 'logprecip_lag14', 'logprecip_lag30')
# plot_components(df.pca, pca.var.list)
# predictive_power(df.pca, pca.var.list)
df.pca <- df.pca %>% select(-logprecip_lag07, -logprecip_lag30)

## report final VIF
vif(glm(y ~ ., data = df.pca, family = 'binomial')) %>% round(2) %>% sort(decreasing = TRUE) 

## report final correlation coefficients
cor(df.pca %>% select(-y)) %>%
  as.data.frame %>% rownames_to_column %>% pivot_longer(-rowname) %>% 
  filter(rowname != name) %>% arrange(desc(abs(value))) %>% 
  mutate(pair = rep(c(1,2), nrow(.)/2)) %>% filter(pair == 1) %>% select(-pair) %>% 
  head


#### model fit ####################################################################################
cat('fitting random forest model...\n')  # 20-30 minutes

## define train control parameters
ctrl <- trainControl(method = 'cv', number = 10, verboseIter = TRUE)

## fit random forest model
timer <- Sys.time()
cl <- makeCluster(cores)
registerDoSNOW(cl)
f.rf <- train(
  y ~ ., data = df.pca,
  method = 'parRF',
  replace = FALSE, localImp = TRUE,
  tuneGrid = expand.grid(mtry = 1:10),
  trControl = ctrl,
  ntree = 2500)
stopCluster(cl)
Sys.time() - timer
save(
  df.train, df.test, obs.id, df.smote, df.pca, f.rf, 
  file = paste0('_results/county_rfmodel_', id, '.Rdata'))

## report accuracy statistics
confusionMatrix(data = predict(f.rf, df.test), reference = df.test$y, positive = 'yes')


###################################################################################################
