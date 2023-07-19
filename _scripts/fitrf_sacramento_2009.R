
###################################################################################################
## fitrf_sacramento_2009.R
## Created 4/24/23 by Corinne Bowers
## latest modification (7/5): updated PCA to include new variables
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
load('_data/_checkpoints/sac2009_catalog_0705.Rdata')

## load tract geometries
caltracts <- st_read('_data/caltracts/cb_2022_06_tract_500k.shp', quiet = TRUE)
polys <- caltracts %>% 
  filter(COUNTYFP == '067') %>% 
  arrange(TRACTCE) %>%
  transmute(id = 1:nrow(.), name = TRACTCE, geometry)


#### initalize dataframe ##########################################################################
cat('initializing dataframe...\n')

## create predictor & outcome variables
df <- catalog.df %>% 
  mutate(y = factor(ifelse(claims_num > 0, 'yes', 'no'))) %>% 
  select(-count, -start, -end, -claims_num, -claims_value)

## report results
cat('  number of records:\n')
nrow(df)
table(df$y)
(table(df$y)/nrow(df)*100) %>% round(2)


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
  mutate(
    cat = factor(cat, levels = 5:1), 
    CRS = factor(CRS, levels = 10:2), 
    riverine = factor(riverine),
    disasters = factor(disasters),
    pct_dac = factor(pct_dac)) %>% 
  as.data.frame %>% 
  select(-id)
ratio <- sum(data$y=='yes')/nrow(data)
df.smote <- SMOTE(
  y ~ ., data = data, 
  perc.over = 25/ratio,  # minority oversampling percentage
  perc.under = 100,       # majority undersampling percentage
  k = 10)   # KNN used to generate new examples of minority class

df.smote <- df.smote %>%
  mutate(
    cat = toNumber(cat),
    CRS = toNumber(CRS),
    riverine = as.numeric(riverine)-1,
    disasters = toNumber(disasters),
    pct_dac = toNumber(pct_dac)) %>%
  filter(complete.cases(.))
obs.id <- data.frame(i = 1:nrow(df.smote), id = NA) #df.smote$id)

cat('  number of records:\n')
nrow(df.smote)
table(df.smote$y)/nrow(df.smote)


#### variable selection ###########################################################################
cat('performing variable selection...\n')

## iteratively remove variables based on PCA and VIF
df.pca <- df.smote

##
# vif(glm(y ~ ., data = df.pca, family = 'binomial')) %>% round(2) %>% sort(decreasing = TRUE) %>% head
# predictive_power(df.pca, c('pct_floodplain', 'pop_pct_floodplain', 'hu_pct_floodplain'))
df.pca <- df.pca %>% select(-hu_pct_floodplain, -pct_floodplain)

##
# vif(glm(y ~ ., data = df.pca, family = 'binomial')) %>% round(2) %>% sort(decreasing = TRUE) %>% head
# pca.var.list <- c('sm_lag03', 'sm_lag07', 'sm_lag14', 'sm_lag30')
# plot_components(df.pca, pca.var.list)
# predictive_power(df.pca, pca.var.list)
df.pca <- df.pca %>% select(-sm_lag07, -sm_lag14, -sm_lag30)

##
# cor(df.pca %>% select(-y)) %>%
#   as.data.frame %>% rownames_to_column %>% pivot_longer(-rowname) %>%
#   filter(rowname != name) %>% arrange(desc(abs(value))) %>%
#   mutate(pair = rep(c(1,2), nrow(.)/2)) %>% filter(pair == 1) %>% select(-pair) %>%
#   filter(grepl('logprecip',rowname) & grepl('logprecip',name)) %>% head
# predictive_power(df.pca, c('logprecip_total', 'logprecip_max'))
df.pca <- df.pca %>% select(-logprecip_max)
# pca.var.list <- c('logprecip_lag03', 'logprecip_lag07', 'logprecip_lag14', 'logprecip_lag30')
# plot_components(df.pca, pca.var.list)
# predictive_power(df.pca, pca.var.list)
df.pca <- df.pca %>% select(-logprecip_lag07, -logprecip_lag30)

##
# vif(glm(y ~ ., data = df.pca, family = 'binomial')) %>% round(2) %>% sort(decreasing = TRUE) %>% head
# cor(df.pca %>% select(-y)) %>%
#   as.data.frame %>% rownames_to_column %>% pivot_longer(-rowname) %>%
#   filter(rowname != name) %>% arrange(desc(abs(value))) %>%
#   mutate(pair = rep(c(1,2), nrow(.)/2)) %>% filter(pair == 1) %>% select(-pair) %>%
#   filter(grepl('cdc',rowname)) %>% head
# pca.var.list <- c('cdc_theme1', 'cdc_theme3', 'pct_white_nonhisp', 'cal_popchar', 'hhincome22')
# plot_components(df.pca, pca.var.list)
# predictive_power(df.pca, c('cdc_theme1', 'cal_popchar'))
df.pca <- df.pca %>% select(-cdc_theme1, -cdc_theme3)

##
# cor(df.pca %>% select(-y)) %>%
#   as.data.frame %>% rownames_to_column %>% pivot_longer(-rowname) %>%
#   filter(rowname != name) %>% arrange(desc(abs(value))) %>%
#   mutate(pair = rep(c(1,2), nrow(.)/2)) %>% filter(pair == 1) %>% select(-pair) %>% head
# predictive_power(df.pca, c('logpop', 'loghu'))
df.pca <- df.pca %>% select(-logpop)
# predictive_power(df.pca, c('med_struct_age', 'pct_over40'))
df.pca <- df.pca %>% select(-med_struct_age)
# predictive_power(df.pca, c('pct_imperv', 'pct_developed'))
df.pca <- df.pca %>% select(-pct_developed)

##
# vif(glm(y ~ ., data = df.pca, family = 'binomial')) %>% round(2) %>% sort(decreasing = TRUE) %>% head
# pca.var.list <- c('pct_sfh', 'pct_ownocc', 'cal_popchar', 'hhincome22')
# plot_components(df.pca, pca.var.list)
# predictive_power(df.pca, c('pct_sfh', 'pct_ownocc'))
df.pca <- df.pca %>% select(-pct_ownocc)

## 
# cor(df.pca %>% select(-y)) %>%
#   as.data.frame %>% rownames_to_column %>% pivot_longer(-rowname) %>%
#   filter(rowname != name) %>% arrange(desc(abs(value))) %>%
#   mutate(pair = rep(c(1,2), nrow(.)/2)) %>% filter(pair == 1) %>% select(-pair) %>% head
df.pca <- df.pca %>% select(-cat)

## report final VIF
vif(glm(y ~ ., data = df.pca, family = 'binomial')) %>% round(2) %>% sort(decreasing = TRUE) %>% head

## report final correlation coefficients
cor(df.pca %>% select(-y)) %>%
  as.data.frame %>% rownames_to_column %>% pivot_longer(-rowname) %>% 
  filter(rowname != name) %>% arrange(desc(abs(value))) %>% 
  mutate(pair = rep(c(1,2), nrow(.)/2)) %>% filter(pair == 1) %>% select(-pair) %>% head


#### model fit ####################################################################################
cat('fitting random forest model...\n')

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
  df.train, df.test, obs.id, df.pca, f.rf, 
  file = paste0('_results/sac_rfmodel_', id, '.Rdata'))

## report accuracy statistics
confusionMatrix(data = predict(f.rf, df.test), reference = df.test$y, positive = 'yes')


###################################################################################################
