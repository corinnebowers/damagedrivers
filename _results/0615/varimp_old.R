
###################################################################################################
## MLDD_results_county.R
## Created 4/24/23 by Corinne Bowers
###################################################################################################

#### user-defined variables #######################################################################

## set working directory
setwd('/scratch/users/cbowers/6-MLDD/')
# setwd('D:/6-MLDD')

## set random seed
set.seed(2023)

## define which analysis to run from array number
id <- commandArgs(trailingOnly=TRUE)[1]

## define which geometry to consider
geom <- commandArgs(trailingOnly=TRUE)[2]
print(geom)

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
  require(rfPermute)
})

## set up parallel backend
cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))-1
cat('  number of cores:\n')
print(cores)

## load fitted random forest model
load(paste0('_results/', geom, '_rfmodel_', id, '.Rdata'))


#### refit random forest ##########################################################################
# cat('refitting best-fit model with permutations...\n')

# timer <- Sys.time()
# f.bestfit <- randomForest(
#   y ~ ., data = df.pca, 
#   replace = FALSE, importance = TRUE,
#   ntree = 2500, mtry = f.rf$bestTune$mtry)
# f.permute <- rfPermute(
#   y ~ ., data = df.pca, 
#   replace = FALSE, proximity = TRUE,
#   ntree = 2500, mtry = f.rf$bestTune$mtry, 
#   nrep = 100, num.cores = cores)
# Sys.time() - timer

## save results
# save(
#   df.train, df.test, obs.id, df.pca, f.rf,
#   # f.bestfit, f.permute,
#   file = paste0('_results/county_rfmodel_', filename, '.Rdata'))


#### partial dependence plots #####################################################################
cat('calculating PDP variable impact metrics...\n')

importance <- f.rf$finalModel$importance %>%
  as.data.frame %>%
  rownames_to_column(var = 'varnames') %>%
  arrange(desc(MeanDecreaseAccuracy))
vars <- importance$varnames

## partial dependence plots (PDPs)
cat(' ...partial dependence plots\n')  # 3 minutes
pb <- txtProgressBar(min = 0, max = 15, style = 3)
timer <- Sys.time()
cl <- makeCluster(min(cores,15))
registerDoSNOW(cl)
partials <-
  foreach (
    i = 1:15,
    .packages = c('randomForest', 'pdp', 'BAMMtools', 'tidyverse'),
    .options.snow = opts) %dopar% {
      breaks <- df.pca %>%
        slice(sample(1:nrow(.), 1e4)) %>%
        pull(get(vars[i])) %>%
        getJenksBreaks(k = 20)
      pdp::partial(
        object = f.rf$finalModel,
        train = df.pca,
        pred.var = vars[i],
        pred.grid = data.frame(breaks) %>% setNames(vars[i]),
        which.class = 'yes',
        prob = TRUE, progress = FALSE, parallel = FALSE)
    }
stopCluster(cl)
cat('\n')
Sys.time() - timer

## individual conditional expectations (ICE)
cat(' ...individual conditional expectations\n')  # 3 minutes
pb <- txtProgressBar(min = 0, max = 15, style = 3)
timer <- Sys.time()
cl <- makeCluster(min(cores,15))
registerDoSNOW(cl)
partials.ice <-
  foreach (
    i = 1:15,
    .packages = c('randomForest', 'pdp', 'BAMMtools', 'tidyverse'),
    .options.snow = opts) %dopar% {
      breaks <- df.pca %>%
        slice(sample(1:nrow(.), 1e4)) %>%
        pull(get(vars[i])) %>%
        getJenksBreaks(k = 20)
      pdp::partial(
        object = f.rf$finalModel,
        train = df.pca,
        pred.var = vars[i],
        pred.grid = data.frame(breaks) %>% setNames(vars[i]),
        ice = TRUE,
        which.class = 'yes',
        prob = TRUE, progress = FALSE, parallel = FALSE)
    }
stopCluster(cl)
cat('\n')
Sys.time() - timer

## save results
save(
  partials, partials.ice,
  file = paste0('_results/', geom, '_pdp_', id, '.Rdata'))


#### accumulated local effects ####################################################################
# cat('calculating IML importance metrics...\n')
# 
# ## create prediction object
# predictor <- Predictor$new(
#   f.rf$finalModel,
#   data = df.pca,
#   y = df.pca$y=='yes',
#   class = 'yes',
#   type = 'prob')
# 
# timer <- Sys.time()
# cat(' ...permutation feature importance\n')  # 45 minutes
# iml.importance <- FeatureImp$new(predictor, loss = "logLoss")
# Sys.time() - timer
# 
# timer <- Sys.time()
# cat(' ...accumulated local effects\n')  # 5-10 minutes
# iml.ale <- FeatureEffects$new(predictor, method = 'ale')
# iml.ale <- iml.ale$results %>% reduce(rbind)
# Sys.time() - timer
# 
# ## save results
# save(
#   predictor, iml.ale, #iml.importance,  
#   file = paste0('_results/', geom, '_iml_', id, '.Rdata'))


#### SHAP #########################################################################################
# cat('calculating Shapley values...\n')
# 
# pb <- txtProgressBar(min = 0, max = nrow(df.pca), style = 3)
# timer <- Sys.time()
# cl <- makeCluster(cores)
# registerDoSNOW(cl)
# shap <-
#   foreach (
#     i = 1:nrow(df.pca),
#     .packages = c('iml', 'randomForest', 'tidyverse'),
#     .combine = 'rbind',
#     .options.snow = opts) %dopar% {
#     Shapley$new(predictor, x.interest = df.pca[i,])$results %>%
#         data.frame %>% column_to_rownames('feature') %>% select(phi) %>% t
#   }
# stopCluster(cl)
# cat('\n')
# Sys.time() - timer
# save(shap, file = paste0('_results/', geom, '_shap_', id, '.Rdata'))


###################################################################################################
