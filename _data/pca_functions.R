## define functions ####
plot_contribution <- function(data, plot.dim = 0) {
  ## get baseline PCA results 
  pca.results <- princomp(data %>% dplyr::select(-y), cor = TRUE)
  pca.var <- get_pca_var(pca.results)
  pca.eig <- get_eigenvalue(pca.results)
  
  print('Cumulative Variance Explained:')
  print(pca.eig$cumulative.variance.percent)
  
  ## get k-nearest neighbor clustering results 
  knn <- kmeans(pca.var$coord, centers = 5, nstart = 10, iter.max = 1e3)
  clusters <- data.frame(varnames = names(knn$cluster), cluster = unname(knn$cluster))
  
  ## plot raster contribution maps
  if (plot.dim > 0) {
    varnames <- names(data %>% dplyr::select(-y))
    temp <- pca.var$contrib %>%
      as.data.frame %>%
      cbind(varnames = row.names(.), .) %>%
      pivot_longer(cols = -varnames, names_to = 'dim', values_to = 'contrib') %>%
      full_join(data.frame(dim = row.names(pca.eig),
                           weight = pca.eig$variance.percent), by = 'dim') %>%
      mutate(dim = toNumber(gsub('Dim.', '', dim)),
             varnames = factor(varnames, levels = arrange(clusters, cluster)$varnames)) %>%
      filter(dim <= plot.dim)
    g1 <- ggplot(temp) +
      geom_tile(aes(y = varnames, x = dim, fill = contrib*weight)) +
      scale_x_continuous(breaks = 1:plot.dim) +
      scale_fill_viridis_c(name = 'Weighted \nContribution') + 
      theme(axis.title.y = element_blank())
    g2 <- ggplot(temp) +
      geom_tile(aes(y = varnames, x = dim, fill = contrib>10), color = 'grey50') +
      geom_tile(aes(y = varnames, x = dim, fill = contrib>7.5), alpha = 0.5, color = 'grey50') +
      scale_x_continuous(breaks = 1:plot.dim) + 
      scale_fill_manual(name = 'Contribution \n < 10%', values = c('grey25', 'orange')) +
      theme(axis.title.y = element_blank())
    (grid.arrange(g1, g2, ncol = 2))
  }
}
plot_components <- function(data, pca.var.list) {
  ## get PCA results for the variables in question 
  pca.results <- princomp(data[,pca.var.list], cor = TRUE)
  pca.var <- get_pca_var(pca.results)
  pca.eig <- get_eig(pca.results)
  
  print('Cumulative Variance Explained:')
  print(pca.eig$cumulative.variance.percent)
  
  ## plot component vector diagrams
  g <- ggplot(data.frame(pca.var$coord)) + 
    geom_circle(aes(x0 = 0, y0 = 0, r = 0.5), color = 'grey70') + 
    geom_circle(aes(x0 = 0, y0 = 0, r = 0.25), color = 'grey90') + 
    geom_circle(aes(x0 = 0, y0 = 0, r = 0.75), color = 'grey90') + 
    geom_circle(aes(x0 = 0, y0 = 0, r = 1)) + 
    geom_hline(yintercept = 0, color = 'grey70') + 
    geom_vline(xintercept = 0, color = 'grey70') + 
    geom_segment(aes(x = 0, xend = Dim.1, y = 0, yend = Dim.2, 
                     color = pca.var$cos2[,1]), size = 1) +
    geom_point(aes(x = Dim.1, y = Dim.2, color = pca.var$cos2[,1])) +
    # scale_color_distiller(name = '', palette = 'Spectral', direction = 1, limits = c(0,1)) + 
    scale_color_viridis_c(option = 'magma', limits = c(0,1), direction = -1) + 
    geom_text_repel(aes(x = Dim.1, y = Dim.2, label = row.names(pca.var$coord))) + 
    labs(x = 'Dim1', y = 'Dim2') +
    theme(panel.border = element_blank(), panel.grid = element_blank(),
          axis.text = element_blank(), axis.ticks = element_line(color = 'white'),
          axis.title = element_text(face = 'bold', size = 14)) + 
    scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) + 
    coord_fixed(ratio = 1)
  print(g)
}
predictive_power <- function(data, pca.var.list) {
  ## get PCA results for the variables in question 
  pca.results <- princomp(data[,pca.var.list], cor = TRUE)
  pca.var <- get_pca_var(pca.results)
  
  var.AIC <- rep(NA, length(pca.var.list))
  for (i in 1:length(pca.var.list)) {
    logit <- glm(factor(y) ~ get(pca.var.list[i]), data = data, family = 'binomial')
    if (summary(logit)$coefficients[2,4] < 0.05) { var.AIC[i] <- AIC(logit) }
  }
  comp <- glm(factor(data$y) ~ pca.results$scores[,1], family = 'binomial')
  # return(ifelse(AIC(comp) < Min(var.AIC), 'dim1', pca.var.list[which.min(var.AIC)]))
  return(pca.var.list[which.min(var.AIC)])
}
