#PLSDA code for later

# ── M6: PLS-DA (Partial Least Squares Discriminant Analysis) ──
# Uses same features as M4/M5 for direct comparison

if (!requireNamespace("systemsseRology", quietly = TRUE)) install.packages("systemsseRology")
library(systemsseRology)
  
  #visualization
  df_features <- data.frame(name = colnames(X_combined))
  df_features$label <- factor(df_features$name)

  #predictors
  X_sel <- X_combined
  X_sel <- scale(X_sel, scale = T, center = T)
  y_sel <- y_factor_by_id[rownames(X_sel)]

  #set visualization parameters
  opts_plot <- list(df_id = df_features,
                    loading_alpha = 1, # transparency for the loadings
                    score_alpha = 1, # transparency for the scores
                    LV_ind = c(1,2), # which LVs to plot
                    size = 2.5,
                    color_features = y_sel, #features should be color-coded
                    y = y_sel,
                    y_name = "HER2 status") 
  
  opts_plot2 <- list(df_features = df_features,
                    loading_alpha = 1, # transparency for the loadings
                    score_alpha = 1, # transparency for the scores
                    color_features = as.factor(y_sel), #features should be color-coded
                    LV_ind = c(1,2), # which LVs to plot
                    y_name = "HER2 status")

  # Perform a PLS-DA using the selected features and plot the scores and loadings
  # Check number of latent variables and increase to 2 if <2 
          #(for visualization purposes)
  opts_model <- list(n_LV = 2)
  
  #train PLS-DA model
  model <- train_ropls(X_sel, as.factor(y_sel), 
                       options = opts_model)
  ropls::getSummaryDF(model)
  plt_scores <- visualize_ropls_scores(model, as.factor(y_sel), 
                                       options = opts_plot)
  print(plt_scores)

  
  # set additional options required to color code enrichment in the bar plot of the loadings
  opts_plot2$X <- X_sel
  opts_plot2$y <- as.factor(df$cytokine)
  opts_plot2$LV_ind <- 1
  opts_plot2$mark_enrichment <- TRUE
  plt_loadings_bar1 <- visualize_ropls_loadings_bar(model, options = opts_plot2)
  print(plt_loadings_bar1)
  opts_plot2$LV_ind <- 2
  plt_loadings_bar2 <- visualize_ropls_loadings_bar(model, options = opts_plot2)
  print(plt_loadings_bar2)
  
  pdf(paste0("14d_plsda_",cytokine,"_LV1.pdf"), height = 4, width = 4)
  print(plt_loadings_bar1)
  dev.off()
  