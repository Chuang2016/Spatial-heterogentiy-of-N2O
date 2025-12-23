gc()
rm(list = ls())
# load required packages
library(rio)
require(ggplot2)
require(dplyr)
require(lubridate)
require(reshape2)
require(randomForest)
require(rfPermute)
require(A3)
require(kernelshap)
require(shapviz)
require(caret)

##load dataset of N2O and environmental factors 
mydata.gas <- import(paste0("Figure 3 dataset for spatial heteorogentity of N2O.xlsx"), sheet = "winter wheat")
head(mydata.gas)
########################################################################################################################################
# RF-SHAP model
########################################################################################################################################


dat.all.topsoil <- mydata.gas 

colnames(dat.all.topsoil);head(dat.all.topsoil)

variable.labels <- c(
  "NH4_N_mg_kg_0_20"= expression(NH[4]^'+'*"-N"),
  "NO3_N_mg_kg_0_20"= expression(NO[3]^'-'*"-N"),
  "MBC_mg_kg_0_20" = expression(MBC), 
  "MBN_mg_kg_0_20" =expression(MBN), 
  "CN_0_20"  =expression(CNR),              
  "SWC_0_20" =expression(SWC),
  "SOM_0_20"=expression(SOM),  
  "pH_0_20"=expression(pH) , 
  "Urease_mg_g_0_20"=expression(Urease) , 
  "Catalase_ml_g_0_20"=expression(Catalase) ,
  "β_glucosidase_mg_g_0_20" =expression("β-"~glucosidase),
  "Nitrate_reductase_mg_g_0_20"=expression(Nitrate~reductase) ,  
  "Nitrite_reductase_mg_g_0_20"=expression(Nitrite~reductase),
  "Xylanase_umol_d_0_20" =expression(Xylanase), 
  "NO2_N_0_20"=expression(NO[2]^'-'*"-N"),
  "DOC_0_20_mg_kg"=expression(DOC), 
  "soil_tem"= expression(TEM[soil]),
  "soil_EC"=expression(EC),
  "soil_O2"=expression(O[2*soil])
)

# add group color
fill_colors <- c("NH4_N_mg_kg_0_20" = "#a3a500", "NO3_N_mg_kg_0_20" ="#a3a500","DOC_0_20_mg_kg"="#a3a500",
                 "MBC_mg_kg_0_20"="#00b0f6", "MBN_mg_kg_0_20"="#00b0f6", "CN_0_20"="#00b0f6" ,
                 "SWC_0_20" ="#e76bf3","SOM_0_20" ="#F8766d", "pH_0_20" ="#F8766d",  "NO2_N_0_20" ="#a3a500",
                 "Urease_mg_g_0_20"="#00bf7d", "Catalase_ml_g_0_20"="#00bf7d", "β_glucosidase_mg_g_0_20"="#00bf7d" ,
                 "Nitrate_reductase_mg_g_0_20"="#00bf7d" ,  "Nitrite_reductase_mg_g_0_20"="#00bf7d","Xylanase_umol_d_0_20"="#00bf7d",
                 "soil_tem"= "#e76bf3", "soil_EC" = "#F8766d", "soil_O2" ="#e76bf3" )



##training model
####Cross-Validation settings
set.seed(123)
train_control <- trainControl(method = "cv", number = 2, savePredictions = "final",  # Save predictions for evaluation
                              verboseIter = TRUE )   # Print progress
rf<-train(N2O_ug_N_m2_hour ~., data = dat.all.topsoil.2[
  ,
  c( 3:22)], method = 'rf', trControl = train_control,
  metric = "RMSE",tuneLength = 5
)

##define SHAP model prection function

pred_fun<-function(model, newdata) {
  as.numeric(predict(model, newdata, type = "raw"))
}

# calculate
system.time(
  ps <- kernelshap(
    rf,
    dat.all.topsoil.2[
      ,
      c( 4:22)] ,  # selected variables 
    pred_fun = pred_fun
  )
)

#Feature importance was quantified using SHAP values, calculated as the mean absolute SHAP value across all samples and normalized to percentage contributions.
shap_mat <- ps$S              
mean_abs_shap <- colMeans(abs(shap_mat))  

pct_contrib <- mean_abs_shap / sum(mean_abs_shap) * 100


#------------------------ Extract Performance Metrics ---

# 1. Cross-Validation Metrics from RF Model
best_metrics <- rf$results[rf$results$mtry == rf$bestTune$mtry, ]
cv_results <- rf$resample

cat("Cross-Validation Performance Metrics (Best Model):\n")
cat("R²: ", best_metrics$Rsquared, "\n")
cat("RMSE: ", best_metrics$RMSE, "\n")
cat("MAE: ", best_metrics$MAE, "\n")
cat("\nCross-Validation (Averaged Across Folds):\n")
cat("Mean R²: ", mean(cv_results$Rsquared), " (SD: ", sd(cv_results$Rsquared), ")\n")
cat("Mean RMSE: ", mean(cv_results$RMSE), " (SD: ", sd(cv_results$RMSE), ")\n")
cat("Mean MAE: ", mean(cv_results$MAE), " (SD: ", sd(cv_results$MAE), ")\n")

# 2. Metrics for SHAP Predictions
# Get actual values (newN2O_mean) and predicted values using pred_fun
actual <- dat.all.topsoil$newN2O_mean
predicted <- pred_fun(rf, dat.all.topsoil[, c(4:22)])
#mean_N2O_obseved <- dat.all.topsoil.2$newN2O_mean %>% mean()
# Calculate R², RMSE, and MAE
residuals <- actual - predicted
rmse <- sqrt(mean(residuals^2))/mean(actual )  
mae <- mean(abs(residuals))/mean(actual ) 
ss_total <- sum((actual - mean(actual))^2)
ss_residual <- sum(residuals^2)
r_squared <- 1 - ss_residual / ss_total

cat("\nPerformance Metrics for SHAP Predictions:\n")
cat("R²: ", r_squared, "\n")
cat("RMSE: ", rmse, "\n")
cat("MAE: ", mae, "\n")
#------------------------#------------------------


#------------------------Visualization------------
library(dplyr)
library(ggplot2)

# create a dataframe for visualization
df_shap <- tibble(
  feature = names(pct_contrib),
  pct = pct_contrib
) %>%
  arrange(pct) %>%  
  mutate(feature = factor(feature, levels = feature))  


p_pct <- ggplot(df_shap, aes(x = pct, y = feature, fill = feature)) +
  geom_col(show.legend = T) +
  scale_fill_manual(values = fill_colors) +
  scale_x_continuous(
    labels = function(x) paste0(formatC(x, format = "f", digits = 0))
  ) +
  labs(
    x = "Average SHAP contribution (%)",
    y = ""
  ) +
  theme_bw(base_family = "serif", base_size = 12) +
  scale_y_discrete(
    labels = variable.labels

  ) +
  theme_bw(base_family = "serif", base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, size = 12, family = "serif"),
    axis.text.y = element_text(size = 12, family = "serif"),
    axis.title = element_text(size = 12, family = "serif"),
    plot.title = element_text(size = 12, family = "serif"),
    panel.grid.major.y = element_blank(),
    strip.text = element_text(size = 12, family = "serif"),
    strip.text.x.top = element_blank(),
    strip.background = element_blank(),
    legend.background = element_blank(),
    legend.title = element_text(size = 12, family = "serif"),
    legend.text = element_text(size = 12, family = "serif")
  )
p_pct+ canvas(width = 16, height = 16, units = "cm", dpi = 200)


ggsave(p_pct, filename = paste0( "Figure 3a.pdf"),
       width = 24, height = 16, units = "cm")

