#File 2 - Internal-External Cross validation (given our clustered dataset) for the Polypharmacy Assessment Score construction.
#Done in 2 stages - 1. Loops to create separate models for all practices (n=1495) except for one, then replacement 2. Meta-analyse performance statistics
#NB dataset is named - foldall, see other scripts for covariates details

# Get libraries needed
library(pscl)
library(dplyr)
library(boot)
library(meta) 

set.seed(12345)  # For reproducibility

#Step1 - Loop through each practice ID for separate models
for (id in 1:1495) {
  train <- foldall[foldall$prid != id, ]  # Exclude the current practice
  test <- foldall[foldall$prid == id, ]   # Use the current practice as test
  model <- zeroinfl(d90count~age+gender+diab_cnd+bro_cnd+ctd_cnd+ibd_cnd+chd_cnd+hyp_cnd+psm_cnd+sin_cnd+ms_cnd+strk_cnd+div_cnd+cld_cnd+blv_cnd+ld_cnd+hf_cnd+anobul_cnd+pros_cnd+thy_cnd+pvd_cnd+hear_cnd+pd_cnd+dem_cnd+atr_cnd+pep_cnd+copd_cnd+epiMM+astMM+ckdMM+con_ther4+mig_ther4+ibsMM+psoeczMM+anxdepMM+pncMM+sczbipMM+alc_high+ca_excl+cvd_gi+cvd_resp +cvd_mh +cvd_neur +cvd_sen +cvd_oth +resp_oth +gi_oth +mh_oth +neur_oth +sen_oth +mh_neur +mh_gi +mh_resp +mh_sen +neur_resp +neur_gi +neur_sen+gi_resp+gi_sen|age+gender+diab_cnd+bro_cnd+ctd_cnd+ibd_cnd+ chd_cnd+hyp_cnd+psm_cnd+sin_cnd+ms_cnd+strk_cnd+div_cnd+cld_cnd+blv_cnd+ld_cnd+hf_cnd+anobul_cnd+pros_cnd+thy_cnd+pvd_cnd+hear_cnd+pd_cnd+dem_cnd+atr_cnd+pep_cnd+copd_cnd+epiMM+astMM+ckdMM+con_ther4+mig_ther4+ibsMM+psoeczMM+anxdepMM+pncMM+sczbipMM+alc_high+ca_excl, data=train, dist="negbin")  
  # Predict on the test set
  preds <- predict(model, newdata = test, type = "response")
  ddiff <- test$d90count - preds
  # Store results 
  results[[as.character(id)]] <- data.frame(predicted = preds, diff = ddiff)
  results2[[as.character(id)]] <- data.frame(MSE = mean(ddiff^2), MAE = mean(abs(ddiff)), AIC = AIC(model), BIC = AIC(model, k=log(dim(test)[1])), logLik(model))
}
# Combine and save all results if needed
save(results, file = "resultscv.RData")
save(results2, file = "results2cv.RData")


####Step 2: Combine and meta-analyse performance estimations#### 

#core functions (mae = mean absolute error and mse=mean squared error)
calc_mae_boot <- function(data,indices) {
    d <- data[indices,]
    mean(abs(d$obs-d$pred))
  }

calc_mse_boot <- function(data,indices) {
  d <- data[indices,]
  mean((d$obs-d$pred)^2)
  }

#sample size per centre
  patients_per_centre <- 0
  for (id in 1:1475) {
    load(paste0("resdat", id, ".RData"))
    patients_per_centre[id] <- nrow(resdat[[as.character(id)]])
  }
  
#create empty vectors to store data
mse_est_vec <- c()
mse_se_vec <- c()
mae_est_vec <- c()
mae_se_vec <- c()

# Combine results for practices by extracting data
for (id in 1:1475) {
  load(paste0("reso", id, ".RData"))
  load(paste0("resdat", id, ".RData"))
  
  #Extract predictions 
  r<- c()
  r$obs <- resdat[[as.character(id)]][["predicted"]] + resdat[[as.character(id)]][["diff"]]
  r$pred <- resdat[[as.character(id)]][["predicted"]]
  r <- as.data.frame(testr)
  
  # bootstrap MAE mean and standard error
  mae_boot <- boot(testr,calc_mae_boot, 1000)
  reso[[as.character(id)]]$mae_est_vec <- c(mae_est_vec,mae_boot$t0)
  reso[[as.character(id)]]$mae_se_vec <- c(mae_se_vec,sd(mae_boot$t))
  
  # bootstrap MSE mean and standard error
  mse_boot <- boot(testr,calc_mse_boot, 1000)
  reso[[as.character(id)]]$mse_est_vec <- c(mse_est_vec,mse_boot$t0)
  reso[[as.character(id)]]$mse_se_vec <- c(mse_se_vec,sd(mse_boot$t))
  
  #combine results
  dftemp <- do.call(rbind, reso)
  results <- rbind(results, dftemp)
}

meta_mse <- metagen(
  TE = mse_est_vec,
  seTE = mse_se_vec,
  studlab = seq(1,10,1),
  random = TRUE,
  n.e = patients_per_centre
)

# Results
print(meta_mse)
forest(meta_mse)
print(meta_mae)
forest(meta_mae)

#check for normal distribution
hist(mse_est_vec)
qqnorm(mse_est_vec, main = "Q-Q Plot of MSE Values")
qqline(mse_est_vec, col = "red") # Add a reference line
shapiro.test(mse_est_vec)
#W = 0.9141, p-value = 0.3104
plot(mse_est_vec)
