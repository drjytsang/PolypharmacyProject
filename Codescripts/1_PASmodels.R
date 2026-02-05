#File1 - Prediction models for fitting the 'expected' medication count for calculation of the Polypharmacy Assessment Score (along with sensitivity analyses)

#Notes: dataset is named data, with covariates age, gender, 37 long-term conditions and 20 disease-disease interactions in the optimal model.

require("mpath")
require("zic")
require("dplyr")
require("pscl")
require("MASS")

#zero-inflated negative binomial (with interactions for count model only) - Optimal model
zinb2 <- zeroinfl(d90count~age+gender+diab_cnd+bro_cnd+ctd_cnd+ibd_cnd+chd_cnd+hyp_cnd+psm_cnd+sin_cnd+ms_cnd+strk_cnd+div_cnd+cld_cnd+blv_cnd+ld_cnd+hf_cnd+anobul_cnd+pros_cnd+thy_cnd+pvd_cnd+hear_cnd+pd_cnd+dem_cnd+atr_cnd+pep_cnd+copd_cnd+epiMM+astMM+ckdMM+con_ther4+mig_ther4+ibsMM+psoeczMM+anxdepMM+pncMM+sczbipMM+alc_high+ca_excl+cvd_gi+cvd_resp +cvd_mh +cvd_neur +cvd_sen +cvd_oth +resp_oth +gi_oth +mh_oth +neur_oth +sen_oth +mh_neur +mh_gi +mh_resp +mh_sen +neur_resp +neur_gi +neur_sen+gi_resp+gi_sen|age+gender+diab_cnd+bro_cnd+ctd_cnd+ibd_cnd+ chd_cnd+hyp_cnd+psm_cnd+sin_cnd+ms_cnd+strk_cnd+div_cnd+cld_cnd+blv_cnd+ld_cnd+hf_cnd+anobul_cnd+pros_cnd+thy_cnd+pvd_cnd+hear_cnd+pd_cnd+dem_cnd+atr_cnd+pep_cnd+copd_cnd+epiMM+astMM+ckdMM+con_ther4+mig_ther4+ibsMM+psoeczMM+anxdepMM+pncMM+sczbipMM+alc_high+ca_excl+cvd_gi+cvd_resp +cvd_mh +cvd_neur +cvd_sen +cvd_oth +resp_oth +gi_oth +mh_oth +neur_oth +sen_oth +mh_neur +mh_gi +mh_resp +mh_sen +neur_resp +neur_gi +neur_sen+gi_resp+gi_sen, data=data, dist="negbin")
summary(zinb2)
zb_counts <- predict(zinb2, type = "response")

#Backward selection (This did not drop any covariates)
zinb2b <- be.zeroinfl(zinb2, data=data, dist="negbin", alpha=0.01, trace=FALSE)
summary(zinb2b)

########### SENSITIVITY ANALYSES ################
#1) LASSO
lasso2 <- zipath(d90count~age+gender+diab_cnd+bro_cnd+ctd_cnd+ibd_cnd+chd_cnd+hyp_cnd+psm_cnd+sin_cnd+ms_cnd+strk_cnd+div_cnd+cld_cnd+blv_cnd+ld_cnd+hf_cnd+anobul_cnd+pros_cnd+thy_cnd+pvd_cnd+hear_cnd+pd_cnd+dem_cnd+atr_cnd+pep_cnd+copd_cnd+epiMM+astMM+ckdMM+con_ther4+mig_ther4+ibsMM+psoeczMM+anxdepMM+pncMM+sczbipMM+alc_high+ca_excl+cvd_gi+cvd_resp +cvd_mh +cvd_neur +cvd_sen +cvd_oth +resp_oth +gi_oth +mh_oth +neur_oth +sen_oth +mh_neur +mh_gi +mh_resp +mh_sen +neur_resp +neur_gi +neur_sen+gi_resp+gi_sen|age+gender+diab_cnd+bro_cnd+ctd_cnd+ibd_cnd+ chd_cnd+hyp_cnd+psm_cnd+sin_cnd+ms_cnd+strk_cnd+div_cnd+cld_cnd+blv_cnd+ld_cnd+hf_cnd+anobul_cnd+pros_cnd+thy_cnd+pvd_cnd+hear_cnd+pd_cnd+dem_cnd+atr_cnd+pep_cnd+copd_cnd+epiMM+astMM+ckdMM+con_ther4+mig_ther4+ibsMM+psoeczMM+anxdepMM+pncMM+sczbipMM+alc_high+ca_excl,data = data, family = "negbin", nlambda=100, lambda.zero.min.ratio=0.001, maxit.em=300, maxit.theta=25, theta.fixed=FALSE, trace=FALSE, penalty="enet", rescale=FALSE)
minBic <- which.min(BIC(lasso2))
coef(lasso2, minBic)
cat("theta estimate", lasso2$theta[minBic])
se(lasso2, minBic, log=FALSE)
AIC(lasso2)[minBic]
BIC(lasso2)[minBic]
logLik(lasso2)[minBic]

foldid <- split(sample(1:n), rep(1:10, length = n)) #initial sensitivity k-fold internal validation for LASSO model.
fitcv <- cv.zipath(d90count~age+gender+diab_cnd+bro_cnd+ctd_cnd+ibd_cnd+chd_cnd+hyp_cnd+psm_cnd+sin_cnd+ms_cnd+strk_cnd+div_cnd+cld_cnd+blv_cnd+ld_cnd+hf_cnd+anobul_cnd+pros_cnd+thy_cnd+pvd_cnd+hear_cnd+pd_cnd+dem_cnd+atr_cnd+pep_cnd+copd_cnd+epiMM+astMM+ckdMM+con_ther4+mig_ther4+ibsMM+psoeczMM+anxdepMM+pncMM+sczbipMM+alc_high+ca_excl+cvd_gi+cvd_resp +cvd_mh +cvd_neur +cvd_sen +cvd_oth +resp_oth +gi_oth +mh_oth +neur_oth +sen_oth +mh_neur +mh_gi +mh_resp +mh_sen +neur_resp +neur_gi +neur_sen+gi_resp+gi_sen|age+gender+diab_cnd+bro_cnd+ctd_cnd+ibd_cnd+ chd_cnd+hyp_cnd+psm_cnd+sin_cnd+ms_cnd+strk_cnd+div_cnd+cld_cnd+blv_cnd+ld_cnd+hf_cnd+anobul_cnd+pros_cnd+thy_cnd+pvd_cnd+hear_cnd+pd_cnd+dem_cnd+atr_cnd+pep_cnd+copd_cnd+epiMM+astMM+ckdMM+con_ther4+mig_ther4+ibsMM+psoeczMM+anxdepMM+pncMM+sczbipMM+alc_high+ca_excl, data = subdata, family = "negbin", nlambda=100,lambda.count=lasso2$lambda.count[1:100],lambda.zero= lasso2$lambda.zero[1:100],maxit.em=300, init.theta=13, theta.fixed=TRUE,penalty="enet", rescale=FALSE, foldid=region)
cat("cross-validated loglik", max(fitcv$cv))
coef(lasso2, which=fitcv$lambda.which)

#2) zero-inflated negative binomial (without any interactions)
zinb1 <- zeroinfl(d90count~age+gender+diab_cnd+bro_cnd+ctd_cnd+ibd_cnd+chd_cnd+hyp_cnd+psm_cnd+sin_cnd+ms_cnd+strk_cnd+div_cnd+cld_cnd+blv_cnd+ld_cnd+hf_cnd+anobul_cnd+pros_cnd+thy_cnd+pvd_cnd+hear_cnd+pd_cnd+dem_cnd+atr_cnd+pep_cnd+copd_cnd+epiMM+astMM+ckdMM+con_ther4+mig_ther4+ibsMM+psoeczMM+anxdepMM+pncMM+sczbipMM+alc_high+ca_excl|age+gender+diab_cnd+bro_cnd+ctd_cnd+ibd_cnd+ chd_cnd+hyp_cnd+psm_cnd+sin_cnd+ms_cnd+strk_cnd+div_cnd+cld_cnd+blv_cnd+ld_cnd+hf_cnd+anobul_cnd+pros_cnd+thy_cnd+pvd_cnd+hear_cnd+pd_cnd+dem_cnd+atr_cnd+pep_cnd+copd_cnd+epiMM+astMM+ckdMM+con_ther4+mig_ther4+ibsMM+psoeczMM+anxdepMM+pncMM+sczbipMM+alc_high+ca_excl, data=data, dist="negbin")
summary(zinb1)

#3) zero-inflated negative binomial (with interactions for both zero and count)
zinb3 <- zeroinfl(d90count~age+gender+diab_cnd+bro_cnd+ctd_cnd+ibd_cnd+chd_cnd+hyp_cnd+psm_cnd+sin_cnd+ms_cnd+strk_cnd+div_cnd+cld_cnd+blv_cnd+ld_cnd+hf_cnd+anobul_cnd+pros_cnd+thy_cnd+pvd_cnd+hear_cnd+pd_cnd+dem_cnd+atr_cnd+pep_cnd+copd_cnd+epiMM+astMM+ckdMM+con_ther4+mig_ther4+ibsMM+psoeczMM+anxdepMM+pncMM+sczbipMM+alc_high+ca_excl+cvd_gi+cvd_resp +cvd_mh +cvd_neur +cvd_sen +cvd_oth +resp_oth +gi_oth +mh_oth +neur_oth +sen_oth +mh_neur +mh_gi +mh_resp +mh_sen +neur_resp +neur_gi +neur_sen+gi_resp+gi_sen|age+gender+diab_cnd+bro_cnd+ctd_cnd+ibd_cnd+ chd_cnd+hyp_cnd+psm_cnd+sin_cnd+ms_cnd+strk_cnd+div_cnd+cld_cnd+blv_cnd+ld_cnd+hf_cnd+anobul_cnd+pros_cnd+thy_cnd+pvd_cnd+hear_cnd+pd_cnd+dem_cnd+atr_cnd+pep_cnd+copd_cnd+epiMM+astMM+ckdMM+con_ther4+mig_ther4+ibsMM+psoeczMM+anxdepMM+pncMM+sczbipMM+alc_high+ca_excl+cvd_gi+cvd_resp +cvd_mh +cvd_neur +cvd_sen +cvd_oth +resp_oth +gi_oth +mh_oth +neur_oth +sen_oth +mh_neur +mh_gi +mh_resp +mh_sen +neur_resp +neur_gi +neur_sen+gi_resp+gi_sen, data=foldall, dist="negbin")
summary(zinb3)
zinb3b <- be.zeroinfl(zinb3, data=data, dist="negbin", alpha=0.01, trace=FALSE)
summary(zinb3b)

#4) zero-inflated poisson (with interactions for count model only)
zip <- zeroinfl(d90count~age+gender+diab_cnd+bro_cnd+ctd_cnd+ibd_cnd+chd_cnd+hyp_cnd+psm_cnd+sin_cnd+ms_cnd+strk_cnd+div_cnd+cld_cnd+blv_cnd+ld_cnd+hf_cnd+anobul_cnd+pros_cnd+thy_cnd+pvd_cnd+hear_cnd+pd_cnd+dem_cnd+atr_cnd+pep_cnd+copd_cnd+epiMM+astMM+ckdMM+con_ther4+mig_ther4+ibsMM+psoeczMM+anxdepMM+pncMM+sczbipMM+alc_high+ca_excl+cvd_gi+cvd_resp +cvd_mh +cvd_neur +cvd_sen +cvd_oth +resp_oth +gi_oth +mh_oth +neur_oth +sen_oth +mh_neur +mh_gi +mh_resp +mh_sen +neur_resp +neur_gi +neur_sen+gi_resp+gi_sen|age+gender+diab_cnd+bro_cnd+ctd_cnd+ibd_cnd+ chd_cnd+hyp_cnd+psm_cnd+sin_cnd+ms_cnd+strk_cnd+div_cnd+cld_cnd+blv_cnd+ld_cnd+hf_cnd+anobul_cnd+pros_cnd+thy_cnd+pvd_cnd+hear_cnd+pd_cnd+dem_cnd+atr_cnd+pep_cnd+copd_cnd+epiMM+astMM+ckdMM+con_ther4+mig_ther4+ibsMM+psoeczMM+anxdepMM+pncMM+sczbipMM+alc_high+ca_excl, data=data, dist="poisson")
summary(zip)
zip2b <- be.zeroinfl(zip, data=data, dist="poisson", alpha=0.01, trace=FALSE)
summary(zip2b)

#5)poisson model (no zero-inflation)
p_mod <- glm(d90count~age+gender+diab_cnd+bro_cnd+ctd_cnd+ibd_cnd+chd_cnd+hyp_cnd+psm_cnd+sin_cnd+ms_cnd+strk_cnd+div_cnd+cld_cnd+blv_cnd+ld_cnd+hf_cnd+anobul_cnd+pros_cnd+thy_cnd+pvd_cnd+hear_cnd+pd_cnd+dem_cnd+atr_cnd+pep_cnd+copd_cnd+epiMM+astMM+ckdMM+con_ther4+mig_ther4+ibsMM+psoeczMM+anxdepMM+pncMM+sczbipMM+alc_high+ca_excl+cvd_gi+cvd_resp +cvd_mh +cvd_neur +cvd_sen +cvd_oth +resp_oth +gi_oth +mh_oth +neur_oth +sen_oth +mh_neur +mh_gi +mh_resp +mh_sen +neur_resp +neur_gi +neur_sen+gi_resp+gi_sen, family="poisson",x=TRUE,y=TRUE, data=data)
summary(p_mod)

#6)negative binomial model (no zero-inflation)
nb_mod <- glm.nb(d90count~age+gender+diab_cnd+bro_cnd+ctd_cnd+ibd_cnd+chd_cnd+hyp_cnd+psm_cnd+sin_cnd+ms_cnd+strk_cnd+div_cnd+cld_cnd+blv_cnd+ld_cnd+hf_cnd+anobul_cnd+pros_cnd+thy_cnd+pvd_cnd+hear_cnd+pd_cnd+dem_cnd+atr_cnd+pep_cnd+copd_cnd+epiMM+astMM+ckdMM+con_ther4+mig_ther4+ibsMM+psoeczMM+anxdepMM+pncMM+sczbipMM+alc_high+ca_excl+cvd_gi+cvd_resp +cvd_mh +cvd_neur +cvd_sen +cvd_oth +resp_oth +gi_oth +mh_oth +neur_oth +sen_oth +mh_neur +mh_gi +mh_resp +mh_sen +neur_resp +neur_gi +neur_sen+gi_resp+gi_sen, data=data)
summary(nb_mod)

#7)Other covariates(including Blood pressure, BMI and eGFR - imputed using MICE) - see 6 steps below 
require(mice)
df <- subset(data, select = c(d90count,age,gender,diab_cnd,bro_cnd,ctd_cnd,ibd_cnd,chd_cnd,hyp_cnd,psm_cnd,sin_cnd,ms_cnd,strk_cnd,div_cnd,cld_cnd,blv_cnd,ld_cnd,hf_cnd,anobul_cnd,pros_cnd,thy_cnd,pvd_cnd,hear_cnd,pd_cnd,dem_cnd,atr_cnd,pep_cnd,copd_cnd,epiMM,astMM,ckdMM,con_ther4,mig_ther4,ibsMM,psoeczMM,anxdepMM,pncMM,sczbipMM,alc_high,ca_excl,cvd_gi,cvd_resp,cvd_mh,cvd_neur,cvd_sen,cvd_oth,resp_oth,gi_oth,mh_oth,neur_oth,sen_oth, mh_neur,mh_gi,mh_resp,mh_sen,neur_resp,neur_gi,neur_sen,gi_resp,gi_sen,BPsys,BPdia,eGFR,BMI))

# Step 1: Initialize mice to get the method and predictorMatrix templates
ini <- mice(df, maxit=0)

# Step 2: Specify imputation methods
meth <- ini$method
meth["BPsys"] <- "norm"   # Impute BP using normal regression imputation
meth["BPdia"] <- "norm"   # Impute BP using normal regression imputation
meth["eGFR"] <- "norm"  # Impute eGFR using normal regression imputation
meth["BMI"] <- "norm"  # Impute BMI similarly

# Step 3: Define predictor matrix
pred <- ini$predictorMatrix
predictor_vars <- c("age", "gender", "diab_cnd","bro_cnd","ctd_cnd","ibd_cnd","chd_cnd","hyp_cnd","psm_cnd","sin_cnd","ms_cnd","strk_cnd","div_cnd","cld_cnd","blv_cnd","ld_cnd","hf_cnd","anobul_cnd","pros_cnd","thy_cnd","pvd_cnd","hear_cnd","pd_cnd","dem_cnd","atr_cnd","pep_cnd","copd_cnd","epiMM","astMM","ckdMM","con_ther4","mig_ther4","ibsMM","psoeczMM","anxdepMM","pncMM","sczbipMM","alc_high","ca_excl")

# Set all to 0 first to avoid unwanted predictors
pred[,] <- 0

# For each variable to impute, assign predictors
impute_vars <- c("BPsys","BPdia", "eGFR", "BMI")
for (var in impute_vars) {
  pred[var, predictor_vars] <- 1
}

# Step 4: Run multiple imputation
imp <- mice(df, m=10, method=meth, predictorMatrix=pred, seed=12345)

# Step 5: Fit zero-inflated negative binomial model on imputed datasets
fit_zinb <- with(imp, zeroinfl(d90count~age+gender+diab_cnd+bro_cnd+ctd_cnd+ibd_cnd+chd_cnd+hyp_cnd+psm_cnd+sin_cnd+ms_cnd+strk_cnd+div_cnd+cld_cnd+blv_cnd+ld_cnd+hf_cnd+anobul_cnd+pros_cnd+thy_cnd+pvd_cnd+hear_cnd+pd_cnd+dem_cnd+atr_cnd+pep_cnd+copd_cnd+epiMM+astMM+ckdMM+con_ther4+mig_ther4+ibsMM+psoeczMM+anxdepMM+pncMM+sczbipMM+alc_high+ca_excl+cvd_gi+cvd_resp+cvd_mh+cvd_neur+cvd_sen+cvd_oth +resp_oth+gi_oth+mh_oth +neur_oth +sen_oth +mh_neur +mh_gi+mh_resp+mh_sen+neur_resp+neur_gi+neur_sen+gi_resp+gi_sen+smok_hx+BPsys+BPdia+eGFR+BMI|age+gender+diab_cnd+bro_cnd+ctd_cnd+ibd_cnd+chd_cnd+hyp_cnd+psm_cnd+sin_cnd+ms_cnd+strk_cnd+div_cnd+cld_cnd+blv_cnd+ld_cnd+hf_cnd+anobul_cnd+pros_cnd+thy_cnd+pvd_cnd+hear_cnd+pd_cnd+dem_cnd+atr_cnd+pep_cnd+copd_cnd+epiMM+astMM+ckdMM+con_ther4+mig_ther4+ibsMM+psoeczMM+anxdepMM+pncMM+sczbipMM+alc_high+ca_excl+smok_hx+BPsys+BPdia+eGFR+BMI, dist="negbin"))

pred_list <- lapply(fit_zinb$analyses, function(model) {
  predict(model, type = "response")
})

# Step 6: Combine predictions across imputations
predicted_counts <- Reduce("+", pred_list) / length(pred_list)
summary(predicted_counts)
range(predicted_counts)
sd(predicted_counts)
mean(predicted_counts)

ddiff <- observed_counts - predicted_counts
mean(ddiff^2) #mean square error
mean(abs(ddiff)) #mean absolute error
