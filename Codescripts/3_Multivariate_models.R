library(survival)
require(dplyr)
#With dataset named data, with rdate=study entry date, IMD=index of multiple deprivation, CMS = Cambridge Multimorbidity Score

#Multivariable Cox models for outcomes (Adverse drug reactions, hospitalisations, death at 1 and 5 years)
library(survival)
cox_model <- coxph(Surv(ADRexit-rdate, ADRany) ~ PASgroup + age + gender + IMD + ethnic_6 + CMS + hx_ADRany + cluster(pracid), data = data)
summary(cox_model)

cox_model2 <- coxph(Surv(hospexit-rdate, hosp) ~ PASgroup + age + gender + IMD + ethnic_6 + CMS +  cluster(pracid), data = data)
summary(cox_model2)

cox_model3 <- coxph(Surv(death1yrexit-rdate, death1yr) ~ PASgroup + age + gender + IMD + ethnic_6 + CMS + cluster(pracid), data = data)
summary(cox_model3)

cox_model4 <- coxph(Surv(death5yrexit-rdate, death5yr) ~ PASgroup + age + gender + IMD + ethnic_6 + CMS + cluster(pracid), data = data)
summary(cox_model4)

#Multivariable logistic models for potentially inappropriate prescribing
library(lme4)
dfa <- subset(data, age >=65) #STOPP/START applied to >=65 years only (ss/sr)
dfb <- subset(data, age <65) #PROMPT applied to <65 years only (pr)

model <- glmer(ss ~ PASgroup + age + gender + IMD + ethnic_6 + CMS + (1 | pracid), data=dfa, family = binomial)
summary(model1)

model2 <- glmer(sr ~ PASgroup + age + gender + IMD + ethnic_6 + CMS + (1 | pracid), data=dfa, family = binomial)
summary(model2)

model3 <- glmer(pr ~ PASgroup + age + gender + IMD + ethnic_6 + CMS + (1 | pracid), data=dfb, family = binomial)
summary(model3)