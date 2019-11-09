#Q2 code
library(mice)

set.seed(1)

load("databp.Rdata")

#add means and standard error = sd/sqrt(n)

#a
databp <- subset(databp, select = -c(R))

databp_cc <- databp[complete.cases(databp),]

cor(databp_cc, method = "pearson")

#b
databp_mi <- complete(mice(data = databp, method = "mean", m = 1))

cor(databp_mi, method = "pearson")

#c
databp_ri <-
  complete(mice(data = databp, method = "norm.predict", m = 1))

cor(databp_ri, method = "pearson")

#d
databp_sri <-
  complete(mice(data = databp, method = "norm.nob", m = 1))
databp_sri$recovtime <-
  as.numeric(lapply(databp_sri$recovtime, max, 0))
cor(databp_sri, method = "pearson")

#e
databp_pmm <-
  complete(mice(
    data = databp,
    method = "pmm",
    m = 1,
    donors = 1
  ))

cor(databp_pmm, method = "pearson")

##########

cc_lm <-
  lm(data = databp_cc, recovtime ~ .) #standard behaviour in R is for lm to drop nas, this just makes it explicit

na_ind <- which(is.na(databp$recovtime))

#complete case analysis
databp_cc <- databp[complete.cases(databp),]


#mean imputation
databp_mi_sr <- databp
databp_mi_sr$recovtime[na_ind] <-
  mean(databp_mi_sr$recovtime, na.rm = TRUE)

#regression imputation
databp_ri_sr <- databp
databp_ri_sr$recovtime[na_ind] <-
  predict(cc_lm, databp_ri_sr[na_ind,])

#stochastic regression imputation
databp_sri_sr <- databp
databp_sri_sr$recovtime[na_ind] <-
  predict(cc_lm, databp_sri_sr[na_ind,]) + rnorm(n = length(na_ind),
                                                 mean = 0,
                                                 sd = sigma(cc_lm))
for (i in 1:length(na_ind)) {
  while (databp_sri_sr$recovtime[na_ind[i]] < 0) {
    databp_sri_sr$recovtime[na_ind[i]] <-
      predict(cc_lm, databp_sri_sr[na_ind[i],]) + rnorm(n = 1,
                                                        mean = 0,
                                                        sd = sigma(cc_lm))
  }
}


databp_sri_sr$recovtime <-
  as.numeric(lapply(databp_sri_sr$recovtime, max, 0))

#predictive mean matching

databp_pmm_sr <- databp
databp_pmm_srnl <- databp

sq_diff <- vector("list", length(na_ind))
to_select <- double(length(na_ind))
for (i in 1:length(sq_diff)) {
  sq_diff[i] <-
    list((cc_lm$fitted.values - predict(cc_lm, databp[na_ind[i],])) ^ 2)
  to_select[i] <- which.min(sq_diff[[i]])
}

for (i in 1:length(na_ind)) {
  databp_pmm_sr$recovtime[na_ind[i]] <-
    databp_cc$recovtime[which.min(sq_diff[[i]])]
}

#more efficient vectorised version of above
databp_pmm_srnl$recovtime[na_ind] <- databp_cc$recovtime[to_select]


cor(databp_cc, method = "pearson")
cor(databp_mi_sr, method = "pearson")
cor(databp_ri_sr, method = "pearson")
cor(databp_sri_sr, method = "pearson")
cor(databp_pmm_sr, method = "pearson")
mean(databp_cc$recovtime)
mean(databp_mi_sr$recovtime)
mean(databp_ri_sr$recovtime)
mean(databp_sri_sr$recovtime)
mean(databp_pmm_sr$recovtime)
mean(databp_pmm_srnl$recovtime)
sd(databp_cc$recovtime) / sqrt(sum(!is.na(databp_cc$recovtime)))
sd(databp_mi_sr$recovtime) / sqrt(sum(!is.na(databp_mi_sr$recovtime)))
sd(databp_ri_sr$recovtime) / sqrt(sum(!is.na(databp_ri_sr$recovtime)))
sd(databp_sri_sr$recovtime) / sqrt(sum(!is.na(databp_sri_sr$recovtime)))
sd(databp_pmm_sr$recovtime) / sqrt(sum(!is.na(databp_pmm_sr$recovtime)))
sd(databp_pmm_srnl$recovtime) / sqrt(sum(!is.na(databp_pmm_srnl$recovtime)))

