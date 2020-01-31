library(jmcm)            # for cattle data
library(bamlss)
library(bamlssMVN)

# load and reshape cattle data
data("cattle", package = "jmcm")
df <- reshape(cattle, direction = "wide", timevar = "day", v.names = "weight", sep = "_")
rownames(df) <- df$id
df$id <- NULL

# Formulas with ridge penalisation 
f_mus <- paste0(paste0(names(df)[2:12], " ~ 0 + s(group, bs = 're')"))
f_lamdiags <- paste0(paste0("lamdiag", seq_len(11)), " ~ 0 + s(group, bs = 're')")
f_lambdas <- paste0(combn(seq_len(11), 2, 
			  function(x) paste0("lambda", x[1], x[2])), 
		    " ~ 0 + s(group, bs = 're')")


# Formulas without penalisation
f_mus <- paste0(paste0(names(df)[2:12], " ~ 0 + group"))
f_lamdiags <- paste0(paste0("lamdiag", seq_len(11)), " ~ 0 + group")
f_lambdas <- paste0(combn(seq_len(11), 2, 
			  function(x) paste0("lambda", x[1], x[2])), 
		    " ~ 0 + group")

# Join formulas
f <- lapply(c(f_mus, f_lamdiags, f_lambdas), FUN = as.formula)

# Estimate model
b <- bamlss(f, family = mvn_chol(k = 11), data = df, 
	    sampler = TRUE, criterion = "BIC", optimizer = bfit, 
	    burnin = 1000, thin = 5, n.iter = 2000)

b <- bamlss(f, family = mvn_chol(k = 11), data = df, 
	    sampler = FALSE, criterion = "BIC", optimizer = bfit)


# Save model
# saveRDS(b, file = "cattle_model_mcmc_burn1000_thin5.RDS")

p <- predict(b, type = "parameter")


