# REPRODUCE PREDICTIONS CALCULATION BY HAND

# generate mock data
mock_ids <- function() {
  return(list("123abc" = 0:39, "123def" = 40:79))
}

mock_config <- function() {
  config <- list(
    auto = TRUE,
    enabled = TRUE,
    filterSettings = list(
      regressionType = 'linear',
      regressionTypeSettings = list(
        linear = list(p.level = 0.1)
      )
    )
  )

  return(config)
}


mock_scdata <- function(with_outlier = FALSE) {
  pbmc_raw <- read.table(
    file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
    as.is = TRUE)

  if (with_outlier) {
    pbmc_raw[, 1] <- 0
    pbmc_raw[1:10, 1] <- 100
  }

  scdata <- Seurat::CreateSeuratObject(counts = pbmc_raw)
  scdata$cells_id <- 0:(ncol(scdata)-1)

  # add samples
  scdata$samples <- rep(c('123abc', '123def'), each = 40)
  return(scdata)
}

cells_id <- mock_ids()

sample_id <- "123abc"
cells_id.sample <- cells_id[[sample_id]]
scdata.sample <- subset_ids(scdata, cells_id.sample)

fit.data <- data.frame(
  log_molecules = log10(scdata.sample$nCount_RNA),
  log_genes = log10(scdata.sample$nFeature_RNA),
  row.names = scdata.sample$cells_id
)

fit.data <- fit.data[order(fit.data$log_molecules), ]

fit <- MASS::rlm(log_genes ~ log_molecules, data = fit.data)

pred.lev <- 0.95
preds <- suppressWarnings(predict(fit, interval = "prediction", level = pred.lev))


# reproduce calculations
# y = mx + b
m <-  fit$coefficients[2]
b <-  fit$coefficients[1]

x <- fit$model$log_molecules

y.hat <- (m*x)+b

# check
all(preds[,"fit"] == y.hat)

# reproduce upr and lwr calculation
#  residuals = y - preds
residuals <- fit$model$log_genes - y.hat
residuals.sq <- residuals^2
n <-  nrow(fit$model)
k <- 1 # number of explanatory variables
mse <- sum(residuals.sq)/(n-(k+1))

x.mean <- sum(fit$model$log_molecules)/length(fit$model$log_molecules)
x.sub <- fit$model$log_molecules - x.mean
x.sub.sq <- x.sub^2
sum.x.sq <- sum(x.sub.sq)

#  confidence interval (x = 1)
point.est <-  (m*1)+b
df <- n-(k+1)
t.crit <- stats::qt(1-pred.lev, df) # only step that needs a library or a pre-loaded table
s <- sqrt(mse)
tmp1 <- (1 - x.mean)^2/sum.x.sq
tmp <- (1/n) + tmp1
s.y.hat <- s*sqrt(tmp)
# point.est +- t.crit * s.y.hat
point.est + (t.crit*s.y.hat)
point.est - (t.crit*s.y.hat)

# prediction interval
# point.est +- t.crit * sqrt(mse + (s.y.hat)^2)
point.est + (t.crit * sqrt(mse + (s.y.hat)^2))
point.est - (t.crit * sqrt(mse + (s.y.hat)^2))

#  check
(y.hat[1] + (t.crit * sqrt(mse + (s.y.hat)^2))) == preds[1,"upr"]
(y.hat[1] - (t.crit * sqrt(mse + (s.y.hat)^2))) == preds[1,"lwr"]

