# Explicit R reference for the Wooldridge wagepan panel FE two-way-clustered-SE case.
#
# Within-transformed OLS with two-way (entity + time) clustered covariance.
# Mirrors the Greeners implementation:
#   V = sandwich(X, meat_1 + meat_2 - meat_12, X) * g/(g-1) * (n-1)/(n-k)
# where g = min(G_entity, G_time).

library(jsonlite)

case_dir <- "validation/cases/panel_fe_twoway_cluster_wagepan"
csv_path <- file.path(case_dir, "data", "wagepan.csv")

if (!file.exists(csv_path)) {
  stop("wagepan.csv is missing; run data/gen.py before the reference script")
}

df <- read.csv(csv_path)
variables <- c("lwage", "union", "married", "d81", "d82", "d83", "d84", "d85", "d86", "d87", "nr", "year")
df <- df[complete.cases(df[, variables]), variables]

y_name <- "lwage"
x_names <- c("union", "married", "d81", "d82", "d83", "d84", "d85", "d86", "d87")
entity_name <- "nr"
time_name <- "year"

within_vector <- function(values, groups) {
  values - ave(values, groups, FUN = mean)
}

y <- within_vector(df[[y_name]], df[[entity_name]])
X <- as.matrix(df[, x_names])
for (j in seq_along(x_names)) {
  X[, j] <- within_vector(X[, j], df[[entity_name]])
}

xtx_inv <- solve(crossprod(X))
beta <- as.numeric(xtx_inv %*% crossprod(X, y))
names(beta) <- x_names

residuals <- as.numeric(y - X %*% beta)
n <- nrow(X)
k <- ncol(X)

clustered_meat <- function(cluster_col) {
  clusters <- df[[cluster_col]]
  unique_clusters <- unique(clusters)
  g <- length(unique_clusters)
  meat <- matrix(0, nrow = k, ncol = k)
  for (cluster in unique_clusters) {
    idx <- clusters == cluster
    x_g <- X[idx, , drop = FALSE]
    u_g <- residuals[idx]
    meat <- meat + t(x_g) %*% u_g %*% t(u_g) %*% x_g
  }
  list(meat = meat, g = g)
}

df$inter <- interaction(df[[entity_name]], df[[time_name]])

m1 <- clustered_meat(entity_name)
m2 <- clustered_meat(time_name)
m12 <- clustered_meat("inter")

meat <- m1$meat + m2$meat - m12$meat
sandwich <- xtx_inv %*% meat %*% xtx_inv

g <- min(m1$g, m2$g)
finite_sample_correction <- (g / (g - 1)) * ((n - 1) / (n - k))
vcov <- finite_sample_correction * sandwich
# Assess admissibility without modifying raw covariance. The declared budget is
# 128*k*eps relative to matrix/spectral scale, matching the Python reference.
roundoff_multiplier <- 128
finite_covariance <- all(is.finite(vcov))
negative_variance_terms <- x_names[diag(vcov) < 0 & is.finite(diag(vcov))]
eigenvalues <- NULL
symmetry_error <- NULL
symmetry_bound <- NULL
eigenvalue_bound <- NULL
reasons <- character(0)
if (!finite_covariance) {
  reasons <- c(reasons, "covariance contains non-finite entries")
} else {
  symmetry_error <- max(abs(vcov - t(vcov)))
  symmetry_bound <- roundoff_multiplier * k * .Machine$double.eps * max(abs(vcov))
  # Only the spectral diagnostic uses the symmetric part; raw vcov is retained.
  eigenvalues <- sort(eigen(0.5 * vcov + 0.5 * t(vcov), symmetric = TRUE, only.values = TRUE)$values)
  eigenvalue_bound <- roundoff_multiplier * k * .Machine$double.eps * max(abs(eigenvalues))
  if (symmetry_error > symmetry_bound) {
    reasons <- c(reasons, "covariance is materially asymmetric")
  }
  if (length(negative_variance_terms) > 0) {
    reasons <- c(reasons, paste("negative diagonal variance for", paste(negative_variance_terms, collapse = ", ")))
  }
  if (any(eigenvalues < -eigenvalue_bound)) {
    reasons <- c(reasons, "covariance has materially negative eigenvalues")
  }
}

json_number <- function(value) {
  # Preserve finite numbers and label IEEE non-finite values in valid JSON.
  if (is.finite(value)) return(value)
  if (is.nan(value)) return("NaN")
  if (value > 0) "Infinity" else "-Infinity"
}

inference_available <- length(reasons) == 0
named_covariance <- setNames(lapply(seq_along(x_names), function(i) {
  setNames(lapply(vcov[i, ], json_number), x_names)
}), x_names)
result <- list(
  coefficients = as.list(beta),
  covariance = named_covariance,
  inference_available = inference_available,
  inference_reason = if (inference_available) NULL else paste(reasons, collapse = "; "),
  covariance_diagnostics = list(
    finite = finite_covariance,
    eigenvalues = if (is.null(eigenvalues)) NULL else as.list(eigenvalues),
    symmetry_error = symmetry_error,
    symmetry_bound = symmetry_bound,
    eigenvalue_bound = eigenvalue_bound,
    roundoff_multiplier = roundoff_multiplier,
    negative_variance_terms = as.list(negative_variance_terms)
  )
)
if (inference_available) {
  se <- sqrt(diag(vcov))
  names(se) <- x_names
  result$standard_errors <- as.list(se)
}

out_dir <- file.path(case_dir, "reference")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_json(
  result,
  file.path(out_dir, "expected.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = 16
)

cat(toJSON(result, auto_unbox = TRUE, digits = 16))
