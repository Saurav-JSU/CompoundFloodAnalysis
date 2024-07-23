# Load necessary libraries
library(dplyr)
library(fitdistrplus)
library(extRemes)
library(copula)
library(ggplot2)
library(gridExtra)

# Read the data
file_path <- "Data/sample_data.csv"
merged_data <- read.csv(file_path)

# Step 1a: Calculate the 98.5th percentile thresholds
precipitation_threshold <- quantile(merged_data$Precipitation, 0.985)
sea_level_threshold <- quantile(merged_data$sea_level, 0.985)

# Step 1b: Identify compound flood events
compound_flood_events <- merged_data %>%
  filter(Precipitation > precipitation_threshold & sea_level > sea_level_threshold)

# Step 2: Setting Thresholds
# Reconfirm the 98.5th percentile thresholds
precipitation_threshold <- quantile(merged_data$Precipitation, 0.985)
sea_level_threshold <- quantile(merged_data$sea_level, 0.985)

# Ensure data completeness for identified compound flood events
compound_flood_events <- merged_data %>%
  filter(Precipitation > precipitation_threshold & sea_level > sea_level_threshold)

# Step 3: Calculation of Empirical Return Periods
# Calculate total number of observation days
total_observation_days <- nrow(merged_data)

# Count the number of compound flood days
number_of_compound_flood_days <- nrow(compound_flood_events)

# Calculate empirical return periods
empirical_return_period <- total_observation_days / number_of_compound_flood_days

# Print the empirical return periods
print(paste("Empirical Return Period of Compound Flood Events: ", empirical_return_period, " days"))

# Step 4a: Measuring Dependence
# Calculate Kendall's tau
kendall_tau <- cor(merged_data$Precipitation, merged_data$sea_level, method = "kendall")

# Print Kendall's tau
print(paste("Kendall's Tau: ", kendall_tau))

# Step 1: Calculate the 95th percentile thresholds
precipitation_threshold_95 <- quantile(merged_data$Precipitation, 0.95, na.rm = TRUE)
sea_level_threshold_95 <- quantile(merged_data$sea_level, 0.95, na.rm = TRUE)

# Select the pairs where both values exceed the 95th percentile thresholds
selected_pairs <- merged_data %>%
  filter(Precipitation > precipitation_threshold_95 & sea_level > sea_level_threshold_95)

# Print the thresholds
print(paste("95th Percentile Threshold for Precipitation: ", precipitation_threshold_95))
print(paste("95th Percentile Threshold for Sea Level: ", sea_level_threshold_95))

# Step 2: Fit and select the best marginal distributions
fit_and_select_best <- function(data) {
  fits <- list()
  aic_values <- c()
  
  fit_lnorm <- fitdist(data, "lnorm")
  fits$lnorm <- fit_lnorm
  aic_values["lnorm"] <- fit_lnorm$aic
  
  fit_norm <- fitdist(data, "norm")
  fits$norm <- fit_norm
  aic_values["norm"] <- fit_norm$aic
  
  fit_exp <- fitdist(data, "exp")
  fits$exp <- fit_exp
  aic_values["exp"] <- fit_exp$aic
  
  fit_weibull <- fitdist(data, "weibull")
  fits$weibull <- fit_weibull
  aic_values["weibull"] <- fit_weibull$aic
  
  fit_gev <- fevd(data, type = "GEV")
  logLik_gev <- -fit_gev$results$value
  aic_gev <- -2 * logLik_gev + 2 * length(fit_gev$results$par)
  fits$gev <- fit_gev
  aic_values["gev"] <- aic_gev
  
  best_dist <- names(which.min(aic_values))
  best_fit <- fits[[best_dist]]
  
  list(best_dist = best_dist, best_fit = best_fit, aic_values = aic_values)
}

# Fit and select best marginal distributions for precipitation and sea level
best_precipitation_fit <- fit_and_select_best(selected_pairs$Precipitation)
best_sea_level_fit <- fit_and_select_best(selected_pairs$sea_level)

# Print the results
print(paste("Best-fitting distribution for Precipitation: ", best_precipitation_fit$best_dist))
print(paste("AIC values for Precipitation: ", best_precipitation_fit$aic_values))

print(paste("Best-fitting distribution for Sea Level: ", best_sea_level_fit$best_dist))
print(paste("AIC values for Sea Level: ", best_sea_level_fit$aic_values))

# Step 3: Transform the data using the best-fitting marginal distributions
get_cdf <- function(data, dist, params) {
  if (dist == "lnorm") {
    return(plnorm(data, meanlog = params["meanlog"], sdlog = params["sdlog"]))
  } else if (dist == "norm") {
    return(pnorm(data, mean = params["mean"], sd = params["sd"]))
  } else if (dist == "exp") {
    return(pexp(data, rate = params["rate"]))
  } else if (dist == "weibull") {
    return(pweibull(data, shape = params["shape"], scale = params["scale"]))
  } else if (dist == "gev") {
    return(pevd(data, loc = params[1], scale = params[2], shape = params[3], type = "GEV"))
  } else {
    stop("Unsupported distribution type")
  }
}

get_params <- function(best_fit, dist) {
  if (dist == "gev") {
    return(c(best_fit$results$par["location"], 
             best_fit$results$par["scale"], 
             best_fit$results$par["shape"]))
  } else {
    return(coef(best_fit))
  }
}

params_precip <- get_params(best_precipitation_fit$best_fit, best_precipitation_fit$best_dist)
params_sea_level <- get_params(best_sea_level_fit$best_fit, best_sea_level_fit$best_dist)

u <- get_cdf(selected_pairs$Precipitation, best_precipitation_fit$best_dist, params_precip)
v <- get_cdf(selected_pairs$sea_level, best_sea_level_fit$best_dist, params_sea_level)

# Print the first few transformed values
print(head(data.frame(u, v)))

# Step 4: Fit copula models and select the best-fitting copula
copulas <- list(
  "gaussian" = normalCopula(dim = 2),
  "clayton" = claytonCopula(dim = 2),
  "gumbel" = gumbelCopula(dim = 2),
  "frank" = frankCopula(dim = 2)
)

best_copula <- NULL
best_aic <- Inf
best_fit <- NULL

for (name in names(copulas)) {
  cop <- copulas[[name]]
  
  tryCatch({
    fit <- fitCopula(cop, data = cbind(u, v), method = "ml", optim.method = "BFGS", 
                     optim.control = list(maxit = 1000, reltol = 1e-8))
    
    aic <- AIC(fit)
    cat("Copula:", name, "AIC:", aic, "\n")
    
    if (aic < best_aic) {
      best_copula <- name
      best_aic <- aic
      best_fit <- fit
    }
  }, warning = function(w) {
    cat("Warning during copula fitting:", conditionMessage(w), "\n")
  }, error = function(e) {
    cat("Error during copula fitting:", conditionMessage(e), "\n")
  })
}

# Output the best-fitting copula and its details
cat("Best-fitting copula:", best_copula, "\n")
cat("AIC:", best_aic, "\n")
summary(best_fit)

# Step 5: Perform a goodness-of-fit test for the best-fitting copula
gof_test <- gofCopula(best_fit@copula, cbind(u, v), N = 1000, ties.method = "average")

print(gof_test)

if (!is.null(gof_test$warning)) {
  warning("Goodness-of-fit test warning: ", gof_test$warning)
}

# Step 6: Calculate copula-based joint return periods
precipitation_threshold_98.5 <- quantile(merged_data$Precipitation, 0.985, na.rm = TRUE)
sea_level_threshold_98.5 <- quantile(merged_data$sea_level, 0.985, na.rm = TRUE)

precipitation_cdf <- ecdf(merged_data$Precipitation)
sea_level_cdf <- ecdf(merged_data$sea_level)

F_R <- precipitation_cdf(precipitation_threshold_98.5)
F_S <- sea_level_cdf(sea_level_threshold_98.5)

joint_cdf <- function(u, v, copula) {
  pCopula(cbind(u, v), copula)
}

F_joint <- joint_cdf(F_R, F_S, best_fit@copula)
joint_return_period <- 1 / (1 - F_R - F_S + F_joint)

print(paste("Copula-based Joint Return Period of Compound Flood Events: ", joint_return_period, " days"))

# Step 7: Save figures and results to CSV
# Create output directories if they do not exist
if (!dir.exists("FIGURES")) {
  dir.create("FIGURES")
}
if (!dir.exists("CSV")) {
  dir.create("CSV")
}

# Save the selected pairs to a CSV file
write.csv(selected_pairs, "CSV/Selected_Pairs.csv", row.names = FALSE)

# Save the compound flood events to a CSV file
write.csv(compound_flood_events, "CSV/Compound_Flood_Events.csv", row.names = FALSE)

# Generate and save histograms for precipitation and sea level
precipitation_histogram <- ggplot(selected_pairs, aes(x = Precipitation)) +
  geom_histogram(binwidth = 5, fill = "blue", color = "black", alpha = 0.7) +
  ggtitle("Histogram of Precipitation") +
  xlab("Precipitation") +
  ylab("Frequency")

sea_level_histogram <- ggplot(selected_pairs, aes(x = sea_level)) +
  geom_histogram(binwidth = 0.1, fill = "green", color = "black", alpha = 0.7) +
  ggtitle("Histogram of Sea Level") +
  xlab("Sea Level") +
  ylab("Frequency")

# Save histograms as PNG files
ggsave("FIGURES/Precipitation_Histogram.png", precipitation_histogram, width = 7, height = 7)
ggsave("FIGURES/Sea_Level_Histogram.png", sea_level_histogram, width = 7, height = 7)

# Generate and save scatter plot of selected pairs
scatter_plot <- ggplot(selected_pairs, aes(x = Precipitation, y = sea_level)) +
  geom_point(color = "red", alpha = 0.7) +
  ggtitle("Scatter Plot of Precipitation vs Sea Level") +
  xlab("Precipitation") +
  ylab("Sea Level")

# Save scatter plot as PNG file
ggsave("FIGURES/Scatter_Plot.png", scatter_plot, width = 7, height = 7)

# Generate and save copula density plot
copula_density <- ggplot(data.frame(u, v), aes(x = u, y = v)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_gradient(low = "blue", high = "red") +
  ggtitle("Copula Density Plot") +
  xlab("Transformed Precipitation (u)") +
  ylab("Transformed Sea Level (v)")

# Save copula density plot as PNG file
ggsave("FIGURES/Copula_Density_Plot.png", copula_density, width = 5, height = 5)

# Combine and save all plots into a single image
combined_plots <- grid.arrange(precipitation_histogram, sea_level_histogram, scatter_plot, copula_density, ncol = 2)

ggsave("FIGURES/Combined_Plots.png", combined_plots, width = 7, height = 7, dpi = 600)

# Print completion message
print("Analysis complete. Results saved in 'CSV' and 'FIGURES' directories.")

