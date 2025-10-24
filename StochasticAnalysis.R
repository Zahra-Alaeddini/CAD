library(dplyr)
library(ggplot2)

setwd("C:/Users/ZAHRA/Desktop/Microbiome")

otu_data <- read.csv("main_OTU.csv", row.names = 1)  
metadata <- read.csv("metadata.csv")  

rownames(otu_data) <- tolower(trimws(rownames(otu_data))) 
metadata$SampleID <- tolower(trimws(metadata$SampleID)) 

threshold <- 0.8
zero_proportion <- colSums(otu_data == 0) / nrow(otu_data)
otu_data <- otu_data[, zero_proportion <= threshold]
cat("Columns removed due to >80% zeros:", names(zero_proportion[zero_proportion > threshold]), "\n")

common_samples <- intersect(rownames(otu_data), metadata$SampleID)
cat("Common samples between OTU table and metadata:", length(common_samples), "\n")

otu_dtm <- otu_data[common_samples, , drop = FALSE]
metadata <- metadata[metadata$SampleID %in% common_samples, ]

control_samples <- metadata$SampleID[metadata$health_status == "Control"]
cad_samples <- metadata$SampleID[metadata$health_status == "CAD"]

taxa_abundances_control <- colMeans(otu_dtm[control_samples, , drop = FALSE], na.rm = TRUE)
taxa_abundances_cad <- colMeans(otu_dtm[cad_samples, , drop = FALSE], na.rm = TRUE)

if (length(taxa_abundances_control) == 0) {
  cat("Warning: No data available for control group.\n")
}
if (length(taxa_abundances_cad) == 0) {
  cat("Warning: No data available for CAD group.\n")
}

if (length(taxa_abundances_control) > 0 && length(taxa_abundances_cad) > 0) {
  
  euler_maruyama <- function(dt, t_steps, initial_conditions, r_base, K_base, sigma) {
    num_taxa <- length(initial_conditions)
    abundances <- matrix(0, nrow = t_steps, ncol = num_taxa)
    abundances[1, ] <- initial_conditions
    
    for (t in 2:t_steps) {
      for (i in 1:num_taxa) {
        abundance <- abundances[t - 1, i]
        drift <- r_base * abundance * (1 - abundance / K_base)  # Logistic growth term
        diffusion <- sigma * abundance  # Stochastic perturbation term
        noise <- rnorm(1, mean = 0, sd = sqrt(dt)) * diffusion
        abundances[t, i] <- max(0, abundances[t - 1, i] + drift * dt + noise)  # Ensure non-negative abundances
      }
    }
    
    return(abundances)
  }
  
  t_steps <- 1000
  dt <- 0.1
  r_base <- 0.01
  K_base <- 1000
  sigma <- 0.001
  
  control_simulation <- euler_maruyama(dt, t_steps, taxa_abundances_control, r_base, K_base, sigma)
  cad_simulation <- euler_maruyama(dt, t_steps, taxa_abundances_cad, r_base, K_base, sigma)
  
  time_steps <- seq(0, t_steps * dt, length.out = t_steps)
  control_df <- data.frame(Time = time_steps, Abundance = rowSums(control_simulation), Group = "Control")
  cad_df <- data.frame(Time = time_steps, Abundance = rowSums(cad_simulation), Group = "CAD")
  plot_data <- rbind(control_df, cad_df)
  
  ggplot(plot_data, aes(x = Time, y = Abundance, color = Group)) +
    geom_line(alpha = 0.8, size = 1) +
    theme_minimal() +
    labs(
         x = "Time Steps",
         y = "Taxa Abundance") +
    scale_color_manual(values = c("Control" = "#1f78b4", "CAD" = "#33a02c")) +
    theme(legend.title = element_blank())
  
} else {
  cat("No valid data available for control or CAD group. Simulation cannot proceed.\n")
}


