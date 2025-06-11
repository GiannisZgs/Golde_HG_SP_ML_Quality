library(jsonlite)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
#data <- readMat("/media/giannis/Extreme SSD/similarity_analysis_results_test.mat")
data <- fromJSON(txt = "/media/giannis/Extreme SSD/similarity_analysis_results.json", simplifyVector = FALSE)

metric_matrices <- c("CorrMat", "NrmseMat", "CosMat", "JsdMat", "EmdMat")

sim_formal_names <- c(
  CorrMat = "Cross-Correlation",
  NrmseMat = "NRMSE",
  CosMat = "Cosine Similarity"
)


for (matrix_name in metric_matrices) {
  if (matrix_name %in% names(data)) {
    data[[matrix_name]][["p7"]] <- NULL
  }
}

if ("SigQual" %in% names(data)) {
  # Get all metrics within SigQual
  sig_metrics <- names(data[["SigQual"]])
  
  # For each metric in SigQual, remove p7
  for (metric_name in sig_metrics) {
    if ("p7" %in% names(data[["SigQual"]][[metric_name]])) {
      data[["SigQual"]][[metric_name]][["p7"]] <- NULL
    }
  }
}

extract_sigqual_metrics <- function(data) {
  if ("SigQual" %in% names(data)) {
    # Get all metrics within SigQual
    sig_metrics <- names(data[["SigQual"]])
    
    # For each metric, extract it to the top level
    for (metric_name in sig_metrics) {
      # Move the metric to top level
      data[[metric_name]] <- data[["SigQual"]][[metric_name]]
    }
    
    # Optionally, you can remove the original SigQual structure if no longer needed
    data[["SigQual"]] <- NULL
  }
  
  return(data)
}

# Apply the extraction
data <- extract_sigqual_metrics(data)

prepare_sigqual_data <- function(data, metrics) {
  df_list <- list()
  
  for (metric_name in metrics) {
    if (metric_name %in% names(data)) {
      temp_df <- data.frame()
      
      for (p_id in names(data[[metric_name]])) {
        for (sensor in c("HG", "AgCl")) {
          if (sensor %in% names(data[[metric_name]][[p_id]])) {
            for (ch in names(data[[metric_name]][[p_id]][[sensor]])) {
              value <- data[[metric_name]][[p_id]][[sensor]][[ch]]
              
              # Skip if value is NULL or empty
              if (is.null(value) || length(value) == 0) {
                next
              }
              # Handle both vector and scalar cases
              if (length(value) > 1) {
                for (v in value) {
                  #Skip NaN
                  if (is.na(v) || is.nan(v)) {
                    next
                 }
                  temp_df <- rbind(temp_df, data.frame(
                    Metric = metric_name,
                    Participant = p_id,
                    Sensor = sensor,
                    Channel = ch,
                    Value = v
                  ))
                }
              } else {
                if (!is.na(value) && !is.nan(value)) {
                  temp_df <- rbind(temp_df, data.frame(
                    Metric = metric_name,
                    Participant = p_id,
                    Sensor = sensor,
                    Channel = ch,
                    Value = value
                  ))
                }
              }
            }
          }
        }
      }
      
      df_list[[metric_name]] <- temp_df
    }
  }
  
  # Combine all metrics into one dataframe
  if (length(df_list) > 0 && sum(sapply(df_list, nrow)) > 0) {
    combined_df <- do.call(rbind, df_list)
    combined_df$Metric <- factor(combined_df$Metric, levels = metrics)
    combined_df$Participant <- factor(combined_df$Participant, levels = c("p2", "p3", "p4", "p5"))
    combined_df$Sensor <- factor(combined_df$Sensor, levels = c("HG", "AgCl"))
    combined_df$Channel <- factor(combined_df$Channel, levels = c("ch1", "ch2", "ch3"))
  
    return(combined_df)
  } else {
    warning("No data found for the specified metrics")
    return(data.frame(
      Metric = character(),
      Participant = character(),
      Sensor = character(),
      Channel = character(),
      Value = numeric()
    ))
  }
}


# Prepare data for signal quality metrics
sigqual_metrics <- c("WaveEntropy", "SpecEntropy", "ZCR1", "ZCR2", "SpecFlatness", "Compress")
sigqual_df <- prepare_sigqual_data(data, sigqual_metrics)

formal_names <- c(
  WaveEntropy = "Waveform Entropy",
  SpecEntropy = "Spectrum Entropy",
  ZCR1 = "1st Order Zero Crossings",
  ZCR2 = "2nd Order Zero Crossings",
  SpecFlatness = "Spectral Detail",
  Compress = "Compression Detail"
)

# Add formal names to the data frame
sigqual_df$MetricFormal <- factor(formal_names[as.character(sigqual_df$Metric)], 
                                  levels = formal_names)

plot1 <- ggplot(sigqual_df, aes(x = Sensor, y = Value, fill = Sensor)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
  facet_grid(Participant ~ MetricFormal + Channel, scales = "free") +  # Use free scales for both x and y
  scale_fill_manual(values = c("HG" = "#1b9e77", "AgCl" = "#d95f02")) +
  labs(title = "Signal Quality Metrics by Participant, Sensor and Channel",
       y = "Value", x = "Sensor Type",
       fill = "Sensor Type") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95"),
    plot.title = element_text(size = 14, face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# 2. Group by Metric and Participant, facet by Sensor and Channel
plot2 <- ggplot(sigqual_df, aes(x = Participant, y = Value, fill = Participant)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
  facet_grid(Sensor ~ MetricFormal + Channel, scales = "free") +  # Use free scales for both x and y
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Signal Quality Metrics by Sensor, Participant and Channel",
       y = "Value", x = "Participant",
       fill = "Participant") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95"),
    plot.title = element_text(size = 14, face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# 3. Group by Metric and Channel, facet by Participant and Sensor
plot3 <- ggplot(sigqual_df, aes(x = Channel, y = Value, fill = Channel)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
  facet_grid(Participant ~ MetricFormal + Sensor, scales = "free") +  # Use free scales for both x and y
  scale_fill_brewer(palette = "Dark2") +
  labs(title = "Signal Quality Metrics by Participant, Sensor and Channel",
       y = "Value", x = "Channel",
       fill = "Channel") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95"),
    plot.title = element_text(size = 14, face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Create individual plots for each metric with white backgrounds
individual_plots <- list()
for (metric in sigqual_metrics) {
  metric_data <- sigqual_df %>% filter(Metric == metric)
  formal_name <- formal_names[metric]
  
  p <- ggplot(metric_data, aes(x = interaction(Sensor, Channel), y = Value, fill = Sensor)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
    facet_wrap(~ Participant, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = c("HG" = "#1b9e77", "AgCl" = "#d95f02")) +
    labs(title = paste(formal_name, "by Participant, Sensor and Channel"),
         y = formal_name, x = "Sensor-Channel",
         fill = "Sensor Type") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 12),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  individual_plots[[metric]] <- p
}

# Create additional plots with integrated participants
integrated_plots <- list()
for (metric in sigqual_metrics) {
  metric_data <- sigqual_df %>% filter(Metric == metric)
  formal_name <- formal_names[metric]
  
  # Plot comparing sensors across all participants with free scales
  p1 <- ggplot(metric_data, aes(x = Sensor, y = Value, fill = Sensor)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
    facet_wrap(~ Channel, scales = "free_y", ncol = 3) + # free_y scale
    scale_fill_manual(values = c("HG" = "#1b9e77", "AgCl" = "#d95f02")) +
    labs(title = paste(formal_name, "by Sensor and Channel (All Participants)"),
         y = formal_name, x = "Sensor Type",
         fill = "Sensor Type") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 12),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  # Plot comparing channels across all participants with free scales
  p2 <- ggplot(metric_data, aes(x = Channel, y = Value, fill = Channel)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
    facet_wrap(~ Sensor, scales = "free_y", ncol = 2) + # free_y scale
    scale_fill_brewer(palette = "Dark2") +
    labs(title = paste(formal_name, "by Channel and Sensor (All Participants)"),
         y = formal_name, x = "Channel",
         fill = "Channel") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 12),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  integrated_plots[[paste0(metric, "_by_sensor")]] <- p1
  integrated_plots[[paste0(metric, "_by_channel")]] <- p2
}

# Create a big integrated plot across all metrics
# All sensors and channels, no participant distinction
all_integrated <- ggplot(sigqual_df, aes(x = Sensor, y = Value, fill = Sensor)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
  facet_grid(Channel ~ MetricFormal, scales = "free_y") +  # free_y ensures each metric has its own scale
  scale_fill_manual(values = c("HG" = "#1b9e77", "AgCl" = "#d95f02")) +
  labs(title = "All Signal Quality Metrics by Sensor and Channel (All Participants)",
       y = "Value", x = "Sensor Type",
       fill = "Sensor Type") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 14, face = "bold")
  )

# Function to save plots with white backgrounds
save_plots <- function(plot_list, prefix = "sigqual_metric") {
  dir.create("plots", showWarnings = FALSE)
  
  # Save group plots
  ggsave(paste0("plots/", prefix, "_by_sensor.png"), plot1, width = 14, height = 10, bg = "white")
  ggsave(paste0("plots/", prefix, "_by_participant.png"), plot2, width = 14, height = 10, bg = "white")
  ggsave(paste0("plots/", prefix, "_by_channel.png"), plot3, width = 14, height = 10, bg = "white")
  
  # Save individual metric plots
  for (metric in names(individual_plots)) {
    formal_name <- tolower(gsub(" ", "_", formal_names[metric]))
    ggsave(paste0("plots/", prefix, "_", formal_name, ".png"), 
           individual_plots[[metric]], width = 10, height = 8, bg = "white")
  }
  
  # Save integrated plots
  for (plot_name in names(integrated_plots)) {
    base_name <- sub("^([^_]+)_.*$", "\\1", plot_name)
    formal_name <- tolower(gsub(" ", "_", formal_names[base_name]))
    suffix <- sub("^[^_]+_(.*)$", "\\1", plot_name)
    
    ggsave(paste0("plots/", prefix, "_", formal_name, "_integrated_", suffix, ".png"),
           integrated_plots[[plot_name]], width = 10, height = 8, bg = "white")
  }
  
  # Save the big integrated plot
  ggsave(paste0("plots/", prefix, "_all_integrated.png"), all_integrated, width = 16, height = 10, bg = "white")
}

# Display plots
print(plot1)
print(plot2)
print(plot3)

# Create a combined plot with all metrics (2x3 grid)
combined_plot <- (individual_plots$WaveEntropy | individual_plots$SpecEntropy | individual_plots$ZCR1) /
                 (individual_plots$ZCR2 | individual_plots$SpecFlatness | individual_plots$Compress) +
  plot_layout(guides = "collect") & 
  theme(
    legend.position = "bottom",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Save plots
save_plots(individual_plots)
ggsave("plots/sigqual_all_metrics.png", combined_plot, width = 18, height = 12, bg = "white")



prepare_similarity_data <- function(data, metrics) {
  df_list <- list()
  
  for (metric_name in metrics) {
    if (metric_name %in% names(data)) {
      temp_df <- data.frame()
      
      for (p_id in names(data[[metric_name]])) {
        # Focus on Xsensor data which contains matrices
        if ("Xsensor" %in% names(data[[metric_name]][[p_id]])) {
          for (ch in names(data[[metric_name]][[p_id]][["Xsensor"]])) {
            # Get the matrix for this channel
            mat_value <- data[[metric_name]][[p_id]][["Xsensor"]][[ch]]
            
            # Skip if value is NULL or empty
            if (is.null(mat_value) || length(mat_value) == 0) {
              next
            }
            
            # Debug print to understand the structure
            cat("Metric:", metric_name, "Participant:", p_id, "Channel:", ch, "\n")
            cat("  Type of mat_value:", class(mat_value), "\n")
            cat("  Structure:", str(mat_value), "\n")
            
            # Handle different types of data
            if (is.list(mat_value)) {
              # If it's a list, we need to extract numeric values recursively
              flat_values <- unlist(mat_value, recursive = TRUE, use.names = FALSE)
            } else {
              # Otherwise, just flatten the matrix or vector
              flat_values <- as.vector(mat_value)
            }
            
            # Make sure we have a numeric vector
            if (is.numeric(flat_values)) {
              # Skip NA/NaN values
              valid_values <- flat_values[!is.na(flat_values)]
              
              if (length(valid_values) > 0) {
                new_rows <- data.frame(
                  Metric = rep(metric_name, length(valid_values)),
                  Participant = rep(p_id, length(valid_values)),
                  Channel = rep(ch, length(valid_values)),
                  Value = valid_values
                )
                temp_df <- rbind(temp_df, new_rows)
              }
            } else {
              cat("Warning: Non-numeric data found in", metric_name, p_id, ch, "\n")
              cat("  Class:", class(flat_values), "\n")
            }
          }
        }
      }
      
      df_list[[metric_name]] <- temp_df
    }
  }
  
  # Combine all metrics into one dataframe
  if (length(df_list) > 0 && sum(sapply(df_list, nrow)) > 0) {
    combined_df <- do.call(rbind, df_list)
    combined_df$Metric <- factor(combined_df$Metric, levels = metrics)
    combined_df$Participant <- factor(combined_df$Participant, levels = c("p2", "p3", "p4", "p5"))
    combined_df$Channel <- factor(combined_df$Channel, levels = c("ch1", "ch2", "ch3"))
    combined_df$MetricFormal <- factor(sim_formal_names[as.character(combined_df$Metric)], 
                                       levels = sim_formal_names)
    
    return(combined_df)
  } else {
    warning("No data found for the specified metrics")
    return(data.frame(
      Metric = character(),
      MetricFormal = character(),
      Participant = character(),
      Channel = character(),
      Value = numeric()
    ))
  }
}

# Prepare similarity metrics data
sim_metrics <- c("CorrMat", "NrmseMat", "CosMat")
sim_df <- prepare_similarity_data(data, sim_metrics)

# Create violin plots for similarity metrics
# 1. Plot by participant and channel
sim_plot1 <- ggplot(sim_df, aes(x = Participant, y = Value, fill = Participant)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
  facet_grid(Channel ~ MetricFormal, scales = "free_y") +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "HG-AgCl Sensor Similarity Metrics by Participant and Channel",
       y = "Value", x = "Participant",
       fill = "Participant") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95"),
    plot.title = element_text(size = 14, face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# 2. Plot by channel comparing participants
sim_plot2 <- ggplot(sim_df, aes(x = Channel, y = Value, fill = Channel)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
  facet_grid(Participant ~ MetricFormal, scales = "free_y") +
  scale_fill_brewer(palette = "Dark2") +
  labs(title = "HG-AgCl Sensor Similarity Metrics by Channel and Participant",
       y = "Value", x = "Channel",
       fill = "Channel") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95"),
    plot.title = element_text(size = 14, face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Create individual plots for each similarity metric
sim_individual_plots <- list()
for (metric in sim_metrics) {
  metric_data <- sim_df %>% filter(Metric == metric)
  formal_name <- sim_formal_names[metric]
  
  p <- ggplot(metric_data, aes(x = Channel, y = Value, fill = Channel)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
    facet_wrap(~ Participant, scales = "free_y", ncol = 2) +
    scale_fill_brewer(palette = "Dark2") +
    labs(title = paste(formal_name, "by Channel and Participant - HG-AgCl Sensor Comparison"),
         y = formal_name, x = "Channel",
         fill = "Channel") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 12),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  sim_individual_plots[[metric]] <- p
}

# Create integrated plots (no participant distinction)
sim_integrated_plots <- list()
for (metric in sim_metrics) {
  metric_data <- sim_df %>% filter(Metric == metric)
  formal_name <- sim_formal_names[metric]
  
  p <- ggplot(metric_data, aes(x = Channel, y = Value, fill = Channel)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
    scale_fill_brewer(palette = "Dark2") +
    labs(title = paste(formal_name, "by Channel (All Participants) - HG-AgCl Sensor Comparison"),
         y = formal_name, x = "Channel",
         fill = "Channel") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 12),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  sim_integrated_plots[[metric]] <- p
}

# Combined plot of all similarity metrics
sim_combined_plot <- (sim_individual_plots$CorrMat | sim_individual_plots$NrmseMat | sim_individual_plots$CosMat) +
  plot_layout(guides = "collect") & 
  theme(
    legend.position = "bottom",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Function to save similarity plots
save_similarity_plots <- function(prefix = "similarity_metric") {
  dir.create("plots", showWarnings = FALSE)
  
  # Save main plots
  ggsave(paste0("plots/", prefix, "_by_participant.png"), sim_plot1, width = 12, height = 8, bg = "white")
  ggsave(paste0("plots/", prefix, "_by_channel.png"), sim_plot2, width = 12, height = 8, bg = "white")
  
  # Save individual metric plots
  for (metric in sim_metrics) {
    formal_name <- tolower(gsub(" ", "_", sim_formal_names[metric]))
    ggsave(paste0("plots/", prefix, "_", formal_name, ".png"), 
           sim_individual_plots[[metric]], width = 10, height = 8, bg = "white")
  }
  
  # Save integrated plots
  for (metric in sim_metrics) {
    formal_name <- tolower(gsub(" ", "_", sim_formal_names[metric]))
    ggsave(paste0("plots/", prefix, "_", formal_name, "_integrated.png"),
           sim_integrated_plots[[metric]], width = 8, height = 6, bg = "white")
  }
  
  # Save combined plot
  ggsave(paste0("plots/", prefix, "_all_combined.png"), sim_combined_plot, width = 16, height = 6, bg = "white")
}

# Display similarity plots
print(sim_plot1)
print(sim_plot2)
print(sim_combined_plot)

# Save similarity plots
save_similarity_plots()