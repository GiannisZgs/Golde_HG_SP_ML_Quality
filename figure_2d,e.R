library(jsonlite)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
library(grid)
library(cowplot)
library(ggrepel)
library(viridis)
library(ggbeeswarm)

# Load data
data_path <- "/home/giannis/Documents/ECG HG paper/results_data/participant_id_results.json"
data <- fromJSON(txt = data_path, simplifyVector = FALSE)

# Create output directory
output_dir <- "R_figures/figures_2d_e/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define sensors and channels
channels <- c("ch1", "ch2", "ch3")
channel_labels <- c("ch1" = "Lead 1", "ch2" = "Lead 2", "ch3" = "Lead 3")
sensors <- c("AgCl", "HG")
sensor_labels <- c("AgCl" = "AgCl", "HG" = "PPHG")

# Define colors for sensors
sensor_colors <- c("AgCl" = "#d9020d", "PPHG" = "#0066CC")

# Function to extract MSE data
extract_mse_data <- function() {
  mse_data <- data.frame()
  
  # Check if MSE data exists
  if (!"mse" %in% names(data)) {
    warning("No MSE data found in the dataset")
    return(mse_data)
  }
  
  # Extract MSE data for each channel and participant
  for (channel in channels) {
    
    if (channel %in% names(data$mse)) {
      channel_data <- data$mse[[channel]]
      
      for (participant in names(channel_data)) {
        if (participant != "_fieldnames") {
          value <- as.numeric(channel_data[[participant]])
          
          if (!is.null(value) && is.finite(value)) {
            mse_data <- rbind(mse_data, data.frame(
              Metric = "MSE",
              Channel = channel,
              ChannelLabel = channel_labels[channel],
              Participant = participant,
              Value = value
            ))
          }
        }
      }
    }
  }
  
  return(mse_data)
}

# Function to extract cross-correlation data
extract_xcorr_data <- function() {
  xcorr_data <- data.frame()
  
  # Check if XCORR data exists
  if (!"xcorr" %in% names(data)) {
    warning("No XCORR data found in the dataset")
    return(xcorr_data)
  }
  
  # Extract XCORR data for each channel and participant
  for (channel in channels) {
    
    if (channel %in% names(data$xcorr)) {
      channel_data <- data$xcorr[[channel]]
      
      for (participant in names(channel_data)) {
        if (participant != "_fieldnames") {
          value <- as.numeric(channel_data[[participant]])
          
          if (!is.null(value) && is.finite(value)) {
            xcorr_data <- rbind(xcorr_data, data.frame(
              Metric = "XCORR",
              Channel = channel,
              ChannelLabel = channel_labels[channel],
              Participant = participant,
              Value = value
            ))
          }
        }
      }
    }
  }
  
  return(xcorr_data)
}

# Function to extract ECG feature data
extract_features_data <- function() {
  features_data <- data.frame()
  
  # Check if features data exists
  if (!"features" %in% names(data)) {
    warning("No features data found in the dataset")
    return(features_data)
  }
  
  # Define feature metrics
  feature_metrics <- c("p_wave_amp", "peak_to_peak", "t_wave_amp")
  
  # Extract features data for each participant, sensor, channel, and metric
  for (participant in names(data$features)) {
    if (participant == "_fieldnames") {
      next  # Skip metadata field
    }
    
    participant_data <- data$features[[participant]]
    
    for (sensor in sensors) {
      if (!sensor %in% names(participant_data)) {
        next
      }
      
      sensor_data <- participant_data[[sensor]]
      
      for (channel in channels) {
        if (!channel %in% names(sensor_data)) {
          next
        }
        
        channel_data <- sensor_data[[channel]]
        
        for (metric in feature_metrics) {
          if (!metric %in% names(channel_data)) {
            next
          }
          
          value <- as.numeric(channel_data[[metric]])
          
          if (!is.null(value) && is.finite(value)) {
            features_data <- rbind(features_data, data.frame(
              MetricType = "Feature",
              Metric = metric,
              Sensor = sensor,
              SensorLabel = sensor_labels[sensor],
              Channel = channel,
              ChannelLabel = channel_labels[channel],
              Participant = participant,
              Value = value
            ))
          }
        }
      }
    }
  }
  
  return(features_data)
}

# Extract all data
mse_data <- extract_mse_data()
xcorr_data <- extract_xcorr_data()
features_data <- extract_features_data()

# Print summary of extracted data
cat("Extracted", nrow(mse_data), "MSE data points\n")
cat("Extracted", nrow(xcorr_data), "XCORR data points\n")
cat("Extracted", nrow(features_data), "feature data points\n")

# Prepare MSE and XCORR data for plotting (add Sensor label for consistency)
mse_data$Sensor <- NA
mse_data$SensorLabel <- NA
xcorr_data$Sensor <- NA
xcorr_data$SensorLabel <- NA

mse_xcorr_color <- "#5D3FD3"  # Deep purple

# Function to create plots with purple color for MSE and XCORR
create_plot <- function(data, plot_type, x_var = "Channel", color_var = "Sensor", 
                        fill_var = "Sensor", title = "", y_label = "Value") {
  
  # Set up the base plot
  p <- ggplot(data, aes_string(x = x_var, y = "Value"))
  
  # Determine if this is MSE/XCORR data (no sensor specified) or feature data
  is_mse_xcorr <- is.null(color_var) || is.na(data$Sensor[1])
  
  # Add appropriate geoms based on plot type
  if (plot_type == "boxplot") {
    if (is_mse_xcorr) {
      # For MSE and XCORR, use the deep purple color
      p <- p + 
        geom_boxplot(fill = mse_xcorr_color, alpha = 0.7, outlier.shape = NA)
      
      # Add jitter with matching color
      if (nrow(data) < 500) {
        p <- p + geom_jitter(color = mse_xcorr_color, width = 0.2, height = 0, alpha = 0.6, size = 1.5)
      }
    } else {
      # For feature data, use the sensor colors
      p <- p + 
        geom_boxplot(aes_string(fill = fill_var), alpha = 0.7, outlier.shape = NA)
      
      if (nrow(data) < 500) {
        p <- p + geom_jitter(aes_string(color = color_var), width = 0.2, height = 0, alpha = 0.6, size = 1.5)
      }
      
      # Add fill scale if fill variable is provided
      if (!is.null(fill_var)) {
        if (fill_var == "Sensor" || fill_var == "SensorLabel") {
          p <- p + scale_fill_manual(values = sensor_colors)
        }
      }
    }
    
  } else if (plot_type == "violin") {
    if (is_mse_xcorr) {
      # For MSE and XCORR, use the deep purple color
      p <- p + 
        geom_violin(fill = mse_xcorr_color, alpha = 0.7, scale = "width", trim = TRUE)
      
      if (nrow(data) < 500) {
        p <- p + geom_jitter(color = mse_xcorr_color, width = 0.1, height = 0, alpha = 0.6, size = 1)
      }
    } else {
      # For feature data, use the sensor colors
      p <- p + 
        geom_violin(aes_string(fill = fill_var), alpha = 0.7, scale = "width", trim = TRUE)
      
      if (nrow(data) < 500) {
        p <- p + geom_jitter(aes_string(color = color_var), width = 0.1, height = 0, alpha = 0.6, size = 1)
      }
      
      if (!is.null(fill_var)) {
        if (fill_var == "Sensor" || fill_var == "SensorLabel") {
          p <- p + scale_fill_manual(values = sensor_colors)
        }
      }
    }
    
  } else if (plot_type == "beeswarm") {
    if (is_mse_xcorr) {
      # For MSE and XCORR, use the deep purple color
      p <- p + 
        geom_beeswarm(color = mse_xcorr_color, size = 1.5, alpha = 0.7, cex = 1.5)
    } else {
      # For feature data, use the sensor colors
      p <- p + 
        geom_beeswarm(aes_string(color = color_var), size = 1.5, alpha = 0.7, cex = 1.5)
      
      if (!is.null(color_var)) {
        if (color_var == "Sensor" || color_var == "SensorLabel") {
          p <- p + scale_color_manual(values = sensor_colors)
        }
      }
    }
  }
  
  # Complete the plot with labels and theme
  p <- p +
    labs(
      title = title,
      x = "",
      y = y_label
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 25),
      axis.text.y = element_text(size = 25),
      axis.title = element_text(size = 25),
      strip.text = element_text(size = 25),
      plot.title = element_text(size = 25),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  return(p)
}

# Create MSE plots
if (nrow(mse_data) > 0) {
  # Format channel labels
  mse_data$ChannelLabel <- factor(mse_data$ChannelLabel, levels = channel_labels)
  
  mse_boxplot <- create_plot(mse_data, "boxplot", "ChannelLabel", NULL, NULL, 
                            "Normalized Root Mean Squared Error", "NRMSE")
  mse_violin <- create_plot(mse_data, "violin", "ChannelLabel", NULL, NULL,
                           "Normalized Root Mean Squared Error", "NRMSE")
  mse_beeswarm <- create_plot(mse_data, "beeswarm", "ChannelLabel", NULL, NULL,
                             "Normalized Root Mean Squared Error", "NRMSE")
  
  # Save plots
  ggsave(paste0(output_dir, "mse_boxplot.png"), mse_boxplot, width = 8, height = 6, dpi = 300, bg = "white")
  ggsave(paste0(output_dir, "mse_violin.png"), mse_violin, width = 8, height = 6, dpi = 300, bg = "white")
  ggsave(paste0(output_dir, "mse_beeswarm.png"), mse_beeswarm, width = 8, height = 6, dpi = 300, bg = "white")
  
  cat("Created MSE plots with deep purple color\n")
} else {
  cat("No MSE data available for plotting\n")
}

# Create XCORR plots
if (nrow(xcorr_data) > 0) {
  # Format channel labels
  xcorr_data$ChannelLabel <- factor(xcorr_data$ChannelLabel, levels = channel_labels)
  
  xcorr_boxplot <- create_plot(xcorr_data, "boxplot", "ChannelLabel", NULL, NULL,
                              "Cross-Correlation", "Cross-Correlation")
  xcorr_violin <- create_plot(xcorr_data, "violin", "ChannelLabel", NULL, NULL,
                             "Cross-Correlation", "Cross-Correlation")
  xcorr_beeswarm <- create_plot(xcorr_data, "beeswarm", "ChannelLabel", NULL, NULL,
                               "Cross-Correlation", "Cross-Correlation")
  
  # Save plots
  ggsave(paste0(output_dir, "xcorr_boxplot.png"), xcorr_boxplot, width = 8, height = 6, dpi = 300, bg = "white")
  ggsave(paste0(output_dir, "xcorr_violin.png"), xcorr_violin, width = 8, height = 6, dpi = 300, bg = "white")
  ggsave(paste0(output_dir, "xcorr_beeswarm.png"), xcorr_beeswarm, width = 8, height = 6, dpi = 300, bg = "white")
  
  cat("Created XCORR plots with deep purple color\n")
} else {
  cat("No XCORR data available for plotting\n")
}

# Create feature plots
if (nrow(features_data) > 0) {
  # Get unique metrics
  unique_metrics <- unique(features_data$Metric)
  
  # Format channel and sensor labels
  features_data$ChannelLabel <- factor(features_data$ChannelLabel, levels = channel_labels)
  features_data$SensorLabel <- factor(features_data$SensorLabel, levels = unname(sensor_labels))
  
  # Create nice labels for metrics
  metric_labels <- c(
    "p_wave_amp" = "P-Wave Amplitude",
    "peak_to_peak" = "Peak-to-Peak Amplitude",
    "t_wave_amp" = "T-Wave Amplitude"
  )
  
  # Create plots for each feature metric
  for (metric in unique_metrics) {
    # Filter data for current metric
    metric_data <- features_data[features_data$Metric == metric, ]
    
    # Skip if no data
    if (nrow(metric_data) == 0) {
      next
    }
    
    # Create display name for the metric
    metric_display <- ifelse(metric %in% names(metric_labels), 
                            metric_labels[metric], 
                            gsub("_", " ", capitalize(metric)))
    
    # Create plots comparing sensors for each channel
    boxplot <- create_plot(metric_data, "boxplot", "ChannelLabel", "SensorLabel", "SensorLabel",
                          paste0(metric_display), "Amplitude (\u00B5V)")
    
    violin <- create_plot(metric_data, "violin", "ChannelLabel", "SensorLabel", "SensorLabel",
                         paste0(metric_display), "Amplitude (\u00B5V)")
    
    beeswarm <- create_plot(metric_data, "beeswarm", "ChannelLabel", "SensorLabel", NULL,
                           paste0(metric_display), "Amplitude (\u00B5V)")
    
    # Add facet by sensor for clearer comparison
    boxplot <- boxplot + facet_wrap(~SensorLabel)
    violin <- violin + facet_wrap(~SensorLabel)
    beeswarm <- beeswarm + facet_wrap(~SensorLabel)
    
    # Save plots
    clean_name <- gsub(" ", "_", tolower(metric))
    ggsave(paste0(output_dir, clean_name, "_boxplot.png"), boxplot, width = 10, height = 6, dpi = 300, bg = "white")
    ggsave(paste0(output_dir, clean_name, "_violin.png"), violin, width = 10, height = 6, dpi = 300, bg = "white")
    ggsave(paste0(output_dir, clean_name, "_beeswarm.png"), beeswarm, width = 10, height = 6, dpi = 300, bg = "white")
    
    cat("Created plots for", metric_display, "\n")
  }
  
  # Create combined features plot
  # First, create a human-readable metric name
  features_data$MetricDisplay <- sapply(features_data$Metric, function(m) {
    ifelse(m %in% names(metric_labels), metric_labels[m], gsub("_", " ", capitalize(m)))
  })
  
  # Create combined boxplot with all features
  combined_boxplot <- ggplot(features_data, aes(x = ChannelLabel, y = Value, fill = SensorLabel)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    scale_fill_manual(values = sensor_colors) +
    labs(
      title = "",
      x = "",
      y = "Amplitude (\u00B5V)",
      fill = "Sensor"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 25),
      axis.text.y = element_text(size = 25),
      axis.title = element_text(size = 25),
      plot.title = element_text(size = 25),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      strip.text = element_text(size = 25)
    ) +
    facet_grid(MetricDisplay ~ SensorLabel, scales = "free_y")
  
  ggsave(paste0(output_dir, "all_features_boxplot.png"), combined_boxplot, 
         width = 10, height = 12, dpi = 300, bg = "white")
  
  cat("Created combined features plot\n")
} else {
  cat("No features data available for plotting\n")
}

# Helper function to capitalize strings
capitalize <- function(x) {
  s <- strsplit(x, "_")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "", collapse = " ")
}

combined_boxplot_with_points <- ggplot(features_data, aes(x = ChannelLabel, y = Value, fill = SensorLabel)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = SensorLabel), width = 0.2, height = 0, alpha = 0.6, size = 1.5) +
  scale_fill_manual(values = sensor_colors) +
  scale_color_manual(values = sensor_colors) +
  labs(
    title = "",
    x = "",
    y = "Amplitude (\u00B5V)",
    fill = "Sensor",
    color = "Sensor"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 25),
    axis.text.y = element_text(size = 25),
    axis.title = element_text(size = 25),
    plot.title = element_text(size = 25),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 25)
  ) +
  facet_grid(MetricDisplay ~ SensorLabel, scales = "free_y")

# Create combined violin plot with all features
combined_violin <- ggplot(features_data, aes(x = ChannelLabel, y = Value, fill = SensorLabel)) +
  geom_violin(alpha = 0.7, scale = "width", trim = TRUE) +
  geom_jitter(aes(color = SensorLabel), width = 0.1, height = 0, alpha = 0.6, size = 1) +
  scale_fill_manual(values = sensor_colors) +
  scale_color_manual(values = sensor_colors) +
  labs(
    title = "",
    x = "",
    y = "Amplitude (\u00B5V)",
    fill = "Sensor",
    color = "Sensor"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 25),
    axis.text.y = element_text(size = 25),
    axis.title = element_text(size = 25),
    plot.title = element_text(size = 25),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 25)
  ) +
  facet_grid(MetricDisplay ~ SensorLabel, scales = "free_y")

# Create combined beeswarm plot with all features
combined_beeswarm <- ggplot(features_data, aes(x = ChannelLabel, y = Value, color = SensorLabel)) +
  geom_beeswarm(size = 1.5, alpha = 0.7, cex = 1.5) +
  scale_color_manual(values = sensor_colors) +
  labs(
    title = "",
    x = "",
    y = "Amplitude (\u00B5V)",
    color = "Sensor"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 25),
    axis.text.y = element_text(size = 25),
    axis.title = element_text(size = 25),
    plot.title = element_text(size = 25),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 25)
  ) +
  facet_grid(MetricDisplay ~ SensorLabel, scales = "free_y")

# Save the combined plots
ggsave(paste0(output_dir, "all_features_boxplot_with_points.png"), combined_boxplot_with_points, 
       width = 10, height = 12, dpi = 300, bg = "white")

ggsave(paste0(output_dir, "all_features_violin.png"), combined_violin, 
       width = 10, height = 12, dpi = 300, bg = "white")

ggsave(paste0(output_dir, "all_features_beeswarm.png"), combined_beeswarm, 
       width = 10, height = 12, dpi = 300, bg = "white")

cat("Created additional combined feature plots (boxplot with points, violin, beeswarm)\n")

cat("All plots have been saved to:", output_dir, "\n")