library(jsonlite)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(scales)
library(grid)
library(ggforce)
library(RColorBrewer)

# Load data
data_path <- "/home/giannis/Documents/ECG HG paper/results_data/participant_id_results.json"
data <- fromJSON(txt = data_path, simplifyVector = FALSE)

# Create output directory
output_dir <- "R_figures/figure_2f/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define sensors and channels
channels <- c("ch1", "ch2", "ch3")
channel_labels <- c("ch1" = "Lead 1", "ch2" = "Lead 2", "ch3" = "Lead 3")
sensors <- c("AgCl", "HG")
sensor_labels <- c("AgCl" = "AgCl", "HG" = "PPHG")

# Define colors for sensors and metrics
sensor_colors <- c("AgCl" = "#d9020d", "PPHG" = "#0066CC")
mse_xcorr_color <- "#5D3FD3"  # Deep purple

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
              Metric = "NRMSE",
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
              Metric = "Cross-Correlation",
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

# Create nice labels for metrics
metric_labels <- c(
  "p_wave_amp" = "P-Wave Amplitude",
  "peak_to_peak" = "Peak-to-Peak",
  "t_wave_amp" = "T-Wave Amplitude",
  "NRMSE" = "NRMSE",
  "Cross-Correlation" = "Cross-Correlation"
)

# Add sensor information to MSE and XCORR
mse_data$Sensor <- "MSE"
mse_data$SensorLabel <- "MSE"
xcorr_data$Sensor <- "XCORR" 
xcorr_data$SensorLabel <- "XCORR"

# Compute summary statistics for features (mean and 95% CI for each metric, sensor, and channel)
features_summary <- features_data %>%
  group_by(Metric, Sensor, SensorLabel, Channel, ChannelLabel) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    N = n(),
    SE = SD / sqrt(N),
    CI_lower = Mean - qt(0.975, N-1) * SE,
    CI_upper = Mean + qt(0.975, N-1) * SE,
    .groups = "drop"
  ) %>%
  mutate(
    MetricLabel = ifelse(Metric %in% names(metric_labels), 
                         metric_labels[Metric], 
                         gsub("_", " ", Metric))
  )

# Compute summary statistics for MSE
mse_summary <- mse_data %>%
  group_by(Metric, Channel, ChannelLabel) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    N = n(),
    SE = SD / sqrt(N),
    CI_lower = Mean - qt(0.975, N-1) * SE,
    CI_upper = Mean + qt(0.975, N-1) * SE,
    .groups = "drop"
  ) %>%
  mutate(
    MetricLabel = ifelse(Metric %in% names(metric_labels), 
                         metric_labels[Metric], 
                         gsub("_", " ", Metric)),
    Sensor = "MSE",
    SensorLabel = "MSE"
  )

# Compute summary statistics for XCORR
xcorr_summary <- xcorr_data %>%
  group_by(Metric, Channel, ChannelLabel) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    N = n(),
    SE = SD / sqrt(N),
    CI_lower = Mean - qt(0.975, N-1) * SE,
    CI_upper = Mean + qt(0.975, N-1) * SE,
    .groups = "drop"
  ) %>%
  mutate(
    MetricLabel = ifelse(Metric %in% names(metric_labels), 
                         metric_labels[Metric], 
                         gsub("_", " ", Metric)),
    Sensor = "XCORR",
    SensorLabel = "XCORR"
  )

# Combine all summaries
all_summaries <- bind_rows(features_summary, mse_summary, xcorr_summary)

normalize_data <- function(data) {
  # For each metric, apply appropriate normalization
  normalized_data <- data %>%
    group_by(Metric) %>%
    mutate(
      # For NRMSE, convert to 1-NRMSE so higher is better (like other metrics)
      Mean = ifelse(Metric == "NRMSE", 1-Mean, Mean),
      CI_lower = ifelse(Metric == "NRMSE", 1-CI_upper, CI_lower),  # Reverse CIs when inverting the metric
      CI_upper = ifelse(Metric == "NRMSE", 1-CI_lower, CI_upper)
    ) %>%
    ungroup()
  
  # For each metric, normalize to 0-1 scale based on maximum value
  normalized_data <- normalized_data %>%
    group_by(Metric) %>%
    mutate(
      max_val = max(Mean, na.rm = TRUE),
      Mean_norm = Mean / max_val,
      CI_lower_norm = CI_lower / max_val,
      CI_upper_norm = CI_upper / max_val
    ) %>%
    ungroup()
  
  return(normalized_data)
}

print("Unique Metrics in all_summaries_norm:")
print(unique(all_summaries_norm$Metric))


create_circular_barplot <- function(sensor_data, sensor_name, palette = NULL) {
  # Define the number of metrics and channels
  # Convert to character first to avoid factor level issues
  sensor_data$Metric <- as.character(sensor_data$Metric)
  unique_metrics <- unique(sensor_data$Metric)
  n_metrics <- length(unique_metrics)
  n_channels <- length(unique(sensor_data$Channel))
  
  # Debug print
  cat("Number of unique metrics in", sensor_name, "data:", n_metrics, "\n")
  cat("Unique metrics in", sensor_name, "data:", paste(unique_metrics, collapse=", "), "\n")
  
  # Create a factor for metric to control the order - include ONLY the 5 correct metrics
  metric_order <- c("p_wave_amp", "peak_to_peak", "t_wave_amp", "NRMSE", "Cross-Correlation")
  
  # Ensure we only use metrics present in the data
  metric_order <- metric_order[metric_order %in% unique_metrics]
  
  # Force Metric to be a factor with ONLY the correct levels
  sensor_data$Metric <- factor(sensor_data$Metric, levels = metric_order)
 
  # Sort data to ensure correct ordering in plot
  sensor_data <- sensor_data[order(sensor_data$Channel, sensor_data$Metric),]
  
  # Define angle for each bar
  total_bars <- nrow(sensor_data)
  sensor_data$id <- seq_len(total_bars)
  
  # Define colors based on metric type
  if (is.null(palette)) {
    # Default colors if no palette is provided
    feature_colors <- viridis(3)
    metric_colors <- c(feature_colors, "purple", "darkblue")
  } else {
    # Match palette colors to metric_order (only use what we need)
    metric_colors <- palette[1:length(metric_order)]
  }
  
  # Set the limit for the radial axis (all metrics normalized to 1)
  y_limit <- 1.1  # Set to 1.1 to leave some space for labels
  
  # Calculate the angles for bars
  sensor_data <- sensor_data %>%
    mutate(
      # Adjust angle to place bars in the right position
      angle = (id - 0.5) * (360 / total_bars),
      # For labels, use the middle of the bar
      label_angle = (id - 0.5) * (360 / total_bars),
      hjust = ifelse(label_angle > 90 & label_angle < 270, 1, 0),
      vjust = ifelse(label_angle > 180 & label_angle < 360, 0, 1),
      label_angle = ifelse(label_angle > 90 & label_angle < 270, label_angle + 180, label_angle)
    )
  
  # Add channel boundary lines
  channel_boundaries <- seq(0, 360, by = 360/n_channels)
  
  # Create channel sectors data for highlighting
  channel_sectors <- data.frame()
  for (i in 1:n_channels) {
    start_angle <- channel_boundaries[i]
    end_angle <- channel_boundaries[i+1]
    if (is.na(end_angle)) end_angle <- 360
    
    # Add sector data
    channel_sectors <- rbind(channel_sectors, data.frame(
      channel = channels[i],
      channel_label = channel_labels[channels[i]],
      start = start_angle,
      end = end_angle,
      # Alternate colors for visual distinction
      fill_color = ifelse(i %% 2 == 0, "gray95", "white")
    ))
  }
  
  # Create the labels that match our actual metrics
  metric_display_labels <- c(
    "p_wave_amp" = "P-Wave Amplitude",
    "peak_to_peak" = "Peak-to-Peak",
    "t_wave_amp" = "T-Wave Amplitude",
    "NRMSE" = "Signal Quality (1-NRMSE)",
    "Cross-Correlation" = "Cross-Correlation"
  )
  
  # Only use labels for metrics we actually have
  metric_labels_to_use <- metric_display_labels[metric_order]
  
  # Create the circular barplot
  p <- ggplot() +
    # First add channel sectors as background
    geom_rect(data = channel_sectors, 
              aes(xmin = start, xmax = end, ymin = 0, ymax = y_limit, fill = fill_color),
              alpha = 0.3) +
    scale_fill_identity() +
    
    # Add the circular barplot bars
    geom_col(data = sensor_data, 
             aes(x = angle, y = Mean_norm, fill = Metric), 
             position = "identity", width = 360/total_bars * 0.9) +
    
    # Add error bars
    geom_errorbar(data = sensor_data,
                 aes(x = angle, ymin = CI_lower_norm, ymax = CI_upper_norm), 
                 width = 360/total_bars * 0.4, color = "black") +
    
    # Add coordinate system and other components
    coord_polar(start = 0) +
    scale_y_continuous(limits = c(0, y_limit), breaks = seq(0, 1, by = 0.25)) +
    scale_fill_manual(values = metric_colors, 
                     name = "Metric",
                     labels = metric_labels_to_use) +
    labs(
      title = paste0(sensor_name, " Sensor Metrics Overview"),
      x = NULL,
      y = NULL,
      caption = "All metrics normalized to 0-1 scale. Bars represent mean values with 95% confidence intervals."
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 10),
      axis.ticks.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "gray80", size = 0.5),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.caption = element_text(size = 10, hjust = 0.5),
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  
  # Add channel labels as annotations at appropriate positions
  channel_labels_df <- data.frame(
    angle = seq(0, 360, length.out = n_channels + 1)[1:n_channels] + 360/(n_channels*2),
    label = channel_labels,
    r = y_limit * 1.05  # Position slightly outside the plot
  )
  
  # Add channel labels
  for (i in 1:n_channels) {
    angle_rad <- channel_labels_df$angle[i] * pi / 180
    x_pos <- sin(angle_rad) * channel_labels_df$r[i]
    y_pos <- cos(angle_rad) * channel_labels_df$r[i]
    
    p <- p + 
      annotate("text", 
               x = channel_labels_df$angle[i], 
               y = y_limit * 1.05, 
               label = channel_labels_df$label[i], 
               size = 5, 
               fontface = "bold")
  }
  
  # Add channel boundary lines
  for (angle in channel_boundaries) {
    p <- p + geom_vline(xintercept = angle, 
                        linetype = "dashed", 
                        color = "gray50", 
                        size = 1)
  }
  
  # Add metric description text at the bottom of the plot
  p <- p + labs(
    subtitle = paste0("NRMSE values inverted (1-NRMSE) so higher values represent better performance")
  ) +
    theme(
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40")
    )
  
  return(p)
}# Filter data for each sensor

# For AgCl data, include only the feature data for AgCl and the NRMSE/Cross-Correlation metrics
agcl_data <- all_summaries_norm %>% 
  filter((Sensor == "AgCl" & Metric %in% c("p_wave_amp", "peak_to_peak", "t_wave_amp")) | 
         (Metric %in% c("NRMSE", "Cross-Correlation")))

# For PPHG data, include only the feature data for HG and the NRMSE/Cross-Correlation metrics  
pphg_data <- all_summaries_norm %>% 
  filter((Sensor == "HG" & Metric %in% c("p_wave_amp", "peak_to_peak", "t_wave_amp")) | 
         (Metric %in% c("NRMSE", "Cross-Correlation")))

# Print the unique metrics for debugging
print("Unique metrics in agcl_data:")
print(unique(agcl_data$Metric))
print("Unique metrics in pphg_data:")
print(unique(pphg_data$Metric))

# Use droplevels to remove any unused factor levels
agcl_data$Metric <- droplevels(as.factor(agcl_data$Metric))
pphg_data$Metric <- droplevels(as.factor(pphg_data$Metric))


# Update metric order to ensure we only have the 5 metrics
metric_order <- c("p_wave_amp", "peak_to_peak", "t_wave_amp", "NRMSE", "Cross-Correlation")

# Create a 5-color palette for the 5 metrics
metric_palette <- c(
  "#E41A1C",  # red - p-wave amplitude
  "#377EB8",  # blue - peak-to-peak
  "#4DAF4A",  # green - t-wave amplitude
  "#984EA3",  # purple - NRMSE
  "#FF7F00"   # orange - Cross-Correlation
)

# Generate the circular barplots
agcl_circular <- create_circular_barplot(agcl_data, "AgCl", metric_palette)
pphg_circular <- create_circular_barplot(pphg_data, "PPHG", metric_palette)

# Save the plots
ggsave(paste0(output_dir, "agcl_circular_barplot.png"), agcl_circular, 
       width = 10, height = 10, dpi = 300, bg = "white")
ggsave(paste0(output_dir, "pphg_circular_barplot.png"), pphg_circular, 
       width = 10, height = 10, dpi = 300, bg = "white")

# Create a combined plot with both circular barplots
combined_plot <- cowplot::plot_grid(
  agcl_circular, pphg_circular, 
  labels = c("A", "B"), 
  ncol = 2, 
  align = "h"
)

# Save the combined plot
ggsave(paste0(output_dir, "combined_circular_barplots.png"), combined_plot, 
       width = 16, height = 8, dpi = 300, bg = "white")

cat("Circular barplots created and saved to:", output_dir, "\n")