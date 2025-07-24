library(jsonlite)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(scales)
library(grid)
library(ggforce)
library(RColorBrewer)
library(cowplot)

# Load data
data_path <- "/home/giannis/Documents/ECG HG paper/results_data/participant_id_results.json"
data <- fromJSON(txt = data_path, simplifyVector = FALSE)

participant_id_path <- "/home/giannis/Documents/ECG HG paper/results_data/participant_id_results.json"
channel_id_path <- "/home/giannis/Documents/ECG HG paper/results_data/channel_id_results.json"

# Load participant identification data
participant_id_data <- fromJSON(txt = participant_id_path, simplifyVector = FALSE)

# Load channel identification data
channel_id_data <- fromJSON(txt = channel_id_path, simplifyVector = FALSE)

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

extract_participant_id_metrics <- function() {
  metrics_data <- data.frame()
  
  for (sensor in c("AgCl", "HG")) {
    if (!sensor %in% names(participant_id_data$tsne_vis)) {
      next
    }
    
    sensor_metrics <- participant_id_data$tsne_vis[[sensor]]
    
    for (channel in c("ch1", "ch2", "ch3")) {
      if (!channel %in% names(sensor_metrics)) {
        next
      }
      
      channel_metrics <- sensor_metrics[[channel]]
      
      # Extract the metrics
      acc <- as.numeric(channel_metrics$acc)
      micro_f1 <- as.numeric(channel_metrics$micro_f1)
      macro_f1 <- as.numeric(channel_metrics$macro_f1)
      
      # Add to dataframe
      metrics_data <- rbind(metrics_data, data.frame(
        Sensor = sensor,
        Channel = channel,
        Metric = "participant_acc",
        Mean = acc,
        CI_lower = acc ,  # Approximation for CI
        CI_upper = acc  # Approximation for CI, capped at 1.0
      ))
      
      metrics_data <- rbind(metrics_data, data.frame(
        Sensor = sensor,
        Channel = channel,
        Metric = "participant_micro_f1",
        Mean = micro_f1,
        CI_lower = micro_f1,  # Approximation for CI
        CI_upper = micro_f1  # Approximation for CI, capped at 1.0
      ))
      
      metrics_data <- rbind(metrics_data, data.frame(
        Sensor = sensor,
        Channel = channel,
        Metric = "participant_macro_f1",
        Mean = macro_f1,
        CI_lower = macro_f1,  # Approximation for CI
        CI_upper = macro_f1  # Approximation for CI, capped at 1.0
      ))
    }
  }
  
  return(metrics_data)
}

# Extract channel identification metrics for each sensor (not channel-specific)
extract_channel_id_metrics <- function() {
  metrics_data <- data.frame()
  
  for (sensor in c("AgCl", "HG")) {
    if (!is.null(channel_id_data$tsne_vis[[sensor]])) {
      # Extract the metrics
      acc <- as.numeric(channel_id_data$tsne_vis[[sensor]]$acc)
      micro_f1 <- as.numeric(channel_id_data$tsne_vis[[sensor]]$micro_f1)
      macro_f1 <- as.numeric(channel_id_data$tsne_vis[[sensor]]$macro_f1)
      
      # Add to dataframe - use "ch4" to represent the fourth sector
      metrics_data <- rbind(metrics_data, data.frame(
        Sensor = sensor,
        Channel = "ch4",  # Fourth sector
        Metric = "channel_acc",
        Mean = acc,
        CI_lower = acc ,  # Approximation for CI
        CI_upper = acc  # Approximation for CI, capped at 1.0
      ))
      
      metrics_data <- rbind(metrics_data, data.frame(
        Sensor = sensor,
        Channel = "ch4",  # Fourth sector
        Metric = "channel_micro_f1",
        Mean = micro_f1,
        CI_lower = micro_f1,  # Approximation for CI
        CI_upper = micro_f1  # Approximation for CI, capped at 1.0
      ))
      
      metrics_data <- rbind(metrics_data, data.frame(
        Sensor = sensor,
        Channel = "ch4",  # Fourth sector
        Metric = "channel_macro_f1",
        Mean = macro_f1,
        CI_lower = macro_f1,  # Approximation for CI
        CI_upper = macro_f1   # Approximation for CI, capped at 1.0
      ))
    }
  }
  
  return(metrics_data)
}

metric_display_labels <- c(
  "p_wave_amp" = "P-Wave Amplitude",
  "peak_to_peak" = "Peak-to-Peak",
  "t_wave_amp" = "T-Wave Amplitude",
  "NRMSE" = "Signal Quality (1-NRMSE)",
  "Cross-Correlation" = "Cross-Correlation",
  "participant_acc" = "Participant ID Accuracy",
  "participant_micro_f1" = "Participant ID Micro-F1",
  "participant_macro_f1" = "Participant ID Macro-F1",
  "channel_acc" = "Channel ID Accuracy",
  "channel_micro_f1" = "Channel ID Micro-F1",
  "channel_macro_f1" = "Channel ID Macro-F1"
)

# Extract and combine all metrics
participant_id_metrics <- extract_participant_id_metrics()
channel_id_metrics <- extract_channel_id_metrics()

# Combine with existing summaries
additional_metrics <- rbind(participant_id_metrics, channel_id_metrics)
additional_metrics$ChannelLabel <- channel_labels[additional_metrics$Channel]

additional_metrics$MetricLabel <- metric_display_labels[additional_metrics$Metric]

# Set SensorLabel based on Sensor
additional_metrics$SensorLabel <- ifelse(additional_metrics$Sensor == "AgCl", 
                                        "AgCl", "PPHG")

if(!"SD" %in% colnames(additional_metrics)) {
  additional_metrics$SD <- (additional_metrics$CI_upper - additional_metrics$CI_lower) / 3.92  # Approximation
}

if(!"N" %in% colnames(additional_metrics)) {
  additional_metrics$N <- 30  # Placeholder value
}

if(!"SE" %in% colnames(additional_metrics)) {
  additional_metrics$SE <- additional_metrics$SD / sqrt(additional_metrics$N)
}

# Convert to tibble for consistency
additional_metrics <- as_tibble(additional_metrics)


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

all_summaries <- bind_rows(all_summaries, additional_metrics)

channels <- c("ch1", "ch2", "ch3", "ch4")
channel_labels <- c("ch1" = "Lead 1", "ch2" = "Lead 2", "ch3" = "Lead 3", "ch4" = "Channel Identification")

normalize_data <- function(data) {
  # For each metric, apply appropriate normalization
  normalized_data <- data %>%
    group_by(Metric) %>%
    mutate(
      # For NRMSE, convert to 1-NRMSE so higher is better (like other metrics)
      # When inverting, the upper CI becomes the lower CI and vice versa
      Mean = ifelse(Metric == "NRMSE", 1-Mean, Mean),
      
      # Store original confidence intervals before swapping
      original_CI_lower = CI_lower,
      original_CI_upper = CI_upper,
      
      # For NRMSE, correctly invert confidence intervals
      # When inverting a metric, the upper bound becomes 1 minus the lower bound
      # and the lower bound becomes 1 minus the upper bound
      CI_lower = ifelse(Metric == "NRMSE", 1-original_CI_upper, original_CI_lower),
      CI_upper = ifelse(Metric == "NRMSE", 1-original_CI_lower, original_CI_upper)
    ) %>%
    select(-original_CI_lower, -original_CI_upper) %>%
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

all_summaries_norm <- normalize_data(all_summaries)


print("Unique Metrics in all_summaries_norm:")
print(unique(all_summaries_norm$Metric))


create_circular_barplot <- function(sensor_data, sensor_name, palette = NULL, show_legend = TRUE) {
  # Define the number of metrics and channels
  # Convert to character first to avoid factor level issues
  sensor_data$Metric <- as.character(sensor_data$Metric)
  unique_metrics <- unique(sensor_data$Metric)
  n_metrics <- length(unique_metrics)
  n_channels <- length(unique(sensor_data$Channel))
  
  # Debug print
  cat("Number of unique metrics in", sensor_name, "data:", n_metrics, "\n")
  cat("Unique metrics in", sensor_name, "data:", paste(unique_metrics, collapse=", "), "\n")
  
  metric_order <- c("p_wave_amp", "peak_to_peak", "t_wave_amp", "NRMSE", "Cross-Correlation",
                "participant_acc", "participant_micro_f1", "participant_macro_f1",
                "channel_acc", "channel_micro_f1", "channel_macro_f1")

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
  y_limit <- 1.05  # Set to 1.1 to leave some space for labels
  
  # Cap the confidence intervals at 100% (1.0) after normalization
  sensor_data <- sensor_data %>%
    mutate(
      CI_upper_norm = pmin(CI_upper_norm, 1.0),  # Cap upper CI at 1.0 (100%)
      CI_lower_norm = pmax(CI_lower_norm, 0.0)   # Ensure lower CI doesn't go below 0
    )

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
  
  if ("ch4" %in% unique(sensor_data$Channel)) {
    # Count bars per channel to determine proportional sector sizes
    bars_per_channel <- table(sensor_data$Channel)
    total_bars <- sum(bars_per_channel)
    
    # Calculate proportional angles
    channel_angles <- c()
    current_angle <- 0
    
    for (ch in names(bars_per_channel)) {
        # Add start angle of this sector
        channel_angles <- c(channel_angles, current_angle)
        # Calculate sector size proportionally to its bar count
        sector_size <- (bars_per_channel[ch] / total_bars) * 360
        current_angle <- current_angle + sector_size
    }
    
    # Add final boundary (360 degrees)
    channel_boundaries <- c(channel_angles, 360)
  } else {
    # If ch4 is not present, use equal distribution
    channel_boundaries <- seq(0, 360, by = 360/n_channels)
  }
  
  # Create channel sectors data for highlighting - MOVED THIS UP BEFORE USING IT
  channel_sectors <- data.frame()
  for (i in 1:n_channels) {
    start_angle <- channel_boundaries[i]
    end_angle <- ifelse(i == n_channels, 360, channel_boundaries[i+1])
    
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
    "NRMSE" = "Similarity (1-NRMSE)",
    "Cross-Correlation" = "Cross-Correlation",
    "participant_acc" = "Participant ID Accuracy",
    "participant_micro_f1" = "Participant ID Micro-F1",
    "participant_macro_f1" = "Participant ID Macro-F1",
    "channel_acc" = "Channel ID Accuracy",
    "channel_micro_f1" = "Channel ID Micro-F1",
    "channel_macro_f1" = "Channel ID Macro-F1"
  )
  
  # Only use labels for metrics we actually have
  metric_labels_to_use <- metric_display_labels[metric_order]
  
  # Create the circular barplot
  p <- ggplot() +
    # First add channel sectors as background
    geom_rect(data = channel_sectors, 
              aes(xmin = start, xmax = end, ymin = 0, ymax = y_limit),
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
    scale_y_continuous(
      limits = c(0, y_limit), 
      breaks = seq(0, 1, by = 0.25),
      labels = NULL  # Remove default labels as we'll add custom ones
    ) +
    scale_fill_manual(values = metric_colors, 
                     name = "Metric",
                     labels = metric_labels_to_use) +
    labs(
      title = paste0(sensor_name, " Performance Overview"),
      x = NULL,
      y = NULL
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),  # Remove default y-axis labels
      axis.ticks.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "gray80", size = 0.5),
      plot.title = element_text(size = 50, hjust = 0.5),
      plot.caption = element_text(size = 10, hjust = 0.5),
      legend.position = ifelse(show_legend, "bottom", "none"),
      legend.title = element_blank(),
      legend.text = element_text(size = 20),
      panel.border = element_blank(),
      axis.line = element_blank()
    )
  
  # Add custom y-axis tick labels - JUST ONCE instead of repeating
  # Add labels only at fixed positions (right side of the plot)
  tick_levels <- seq(0, 1, by = 0.25)
  tick_labels <- scales::percent_format(accuracy = 1)(tick_levels)
  
  # Add circular grid lines
  for (level in tick_levels) {
    if (level > 0) {  # Skip the innermost circle (0%)
      p <- p + geom_hline(
        yintercept = level,
        color = "gray80",
        size = 0.5,
        alpha = 0.8
      )
    }
  }
  
  # Add labels only on the right side of the circle (at 0 degrees)
  for (i in 1:length(tick_levels)) {
    p <- p + annotate(
      "text",
      x = 0,  # Right side of the circle
      y = tick_levels[i],
      label = tick_labels[i],
      size = 5,
      color = "gray30",
      hjust = -0.2  # Position text to the right of the point
    )
  }

  # Add channel boundary lines
  for (angle in channel_boundaries) {
    p <- p + geom_vline(xintercept = angle, 
                        linetype = "dashed", 
                        color = "gray50", 
                        size = 1)
  }
    
  return(p)
}

agcl_data <- all_summaries_norm %>% 
  filter((Sensor == "AgCl" & Metric %in% c("p_wave_amp", "peak_to_peak", "t_wave_amp", 
                                         "participant_acc", "participant_micro_f1", "participant_macro_f1", 
                                         "channel_acc", "channel_micro_f1", "channel_macro_f1")) | 
         (Metric %in% c("NRMSE", "Cross-Correlation")))

pphg_data <- all_summaries_norm %>% 
  filter((Sensor == "HG" & Metric %in% c("p_wave_amp", "peak_to_peak", "t_wave_amp", 
                                       "participant_acc", "participant_micro_f1", "participant_macro_f1", 
                                       "channel_acc", "channel_micro_f1", "channel_macro_f1")) | 
         (Metric %in% c("NRMSE", "Cross-Correlation")))



# Print the unique metrics for debugging
print("Unique metrics in agcl_data:")
print(unique(agcl_data$Metric))
print("Unique metrics in pphg_data:")
print(unique(pphg_data$Metric))

# Use droplevels to remove any unused factor levels
agcl_data$Metric <- droplevels(as.factor(agcl_data$Metric))
pphg_data$Metric <- droplevels(as.factor(pphg_data$Metric))

agcl_data <- agcl_data %>%
  mutate(
    CI_upper_norm = pmin(CI_upper_norm, 1.0),  # Cap upper CI at 1.0 (100%)
    CI_lower_norm = pmax(CI_lower_norm, 0.0)   # Ensure lower CI doesn't go below 0
  )

pphg_data <- pphg_data %>%
  mutate(
    CI_upper_norm = pmin(CI_upper_norm, 1.0),  # Cap upper CI at 1.0 (100%)
    CI_lower_norm = pmax(CI_lower_norm, 0.0)   # Ensure lower CI doesn't go below 0
  )

# Update metric order to ensure we only have the 5 metrics
metric_order <- c("p_wave_amp", "peak_to_peak", "t_wave_amp", "NRMSE", "Cross-Correlation")

# Create a 5-color palette for the 5 metrics
metric_palette <- c(
  "#E41A1C",  # red - p_wave_amp
  "#377EB8",  # light blue - peak_to_peak  
  "#4DAF4A",  # green - t_wave_amp
  "#984EA3",  # purple - NRMSE
  "#FF7F00",  # orange - Cross-Correlation
  "#A65628",  # brown - participant_acc
  "#F781BF",  # pink - participant_micro_f1
  "#999999",  # grey - participant_macro_f1
  "#6933ff",  # blue - channel_acc
  "#A6761D",  # darker brown - channel_micro_f1
  "#1B9E77"   # teal - channel_macro_f1
)
agcl_circular <- create_circular_barplot(agcl_data, "AgCl", metric_palette, show_legend = TRUE)
pphg_circular <- create_circular_barplot(pphg_data, "PPHG", metric_palette, show_legend = TRUE)


agcl_with_labels <- ggdraw(agcl_circular) + 
  draw_label("Lead 1", x = 0.8, y = 0.8,  size = 35, 
             hjust = 0.5, vjust = 0.5, color = "black", 
             fontfamily = "sans") +
  draw_label("Lead 2", x = 0.6, y = 0.12, size = 35, 
             hjust = 0.5, vjust = 0.5, color = "black", 
             fontfamily = "sans") +
  draw_label("Lead 3", x = 0.1, y = 0.5,  size = 35, 
             hjust = 0.5, vjust = 0.5, color = "black", 
             fontfamily = "sans") +
  draw_label("Channel ID", x = 0.3, y = 0.9, size = 35, 
             hjust = 0.5, vjust = 0.5, color = "black", 
             fontfamily = "sans")

pphg_with_labels <- ggdraw(pphg_circular) + 
  draw_label("Lead 1", x = 0.8, y = 0.8, size = 35, 
             hjust = 0.5, vjust = 0.5, color = "black", 
             fontfamily = "sans") +
  draw_label("Lead 2", x = 0.6, y = 0.12, size = 35, 
             hjust = 0.5, vjust = 0.5, color = "black", 
             fontfamily = "sans") +
  draw_label("Lead 3", x = 0.1, y = 0.5, size = 35, 
             hjust = 0.5, vjust = 0.5, color = "black", 
             fontfamily = "sans") +
  draw_label("Channel ID", x = 0.3, y = 0.9, size = 35, 
             hjust = 0.5, vjust = 0.5, color = "black", 
             fontfamily = "sans")

ggsave(paste0(output_dir, "agcl_circular_barplot.png"), agcl_with_labels, 
       width = 16, height = 12, dpi = 300, bg = "white")
ggsave(paste0(output_dir, "pphg_circular_barplot.png"), pphg_with_labels, 
       width = 16, height = 12, dpi = 300, bg = "white")

# Extract the legend from the AgCl plot to use in the combined plot
legend <- cowplot::get_legend(
  agcl_circular + 
    theme(legend.position = "bottom",
      legend.box.margin = margin(t = 25),    # Increased top margin
      legend.key.size = unit(1.5, "cm"),     # Increased key size
      legend.title = element_blank(), # Larger title
      legend.text = element_text(size = 19)   # Larger text
    )
)

# Create plots without legends for the combined figure
agcl_circular_no_legend <- create_circular_barplot(agcl_data, "AgCl", metric_palette, show_legend = FALSE)
pphg_circular_no_legend <- create_circular_barplot(pphg_data, "PPHG", metric_palette, show_legend = FALSE)

# Combine the plots with a shared legend at the bottom
combined_plot <- cowplot::plot_grid(
  cowplot::plot_grid(
    agcl_circular_no_legend, 
    pphg_circular_no_legend, 
    labels = c("A", "B"), 
    ncol = 2, 
    align = "h"
  ),
  legend,
  ncol = 1,
  rel_heights = c(0.85, 0.15)  # Adjust these values to control the relative size
)

# Save the combined plot
#ggsave(paste0(output_dir, "combined_circular_barplots.png"), combined_plot, 
#       width = 18, height = 12, dpi = 300, bg = "white")