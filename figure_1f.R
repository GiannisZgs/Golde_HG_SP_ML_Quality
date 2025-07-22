library(jsonlite)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(R.matlab)
library(scales)
library(grid)
library(cowplot)
library(ggrepel)
library(viridis)
library(ggbeeswarm)  

"Plot for the whole dataset instead of cleaned to have more samples"

# Load data
fs <- 200
data_path <- paste0("/home/giannis/Documents/ECG HG paper/results_data/metrics_deviation_from_noise.json")
data <- fromJSON(txt = data_path, simplifyVector = FALSE)

# Create output directory
output_dir <- "R_figures/figure_1f/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define metric names and their display labels
metrics <- c(
  "high_freq_deviation", 
  "low_freq_deviation", 
  "mid_freq_deviation",
  "powerline_deviation", 
  "powerline_harmonic_deviation", 
  "trend_deviation",
  "ultra_high1_freq_deviation",
  "ultra_high2_freq_deviation",
  "vlf_deviation"
)

metric_labels <- c(
  "High Frequency (15-40 Hz)",
  "Low Frequency (0.5-5 Hz)",
  "Mid Frequency (5-15 Hz)",
  "Powerline Noise (50 Hz)",
  "Powerline Noise Harmonic (80 Hz)",
  "Baseline Wander (0-0.3 Hz)",
  "Ultra High Frequency 1 (52-78 Hz)",
  "Ultra High Frequency 2 (82-100 Hz)",
  "Very Low Frequency (0.05-0.5 Hz)"
)

names(metric_labels) <- metrics

# Define channels and sensors
channels <- c("ch1", "ch2", "ch3")
sensors <- c("AgCl", "HG1", "HG2")
sensor_labels <- c("AgCl" = "AgCl", "HG1" = "PPHG", "HG2" = "PPHG")
names(sensor_labels) <- sensors

# Function to extract data for all metrics across participants
extract_metrics_data <- function() {
  all_data <- data.frame()
  
  # Loop through all participants
  for (participant in names(data$metrics)) {
    if (participant == "_fieldnames") {
      next  # Skip the metadata field
    }
    participant_data <- data$metrics[[participant]]
    
    # Loop through sensors
    for (sensor in sensors) {
      if (sensor %in% names(participant_data)) {
        sensor_data <- participant_data[[sensor]]
        
        # Loop through metrics
        for (metric in metrics) {
          if (metric %in% names(sensor_data)) {
            metric_data <- sensor_data[[metric]]
            
            # Loop through channels
            for (channel in channels) {
              if (channel %in% names(metric_data)) {
                value <- metric_data[[channel]]
                
                # Combine HG1 and HG2 into a single PPHG sensor
                sensor_label <- ifelse(sensor %in% c("HG1", "HG2"), "PPHG", sensor_labels[sensor])
                
                # Add to dataframe
                all_data <- rbind(all_data, data.frame(
                  Participant = participant,
                  Sensor = sensor_label,  # Use combined sensor label
                  SensorLabel = sensor_label,
                  Metric = metric,
                  MetricLabel = metric_labels[metric],
                  Channel = channel,
                  Value = as.numeric(value)
                ))
              }
            }
          }
        }
      }
    }
  }
  
  return(all_data)
}

# Extract all metrics data
metrics_data <- extract_metrics_data()

# Clean up channel names for display
metrics_data$Channel <- gsub("ch1", "Lead 1", metrics_data$Channel)
metrics_data$Channel <- gsub("ch2", "Lead 2", metrics_data$Channel)
metrics_data$Channel <- gsub("ch3", "Lead 3", metrics_data$Channel)

metrics_data$SensorChannel <- paste(metrics_data$SensorLabel, metrics_data$Channel)


print("Unique sensor-channel combinations:")
print(unique(metrics_data$SensorChannel))

# Define colors based on the actual combinations
sensor_colors <- c("AgCl" = "#d9020d", "PPHG" = "#0066CC")

# Create a new vector to map colors to the actual SensorChannel values
channel_colors <- vector()

# Map colors to each unique SensorChannel value
unique_combos <- unique(metrics_data$SensorChannel)
for(combo in unique_combos) {
  if(grepl("^AgCl", combo)) {
    # Red colors for AgCl with different shades per channel
    if(grepl("ch1$", combo)) {
      channel_colors[combo] <- "#d9020d"  # Full red for AgCl ch1
    } else if(grepl("ch2$", combo)) {
      channel_colors[combo] <- "#e33e3e"  # Slightly lighter red for AgCl ch2
    } else {
      channel_colors[combo] <- "#eb7070"  # Even lighter red for AgCl ch3
    }
  } else {
    # Blue colors for PPHG with different shades per channel
    if(grepl("ch1$", combo)) {
      channel_colors[combo] <- "#0066CC"  # Full blue for PPHG ch1
    } else if(grepl("ch2$", combo)) {
      channel_colors[combo] <- "#3385d6"  # Slightly lighter blue for PPHG ch2
    } else {
      channel_colors[combo] <- "#66a3e0"  # Even lighter blue for PPHG ch3
    }
  }
}

# Convert SensorChannel to factor after colors are defined
metrics_data$SensorChannel <- factor(metrics_data$SensorChannel)
print("Sensor channel levels after factoring:")
print(levels(metrics_data$SensorChannel))

# Print the color mapping to confirm it matches the factor levels
print("Final color mapping:")
print(channel_colors)

# Function to create a single metric plot (for the combined figure)
create_metric_plot <- function(data, metric_name, plot_type = "boxplot") {
  metric_data <- data[data$Metric == metric_name,]
  
  # Set up the base plot
  p <- ggplot(metric_data, aes(x = SensorChannel, y = Value))
  
  # Add appropriate geoms based on plot type
  if (plot_type == "boxplot") {
    p <- p + 
      geom_boxplot(aes(fill = SensorChannel), alpha = 0.7, outlier.shape = NA) +
      geom_jitter(width = 0.2, height = 0, alpha = 0.6, size = 1.5) +
      scale_fill_manual(values = channel_colors)
  } else if (plot_type == "violin") {
    p <- p + 
      geom_violin(aes(fill = SensorChannel), alpha = 0.7, scale = "width", trim = TRUE) +
      geom_jitter(width = 0.1, height = 0, alpha = 0.6, size = 1) +
      scale_fill_manual(values = channel_colors)
  } else if (plot_type == "beeswarm") {
    p <- p + 
      geom_beeswarm(aes(color = SensorChannel), size = 1.5, alpha = 0.7, cex = 1.5) +
      scale_color_manual(values = channel_colors)
  }
  
  p <- p +
    labs(
      title = metric_labels[metric_name],
      x = "",
      y = "Deviation"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.title = element_text(size = 10),
      plot.title = element_text(size = 10, face = "bold"),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  return(p)
}


# Create and save combined plots
print("Creating combined boxplots...")
metric_plots <- list()
for (metric in metrics) {
  metric_plots[[metric]] <- create_metric_plot(metrics_data, metric, "boxplot")
}

combined_boxplots <- wrap_plots(metric_plots, ncol = 3) +
  plot_annotation(
    title = "",
    subtitle = "",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 14),
      plot.background = element_rect(fill = "white", color = NA)
    )
  ) &
  theme(panel.background = element_rect(fill = "white", color = NA))

ggsave(paste0(output_dir, "combined_metrics_boxplot.png"), combined_boxplots, 
       width = 15, height = 12, dpi = 300, bg = "white")
print("Saved combined boxplots")

print("Creating combined violin plots...")
metric_violin_plots <- list()
for (metric in metrics) {
  metric_violin_plots[[metric]] <- create_metric_plot(metrics_data, metric, "violin")
}

combined_violins <- wrap_plots(metric_violin_plots, ncol = 3) +
  plot_annotation(
    title = "",
    subtitle = "",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 14),
      plot.background = element_rect(fill = "white", color = NA)
    )
  ) &
  theme(panel.background = element_rect(fill = "white", color = NA))

ggsave(paste0(output_dir, "combined_metrics_violin.png"), combined_violins, 
       width = 15, height = 12, dpi = 300, bg = "white")
print("Saved combined violin plots")

# Create and save combined beeswarm plots
print("Creating combined beeswarm plots...")
metric_beeswarm_plots <- list()
for (metric in metrics) {
  metric_beeswarm_plots[[metric]] <- create_metric_plot(metrics_data, metric, "beeswarm")
}

combined_beeswarms <- wrap_plots(metric_beeswarm_plots, ncol = 3) +
  plot_annotation(
    title = "",
    subtitle = "",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 14),
      plot.background = element_rect(fill = "white", color = NA)
    )
  ) &
  theme(panel.background = element_rect(fill = "white", color = NA))

ggsave(paste0(output_dir, "combined_metrics_beeswarm.png"), combined_beeswarms, 
       width = 15, height = 12, dpi = 300, bg = "white")
print("Saved combined beeswarm plots")



# Save individual metric plots
print("Creating and saving individual metric plots...")

# Create a subdirectory for individual metric plots
individual_dir <- paste0(output_dir, "individual_metrics/")
if (!dir.exists(individual_dir)) {
  dir.create(individual_dir, recursive = TRUE)
}

custom_y_breaks <- function(limits) {
  breaks <- pretty(limits)  # Generate default breaks
  if (1 >= limits[1] && 1 <= limits[2]) {
    breaks <- sort(unique(c(breaks, 1)))  # Ensure 1 is included
  }
  return(breaks)
}

# Save individual boxplots for each metric
for (metric in metrics) {
  # Create individual boxplot for this metric
  metric_data <- metrics_data[metrics_data$Metric == metric, ]
  
  # Boxplot
  p_box <- ggplot(metric_data, aes(x = Channel, y = Value, fill = SensorLabel)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, height = 0, alpha = 0.6, size = 2) +
    scale_fill_manual(values = sensor_colors) +
    scale_y_continuous(breaks = function(limits) custom_y_breaks(limits)) +  # Custom y-axis breaks
    labs(
      title = paste(metric_labels[metric]),
      x = "",
      y = expression(frac("Processed Power", "Raw Power")),
      fill = "Sensor"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 25),
      axis.text.y = element_text(size = 25),
      axis.title = element_text(size = 25),
      legend.position = "none",
      legend.title = element_text(size = 12),
      plot.title = element_text(size = 25),
      strip.text = element_text(size = 25),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    facet_grid(~ SensorLabel, scales = "free_x", space = "free_x")
  
  # Save the boxplot
  ggsave(paste0(individual_dir, gsub(" ", "_", metric_labels[metric]), "_boxplot.png"), 
         p_box, width = 8, height = 6, dpi = 300, bg = "white")
  
  # Violin plot
  p_violin <- ggplot(metric_data, aes(x = Channel, y = Value, fill = SensorLabel)) +
    geom_violin(alpha = 0.7, scale = "width", trim = TRUE) +
    geom_jitter(width = 0.1, height = 0, alpha = 0.6, size = 1.5) +
    scale_fill_manual(values = sensor_colors) +
    scale_y_continuous(breaks = function(limits) custom_y_breaks(limits)) +  # Custom y-axis breaks
    labs(
      title = paste(metric_labels[metric]),
      x = "",
      y = expression(frac("Processed Power", "Raw Power")),
      fill = "Sensor"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
      axis.text.y = element_text(size = 18),
      axis.title = element_text(size = 18),
      legend.position = "none",
      legend.title = element_text(size = 12),
      plot.title = element_text(size = 20),
      strip.text = element_text(size = 18),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    facet_grid(~ SensorLabel, scales = "free_x", space = "free_x")
  
  # Save the violin plot
  ggsave(paste0(individual_dir, gsub(" ", "_", metric_labels[metric]), "_violin.png"), 
         p_violin, width = 8, height = 6, dpi = 300, bg = "white")
  
  # Beeswarm plot
  p_beeswarm <- ggplot(metric_data, aes(x = Channel, y = Value, color = SensorLabel)) +
    geom_beeswarm(size = 2.5, alpha = 0.7, cex = 2) +
    scale_color_manual(values = sensor_colors) +
    scale_y_continuous(breaks = function(limits) custom_y_breaks(limits)) +  # Custom y-axis breaks
    labs(
      title = paste(metric_labels[metric]),
      x = "",
      y = expression(frac("Processed Power", "Raw Power")),
      fill = "Sensor"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
      axis.text.y = element_text(size = 18),
      axis.title = element_text(size = 18),
      legend.position = "none",
      legend.title = element_text(size = 12),
      plot.title = element_text(size = 20),
      strip.text = element_text(size = 18),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    facet_grid(~ SensorLabel, scales = "free_x", space = "free_x")
  
  # Save the beeswarm plot
  ggsave(paste0(individual_dir, gsub(" ", "_", metric_labels[metric]), "_beeswarm.png"), 
         p_beeswarm, width = 8, height = 6, dpi = 300, bg = "white")
  
  print(paste("Saved individual plots for metric:", metric_labels[metric]))
}

print("Completed saving all individual metric plots")