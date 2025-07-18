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

nans_removed = TRUE
data_path <- "/home/giannis/Documents/ECG HG paper/results_data/heartbeat_profiles_MA.json"

if (nans_removed) {
  data <- fromJSON(txt = data_path, simplifyVector = FALSE)

} else {
  
  json_text <- readLines(data_path)

  # Replace all occurrences of NaN with null
  json_text <- gsub("\\bNaN\\b", "null", json_text)

  # Parse the modified JSON text
  data <- fromJSON(txt = paste(json_text, collapse = "\n"), simplifyVector = FALSE)

  cleaned_json <- toJSON(data, pretty = TRUE, auto_unbox = TRUE)

  # Write the JSON to a file
  write(cleaned_json, file = data_path)

  print(paste("Cleaned data saved to:", data_path))
}

fs <- 200
all_participants <- names(data$profiling_struct[-1])  # Assuming AgCl has all participants
channels <- c("ch1", "ch2", "ch3")
sensors <- c("AgCl", "HG")

# Create output directory if it doesn't exist
output_dir <- paste0("R_figures/figure_2a/")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Initialize an empty data frame for waveforms
waveform_data <- data.frame()

# Extract waveforms for all participants, sensors, and channels
for (participant in all_participants) {
  for (sensor in sensors) {
    for (channel in channels) {
      # Check if the sensor and channel exist for the participant
      if (!is.null(data$profiling_struct[[participant]][[sensor]][[channel]])) {
        channel_data <- data$profiling_struct[[participant]][[sensor]][[channel]]
        
        # Combine HG1 and HG2 as PPHG
        sensor_label <- ifelse(sensor == "HG", "PPHG", sensor)
        
        if (length(channel_data) == 121) {
          # Use lapply to process all time steps at once
          temp_df <- do.call(rbind, lapply(seq_along(channel_data), function(time_step) {
            amplitudes <- unlist(channel_data[[time_step]])
            data.frame(
              Time = time_step,
              Amplitude = amplitudes,
              Participant = participant,
              Channel = channel,
              Sensor = sensor_label
            )
          }))
            

          # Append to the main data frame
          waveform_data <- rbind(waveform_data, temp_df)
        }
      }
    }
  }
}

# Remove rows with NaN values
waveform_data <- waveform_data %>% filter(!is.na(Amplitude))

# Normalize time to seconds
waveform_data$Time <- (waveform_data$Time - 1) / fs

# Map channel names to "Lead 1", "Lead 2", "Lead 3"
waveform_data$Channel <- gsub("ch1", "Lead 1", waveform_data$Channel)
waveform_data$Channel <- gsub("ch2", "Lead 2", waveform_data$Channel)
waveform_data$Channel <- gsub("ch3", "Lead 3", waveform_data$Channel)

# Combine Sensor and Channel into a single variable
waveform_data$SensorChannel <- paste(waveform_data$Sensor, waveform_data$Channel)

# Create a list of unique SensorChannel combinations
sensor_channels <- unique(waveform_data$SensorChannel)


amp_counts <- waveform_data %>%
  group_by(Participant, Sensor, Channel, Time) %>%
  summarize(count = n(), .groups = 'drop')
print(head(amp_counts))

waveform_data_agg <- waveform_data %>%
  group_by(Time, Participant, Channel, Sensor, SensorChannel) %>%
  summarize(Amplitude = mean(Amplitude, na.rm = TRUE), .groups = 'drop')

# Initialize an empty list to store plots
plot_list <- list()



aggregated_dir <- paste0(output_dir, "aggregated/")
raw_dir <- paste0(output_dir, "raw/")

# Create directories if they don't exist
if (!dir.exists(aggregated_dir)) dir.create(aggregated_dir, recursive = TRUE)
if (!dir.exists(raw_dir)) dir.create(raw_dir, recursive = TRUE)

# OPTION 1: Aggregated waveforms (one line per participant) for clearer visualization
plot_list_agg <- list()

# Generate overlay plots for each SensorChannel using aggregated data
for (sensor_channel in unique(waveform_data_agg$SensorChannel)) {
  # Subset data for the current SensorChannel
  subset_data <- waveform_data_agg[waveform_data_agg$SensorChannel == sensor_channel, ]
  
  # Create the overlay plot
  p <- ggplot(subset_data, aes(x = Time, y = Amplitude, group = Participant)) +
    geom_line(alpha = 0.3, size = 0.5, color = ifelse(grepl("AgCl", sensor_channel), "#d9020d", "#0066CC")) +
    labs(
      title = sensor_channel,
      x = "Time (s)",
      y = "Amplitude"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 20),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 18),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95")
    )
  
  # Add the plot to the list
  plot_list_agg[[sensor_channel]] <- p
}

# Combine all aggregated plots into a single figure using patchwork
if (length(plot_list_agg) > 0) {
  agcl_plots <- plot_list_agg[grep("AgCl", names(plot_list_agg))]
  pphg_plots <- plot_list_agg[grep("PPHG", names(plot_list_agg))]

  # Sort each group by lead number
  agcl_plots <- agcl_plots[order(names(agcl_plots))]
  pphg_plots <- pphg_plots[order(names(pphg_plots))]

  # Combine into a single ordered list
  ordered_plots <- c(agcl_plots, pphg_plots)

  # Create the combined plot with this specific order
  combined_plot_agg <- wrap_plots(ordered_plots, ncol = 3, nrow = 2) + 
    plot_layout(guides = "collect") & 
    theme(legend.position = "none")
  
  # Display the combined plot
  print(combined_plot_agg)
  
  # Save the combined plot
  ggsave(paste0(aggregated_dir, "overlay_waveforms_aggregated.png"), combined_plot_agg, 
         width = 15, height = 10, dpi = 300, bg = "white")
}

# OPTION 2: Raw data with all points to show variability
plot_list_raw <- list()

# Generate overlay plots for each SensorChannel using raw data
for (sensor_channel in unique(waveform_data$SensorChannel)) {
  # Subset data for the current SensorChannel
  subset_data <- waveform_data[waveform_data$SensorChannel == sensor_channel, ]
  
  # Create the overlay plot with raw data
  p <- ggplot(subset_data, aes(x = Time, y = Amplitude)) +
    geom_point(alpha = 0.1, size = 0.3, color = ifelse(grepl("AgCl", sensor_channel), "#d9020d", "#0066CC")) +
    labs(
      title = paste0(sensor_channel),
      x = "Time (s)",
      y = "Amplitude"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 20),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 18),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95")
    )
  
  # Add the plot to the list
  plot_list_raw[[sensor_channel]] <- p
}

# Combine all raw plots into a single figure using patchwork
if (length(plot_list_raw) > 0) {
  # Do the same for raw plots
  agcl_plots_raw <- plot_list_raw[grep("AgCl", names(plot_list_raw))]
  pphg_plots_raw <- plot_list_raw[grep("PPHG", names(plot_list_raw))]

  # Sort each group by lead number
  agcl_plots_raw <- agcl_plots_raw[order(names(agcl_plots_raw))]
  pphg_plots_raw <- pphg_plots_raw[order(names(pphg_plots_raw))]

  # Combine into a single ordered list
  ordered_plots_raw <- c(agcl_plots_raw, pphg_plots_raw)

  # Create the combined plot with this specific order
  combined_plot_raw <- wrap_plots(ordered_plots_raw, ncol = 3, nrow = 2) + 
    plot_layout(guides = "collect") & 
    theme(legend.position = "none")
  
  # Display the combined plot
  print(combined_plot_raw)
  
  # Save the combined plot
  ggsave(paste0(raw_dir, "overlay_waveforms_raw.png"), combined_plot_raw, 
         width = 15, height = 10, dpi = 300, bg = "white")
}


"Example for single participant"

single_participant_dir <- paste0(output_dir, "single_participant/")
if (!dir.exists(single_participant_dir)) dir.create(single_participant_dir, recursive = TRUE)


# You can change this to any specific participant ID you want
single_participant <- "p39"

# Create plots for a single participant
plot_list_single <- list()

# Generate overlay plots for each SensorChannel for the selected participant
for (sensor_channel in unique(waveform_data$SensorChannel)) {
  # Subset data for the current SensorChannel and participant
  subset_data <- waveform_data[waveform_data$SensorChannel == sensor_channel & 
                               waveform_data$Participant == single_participant, ]
  
  # Skip if no data exists for this sensor-channel for the selected participant
  if(nrow(subset_data) == 0) next
  
  # Create the plot with raw data for a single participant
  p <- ggplot(subset_data, aes(x = Time, y = Amplitude)) +
    geom_point(alpha = 0.3, size = 0.5, color = ifelse(grepl("AgCl", sensor_channel), "#d9020d", "#0066CC")) +
    labs(
      title = paste0(sensor_channel),
      x = "Time (s)",
      y = "Amplitude"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 20),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 18),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95")
    )
  
  # Add the plot to the list
  plot_list_single[[sensor_channel]] <- p
}

# Combine all single participant plots into a single figure using patchwork
if (length(plot_list_single) > 0) {
  # Sort by sensor type and lead number
  agcl_plots_single <- plot_list_single[grep("AgCl", names(plot_list_single))]
  pphg_plots_single <- plot_list_single[grep("PPHG", names(plot_list_single))]

  # Sort each group by lead number
  agcl_plots_single <- agcl_plots_single[order(names(agcl_plots_single))]
  pphg_plots_single <- pphg_plots_single[order(names(pphg_plots_single))]

  # Combine into a single ordered list
  ordered_plots_single <- c(agcl_plots_single, pphg_plots_single)

  # Create the combined plot with this specific order
  combined_plot_single <- wrap_plots(ordered_plots_single, ncol = 3, nrow = 2) + 
    plot_layout(guides = "collect") & 
    theme(legend.position = "none") &
    labs(title = paste0("Participant: ", single_participant, ""))
  
  # Display the combined plot
  print(combined_plot_single)
  
  # Save the combined plot
  ggsave(paste0(single_participant_dir, "single_participant_raw_", single_participant, ".png"), 
         combined_plot_single, width = 15, height = 10, dpi = 300, bg = "white")
  
  # Also save individual plots for this participant
  for (sensor_channel in names(plot_list_single)) {
    clean_name <- gsub(" ", "_", sensor_channel)
    ggsave(paste0(single_participant_dir, single_participant, "_", clean_name, ".png"), 
           plot_list_single[[sensor_channel]], width = 6, height = 4, dpi = 300, bg = "white")
  }
}

# Loop through all participants to create individual plots for each
generate_all_participant_plots <- FALSE

if (generate_all_participant_plots) {
  for (current_participant in all_participants) {
    participant_plot_list <- list()
    
    for (sensor_channel in unique(waveform_data$SensorChannel)) {
      subset_data <- waveform_data[waveform_data$SensorChannel == sensor_channel & 
                                   waveform_data$Participant == current_participant, ]
      
      # Skip if no data exists
      if(nrow(subset_data) == 0) next
      
      # Create the plot
      p <- ggplot(subset_data, aes(x = Time, y = Amplitude)) +
        geom_point(alpha = 0.3, size = 0.5, color = ifelse(grepl("AgCl", sensor_channel), "#d9020d", "#0066CC")) +
        labs(
          title = paste0(sensor_channel),
          x = "Time (s)",
          y = "Amplitude"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 16),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12)
        )
      
      participant_plot_list[[sensor_channel]] <- p
    }
    
    # Only create combined plot if we have data for this participant
    if (length(participant_plot_list) > 0) {
      # Sort and combine plots
      agcl_plots_p <- participant_plot_list[grep("AgCl", names(participant_plot_list))]
      pphg_plots_p <- participant_plot_list[grep("PPHG", names(participant_plot_list))]
      
      agcl_plots_p <- agcl_plots_p[order(names(agcl_plots_p))]
      pphg_plots_p <- pphg_plots_p[order(names(pphg_plots_p))]
      
      ordered_plots_p <- c(agcl_plots_p, pphg_plots_p)
      
      combined_plot_p <- wrap_plots(ordered_plots_p, ncol = 3, nrow = 2) + 
        plot_layout(guides = "collect") & theme(legend.position = "none")
      
      # Save the combined plot for this participant
      participant_dir <- paste0(single_participant_dir, current_participant, "/")
      if (!dir.exists(participant_dir)) dir.create(participant_dir, recursive = TRUE)
      
      ggsave(paste0(participant_dir, "raw_waveforms.png"), 
             combined_plot_p, width = 15, height = 10, dpi = 300, bg = "white")
    }
  }
}