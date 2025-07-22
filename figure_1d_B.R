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


#data <- fromJSON(txt = "/home/giannis/Documents/ECG HG paper/results_data/manually_cleaned_participant_id_results.json", simplifyVector = FALSE)

data <- fromJSON(txt = "/home/giannis/Documents/ECG HG paper/results_data/manually_cleaned_channel_id_results.json", simplifyVector = FALSE)

# Create output directory if it doesn't exist
output_dir <- paste0("R_figures/figure_1d_B/without_labels/")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

show_labels <- FALSE

fs <- 200
selected_channels <- c("ch1", "ch2","ch3")

waveform_data <- data.frame()
features_data <- data.frame()
for (channel in selected_channels) {
  # Process AgCl data
  if (channel %in% names(data$centroids$AgCl)) {
    agcl_waveform <- data$centroids$AgCl[[channel]]$kmeans_cluster_mapped_waveform
    
    # Extract features for AgCl if available
    if (!is.null(data$features)) {
      p_wave_amp <- data$features$AgCl[[channel]]$p_wave_amp
      peak_to_peak <- data$features$AgCl[[channel]]$peak_to_peak
      t_wave_amp <- data$features$AgCl[[channel]]$t_wave_amp
      
      # Add to features data frame
      features_data <- rbind(features_data, data.frame(
        Channel = channel,
        Sensor = "AgCl",
        p_wave_amp = ifelse(is.null(p_wave_amp), NA, as.numeric(p_wave_amp)),
        peak_to_peak = ifelse(is.null(peak_to_peak), NA, as.numeric(peak_to_peak)),
        t_wave_amp = ifelse(is.null(t_wave_amp), NA, as.numeric(t_wave_amp))
      ))
    }
    
    if (!is.null(agcl_waveform)) {
      temp_df <- data.frame(
        Time = 1:length(agcl_waveform),
        Amplitude = as.numeric(agcl_waveform),
        Channel = channel,
        Sensor = "AgCl"
      )
      waveform_data <- rbind(waveform_data, temp_df)
    }
  }
  
  # Process HG data
  if (channel %in% names(data$centroids$HG)) {
    hg_waveform <- data$centroids$HG[[channel]]$kmeans_cluster_mapped_waveform
    
    # Extract features for PPHG if available
    if (!is.null(data$features)) {
      p_wave_amp <- data$features$HG[[channel]]$p_wave_amp
      peak_to_peak <- data$features$HG[[channel]]$peak_to_peak
      t_wave_amp <- data$features$HG[[channel]]$t_wave_amp
      
      # Add to features data frame
      features_data <- rbind(features_data, data.frame(
        Channel = channel,
        Sensor = "PPHG",
        p_wave_amp = ifelse(is.null(p_wave_amp), NA, as.numeric(p_wave_amp)),
        peak_to_peak = ifelse(is.null(peak_to_peak), NA, as.numeric(peak_to_peak)),
        t_wave_amp = ifelse(is.null(t_wave_amp), NA, as.numeric(t_wave_amp))
      ))
    }
    
    if (!is.null(hg_waveform)) {
      # Create a data frame for this waveform
      temp_df <- data.frame(
        Time = 1:length(hg_waveform),
        Amplitude = as.numeric(hg_waveform),
        Channel = channel,
        Sensor = "PPHG"
      )
      waveform_data <- rbind(waveform_data, temp_df)
    }
  }
}

# Normalize time to seconds 
waveform_data$Time <- (waveform_data$Time - 1) / fs

# Create annotations for each channel and sensor
annotations <- data.frame()

for (channel in selected_channels) {
  for (sensor in c("AgCl", "PPHG")) {
    subset_data <- waveform_data[waveform_data$Channel == channel &
                               waveform_data$Sensor == sensor, ]
    
    if (nrow(subset_data) > 0) {
      # Estimate P wave position (assuming it's at 1/4 of the waveform)
      p_wave_time <- min(subset_data$Time) + 0.2 * (max(subset_data$Time) - min(subset_data$Time))
      p_wave_idx <- which.min(abs(subset_data$Time - p_wave_time))
      p_wave_amp <- subset_data$Amplitude[p_wave_idx]
      
      # Find R peak (assuming it's the max value)
      r_peak_idx <- which.max(subset_data$Amplitude)
      r_peak_time <- subset_data$Time[r_peak_idx]
      r_peak_amp <- subset_data$Amplitude[r_peak_idx]
      
      q_window_start <- max(1, r_peak_idx - round(0.1 * fs))
      if (q_window_start < r_peak_idx) {
        q_window <- subset_data$Amplitude[q_window_start:r_peak_idx]
        q_offset <- which.min(q_window) - 1
        q_peak_idx <- q_window_start + q_offset
        q_peak_time <- subset_data$Time[q_peak_idx]
        q_peak_amp <- subset_data$Amplitude[q_peak_idx]
      } else {
        # Fallback if window is problematic
        q_peak_time <- r_peak_time - 0.05 # Approximate 50ms before R peak
        q_peak_amp <- min(subset_data$Amplitude) # Use overall minimum as fallback
      }

      global_min_idx <- which.min(subset_data$Amplitude)
      global_min_time <- subset_data$Time[global_min_idx]
      global_min_amp <- subset_data$Amplitude[global_min_idx]

      # Estimate T wave position (assuming it's at 3/4 of the waveform)
      t_wave_time <- min(subset_data$Time) + 0.8 * (max(subset_data$Time) - min(subset_data$Time))
      t_wave_idx <- which.min(abs(subset_data$Time - t_wave_time))
      t_wave_amp <- subset_data$Amplitude[t_wave_idx]
      
      # Get features from the features data frame
      feature_row <- features_data[features_data$Channel == channel &
                                 features_data$Sensor == sensor, ]
      
      p_wave_feature <- if(nrow(feature_row) > 0 && !is.na(feature_row$p_wave_amp)) 
                          feature_row$p_wave_amp else NA
      peak_to_peak_feature <- if(nrow(feature_row) > 0 && !is.na(feature_row$peak_to_peak)) 
                               feature_row$peak_to_peak else NA
      t_wave_feature <- if(nrow(feature_row) > 0 && !is.na(feature_row$t_wave_amp)) 
                         feature_row$t_wave_amp else NA
      
      # Add to annotations data frame
      annotations <- rbind(annotations, data.frame(
        Channel = channel,
        Sensor = sensor,
        p_wave_time = p_wave_time,
        p_wave_amp = p_wave_amp,
        p_wave_feature = p_wave_feature,
        q_peak_time = q_peak_time,
        q_peak_amp = q_peak_amp,
        r_peak_time = r_peak_time,
        r_peak_amp = r_peak_amp,
        global_min_time = global_min_time,
        global_min_amp = global_min_amp,
        peak_to_peak_feature = peak_to_peak_feature,
        t_wave_time = t_wave_time,
        t_wave_amp = t_wave_amp,
        t_wave_feature = t_wave_feature
      ))
    }
  }
}

# Create plots for each channel
plot_list <- list()

for (channel in selected_channels) {
  # Get data for this channel
  subset_data <- waveform_data[waveform_data$Channel == channel, ]
  subset_annotations <- annotations[annotations$Channel == channel, ]
  
  if (nrow(subset_data) > 0) {
    # Create a base plot without shaded areas first
    p <- ggplot(subset_data, aes(x = Time, y = Amplitude, color = Sensor)) +
      # Use different colors for AgCl and PPHG
      scale_color_manual(values = c("AgCl" = "#d9020d", "PPHG" = "#0066CC")) +
      scale_fill_manual(values = c("AgCl" = "#d9020d", "PPHG" = "#0066CC")) +
      labs(
        title = "",
        x = "Time (s)",
        y = "Amplitude",
        color = "",
        fill = "Sensor Type"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 25),
        legend.position = "bottom",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 25),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_line(color = "gray95"),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )
    
    # Add shaded areas - create separate dataframes for each sensor
    agcl_data <- subset_data[subset_data$Sensor == "AgCl", ]
    pphg_data <- subset_data[subset_data$Sensor == "PPHG", ]
    
    # Add shaded areas from curve to baseline (0)
    if (nrow(agcl_data) > 0) {
      p <- p + geom_ribbon(
        data = agcl_data,
        aes(x = Time, ymin = 0, ymax = Amplitude),
        fill = "#d9020d", 
        alpha = 0.3,
        inherit.aes = FALSE
      )
    }
    
    if (nrow(pphg_data) > 0) {
      p <- p + geom_ribbon(
        data = pphg_data,
        aes(x = Time, ymin = 0, ymax = Amplitude),
        fill = "#0066CC", 
        alpha = 0.3,
        inherit.aes = FALSE
      )
    }
    
    # Add the ECG lines on top
    p <- p + geom_line(data = subset_data, 
                      aes(x = Time, y = Amplitude, color = Sensor),
                      size = 1.2)
    
    if (nrow(subset_annotations) > 0) {
      # Add dashed lines at Q peak for reference
      p <- p + 
        geom_hline(
          data = subset_annotations,
          aes(yintercept = q_peak_amp, color = Sensor),
          linetype = "dashed",
          size = 0.5,
          alpha = 0.7
        )
      
      # P wave annotations with colored arrows from Q peak baseline
      if (all(!is.na(subset_annotations$p_wave_feature))&& show_labels) {
        # For PPHG sensor, offset arrow significantly to the left
        pphg_labels <- subset_annotations[subset_annotations$Sensor == "PPHG", ]
        if(nrow(pphg_labels) > 0) {
        p <- p + 
            geom_segment(
            data = pphg_labels,
            aes(x = p_wave_time - 0.04, y = q_peak_amp, 
                xend = p_wave_time - 0.04, yend = p_wave_amp),
            color = "#0066CC",
            arrow = arrow(type = "closed", length = unit(0.08, "inches")),
            linetype = "solid",
            size = 0.7,
            inherit.aes = FALSE
            ) +
            geom_label_repel(
            data = pphg_labels,
            aes(x = p_wave_time - 0.04, 
                y = q_peak_amp + (p_wave_amp - q_peak_amp)/2),
            label = paste0("P: ", round(pphg_labels$p_wave_feature, 1)),
            fill = "#0066CC",
            color = "white",
            fontface = "bold",
            size = 3,
            box.padding = 2.0,          # Increased box padding
            point.padding = 1.5,        # Increased point padding
            force = 30,                 # Increased repulsion force
            max.overlaps = 20,          # Increased max overlaps
            min.segment.length = 0.2,   
            label.padding = unit(0.3, "lines"),  # Increased label padding
            label.r = unit(0.2, "lines"),        # Increased corner radius
            direction = "both",         
            hjust = 1.2,
            nudge_x = -0.1,             # Increased nudge left
            segment.color = "black",  
            segment.size = 0.5,
            inherit.aes = FALSE
            )
        }
        
        # For AgCl sensor, offset arrow significantly to the right
        agcl_labels <- subset_annotations[subset_annotations$Sensor == "AgCl", ]
        if(nrow(agcl_labels) > 0) {
        p <- p + 
            geom_segment(
            data = agcl_labels,
            aes(x = p_wave_time + 0.04, y = q_peak_amp, 
                xend = p_wave_time + 0.04, yend = p_wave_amp),
            color = "#d9020d",
            arrow = arrow(type = "closed", length = unit(0.08, "inches")),
            linetype = "solid",
            size = 0.7,
            inherit.aes = FALSE
            ) +
            geom_label_repel(
            data = agcl_labels,
            aes(x = p_wave_time + 0.04, 
                y = q_peak_amp + (p_wave_amp - q_peak_amp)/2),
            label = paste0("P: ", round(agcl_labels$p_wave_feature, 1)),
            fill = "#d9020d",
            color = "white",
            fontface = "bold",
            size = 3,
            box.padding = 2.0,          # Increased box padding
            point.padding = 1.5,        # Increased point padding
            force = 30,                 # Increased repulsion force
            max.overlaps = 20,          # Increased max overlaps  
            min.segment.length = 0.2,   
            label.padding = unit(0.3, "lines"),  # Increased label padding
            label.r = unit(0.2, "lines"),        # Increased corner radius
            direction = "both",         
            hjust = 1.2,
            nudge_x = 0.1,              # Increased nudge right
            segment.color = "black",  
            segment.size = 0.5,
            inherit.aes = FALSE
            )
        }
      }
      
      # R peak annotations with horizontally offset arrows
      if (all(!is.na(subset_annotations$peak_to_peak_feature)) && show_labels) {
        pphg_labels <- subset_annotations[subset_annotations$Sensor == "PPHG", ]
        if(nrow(pphg_labels) > 0) {
        p <- p + 
            geom_segment(
            data = pphg_labels,
            aes(x = r_peak_time - 0.15, y = global_min_amp,  # Increased offset
                xend = r_peak_time, yend = global_min_amp),
            color = "#0066CC",
            linetype = "dotted",
            size = 0.8,
            inherit.aes = FALSE
            ) +
            geom_segment(
            data = pphg_labels,
            aes(x = r_peak_time - 0.04, y = global_min_amp,  # Increased offset
                xend = r_peak_time - 0.04, yend = r_peak_amp),
            color = "#0066CC",
            arrow = arrow(type = "closed", length = unit(0.08, "inches")),
            linetype = "solid",
            size = 0.7,
            inherit.aes = FALSE
            ) +
            geom_label_repel(
            data = pphg_labels,
            aes(x = r_peak_time - 0.04, 
                y = global_min_amp + (r_peak_amp - global_min_amp)*0.65),
            label = paste0("Peak-to-Peak: ", round(pphg_labels$peak_to_peak_feature, 1)),
            fill = "#0066CC",
            color = "white",
            fontface = "bold",
            size = 3,
            box.padding = 2.0,
            point.padding = 1.5,
            force = 30,
            max.overlaps = 20,
            min.segment.length = 0.2,
            label.padding = unit(0.3, "lines"),
            label.r = unit(0.2, "lines"),
            direction = "both",
            hjust = 1.2,
            nudge_x = -0.15,  # Increased nudge
            segment.color = "black",
            segment.size = 0.5,
            inherit.aes = FALSE
            )
        }
        
        agcl_labels <- subset_annotations[subset_annotations$Sensor == "AgCl", ]
        if(nrow(agcl_labels) > 0) {
        p <- p + 
            geom_segment(
            data = agcl_labels,
            aes(x = r_peak_time, y = global_min_amp,
                xend = r_peak_time + 0.15, yend = global_min_amp),  # Increased offset
            color = "#d9020d",
            linetype = "dotted",
            size = 0.8,
            inherit.aes = FALSE
            ) +
            geom_segment(
            data = agcl_labels,
            aes(x = r_peak_time + 0.04, y = global_min_amp,  # Increased offset
                xend = r_peak_time + 0.04, yend = r_peak_amp),
            color = "#d9020d",
            arrow = arrow(type = "closed", length = unit(0.08, "inches")),
            linetype = "solid",
            size = 0.7,
            inherit.aes = FALSE
            ) +
            geom_label_repel(
            data = agcl_labels,
            aes(x = r_peak_time + 0.04, 
                y = global_min_amp + (r_peak_amp - global_min_amp)*0.35),
            label = paste0("Peak-to-Peak: ", round(agcl_labels$peak_to_peak_feature, 1)),
            fill = "#d9020d",
            color = "white",
            fontface = "bold",
            size = 3,
            box.padding = 2.0,
            point.padding = 1.5,
            force = 30,
            max.overlaps = 20,
            min.segment.length = 0.2,
            label.padding = unit(0.3, "lines"),
            label.r = unit(0.2, "lines"),
            direction = "both",
            hjust = 1.2,
            nudge_x = 0.15,  # Increased nudge
            segment.color = "black",
            segment.size = 0.5,
            inherit.aes = FALSE
            )
        }
    }
      
      # T wave annotations with horizontally offset arrows
      if (all(!is.na(subset_annotations$t_wave_feature)) && show_labels) {
        pphg_labels <- subset_annotations[subset_annotations$Sensor == "PPHG", ]
        if(nrow(pphg_labels) > 0) {
        p <- p + 
            geom_segment(
            data = pphg_labels,
            aes(x = t_wave_time - 0.04, y = q_peak_amp,  # Increased offset
                xend = t_wave_time - 0.04, yend = t_wave_amp),
            color = "#0066CC",
            arrow = arrow(type = "closed", length = unit(0.08, "inches")),
            linetype = "solid",
            size = 0.7,
            inherit.aes = FALSE
            ) +
            geom_label_repel(
            data = pphg_labels,
            aes(x = t_wave_time - 0.04, 
                y = q_peak_amp + (t_wave_amp - q_peak_amp)/2),
            label = paste0("T: ", round(pphg_labels$t_wave_feature, 1)),
            fill = "#0066CC",
            color = "white",
            fontface = "bold",
            size = 3,
            box.padding = 2.0,
            point.padding = 1.5,
            force = 30,
            max.overlaps = 20,
            min.segment.length = 0.2,
            label.padding = unit(0.3, "lines"),
            label.r = unit(0.2, "lines"),
            direction = "both",
            hjust = 1.2,
            nudge_x = -0.15,  # Increased nudge
            segment.color = "black",
            segment.size = 0.5,
            inherit.aes = FALSE
            )
        }
        
        agcl_labels <- subset_annotations[subset_annotations$Sensor == "AgCl", ]
        if(nrow(agcl_labels) > 0) {
        p <- p + 
            geom_segment(
            data = agcl_labels,
            aes(x = t_wave_time + 0.04, y = q_peak_amp,  # Increased offset
                xend = t_wave_time + 0.04, yend = t_wave_amp),
            color = "#d9020d",
            arrow = arrow(type = "closed", length = unit(0.08, "inches")),
            linetype = "solid",
            size = 0.7,
            inherit.aes = FALSE
            ) +
            geom_label_repel(
            data = agcl_labels,
            aes(x = t_wave_time + 0.04, 
                y = q_peak_amp + (t_wave_amp - q_peak_amp)/2),
            label = paste0("T: ", round(agcl_labels$t_wave_feature, 1)),
            fill = "#d9020d",
            color = "white",
            fontface = "bold",
            size = 3,
            box.padding = 2.0,
            point.padding = 1.5,
            force = 30,
            max.overlaps = 20,
            min.segment.length = 0.2,
            label.padding = unit(0.3, "lines"),
            label.r = unit(0.2, "lines"),
            direction = "both",
            hjust = 1.2,
            nudge_x = 0.15,  # Increased nudge
            segment.color = "black",
            segment.size = 0.5,
            inherit.aes = FALSE
            )
        }
      }
    }
    
    plot_list[[channel]] <- p
  }
}

# Combine all plots using patchwork
if (length(plot_list) > 0) {
  combined_plot <- wrap_plots(plot_list, ncol = length(selected_channels)) + 
    plot_layout(guides = "collect") & 
    theme(legend.position = "bottom")
  
  # Create output directory if it doesn't exist
  if (!dir.exists("R_figures/figure_1d")) {
    dir.create("R_figures/figure_1d", recursive = TRUE)
  }
  
  # Display the combined plot
  print(combined_plot)
  
  # Save the combined plot
  ggsave(paste0(output_dir,"ecg_channel_kmeans_with_features.png"), combined_plot, 
         width = 15, height = 6, dpi = 300, bg = "white")
}

# Also create high-quality individual plots
for (channel_name in names(plot_list)) {
  ggsave(paste0(output_dir, channel_name, "_features.png"), 
         plot_list[[channel_name]], width = 8, height = 6, dpi = 300, bg = "white")
}