#
# Power Spectral Density Visualization
# Creates comparative PSD plots between sensors (inter-sensor) or between raw/processed signals (intra-sensor)
# Generates main PSD overview and detailed band-specific visualizations
# Supports customizable frequency bands and layout with multi-panel presentation options

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

data_dir <- '/home/giannis/Documents/ECG HG paper/results_data' #your data_dir here
fs <- 200
participant <- "p10" #c("p1","p5","p10","p39")
sensor1 <- "AgCl"
sensor2 <- "HG1"
same_sensor <- "HG1" #used if plot_type is "intra_sensor"
plot_type <- "intra_sensor"
channel <- "channel_1"
bands_to_use <- "all"  # Can be "all" or a subset like c("low", "high")
#available bands are:
#        'trend',        [0.0 0.3];    
#        'vlf',          [0.05 0.5];
#        'low',          [0.5 5];
#        'mid',          [5 15];
#        'high',         [15 40];
#        'powerline',    [48 52];
#        'powerline_harmonic', [78,82];
#        'ultra_high1',   [52 78];
#        'ultra_high2',   [82 fs/2];

# Load the PSD data
if (plot_type == "inter_sensor") {
    data_path <- file.path(data_dir,paste0("manually_cleaned_power_spectral_density_diff_sensor_data/PSD_", participant, "_", sensor1, "-", sensor2, ".json"))
} else {
    data_path <- file.path(data_dir,paste0("manually_cleaned_power_spectral_density_same_sensor_data/PSD_", participant, "_", same_sensor, ".json"))
}
data <- fromJSON(txt = data_path, simplifyVector = FALSE)

# Create output directory if it doesn't exist
output_dir <- file.path("..", "imgs_figures", "figure_S7_no_axes", participant, channel, plot_type)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

if(plot_type == "inter_sensor") {
  # Check if we have the right structure before renaming
  for(ch_name in names(data)) {
    if("Pxx_raw" %in% names(data[[ch_name]]) && "Pxx_proc" %in% names(data[[ch_name]])) {
      # Rename the fields to sensor names
      data[[ch_name]][[sensor1]] <- data[[ch_name]]$Pxx_raw
      data[[ch_name]][[sensor2]] <- data[[ch_name]]$Pxx_proc
      
      # Save the frequency values separately if needed
      if("freq" %in% names(data[[ch_name]])) {
        freq_values <- data[[ch_name]]$freq
        data[[ch_name]][[sensor1]]$freq <- freq_values
        data[[ch_name]][[sensor2]]$freq <- freq_values
      }

      data[[ch_name]]$Pxx_raw <- NULL
      data[[ch_name]]$Pxx_proc <- NULL

    }
  }
}


process_psd_data <- function(data, channel_name, sensor_name, signal_type = NULL) {
  channel_data <- data[[channel_name]]
  
  # Extract frequency values
  if("freq" %in% names(channel_data)) {
    freq <- as.numeric(channel_data$freq)
  } else if("freq" %in% names(channel_data[[sensor_name]])) {
    freq <- as.numeric(channel_data[[sensor_name]]$freq)
  } else {
    # Create frequency vector based on sampling rate if not available
    freq <- seq(0, fs/2, length.out = length(channel_data[[sensor_name]]))
  }
  
  # For inter_sensor case, sensor_name contains the PSD values directly
  if(plot_type == "inter_sensor") {
    psd <- unlist(channel_data[[sensor_name]])
    psd <- 10*log10(psd + 1e-10)
    n <- length(psd)
    half_length <- ceiling(n/2)
    psd <- psd[1:half_length]

    # Create data frame
    df <- data.frame(
      Frequency = freq,
      Power = psd,
      Sensor = sensor_name,
      Signal_Type = "inter_sensor" 
    )
    
    return(df)
  } 
  # For intra_sensor case, use raw or proc
  else {
    if (signal_type == "raw") {
      psd <- as.numeric(channel_data$Pxx_raw)
    } else {
      psd <- as.numeric(channel_data$Pxx_proc)
    }
    
    psd <- 10*log10(psd + 1e-10)

    df <- data.frame(
      Frequency = freq,
      Power = psd,
      Sensor = sensor_name,
      Signal_Type = signal_type
    )
    
    return(df)
  }
}

if(plot_type == "inter_sensor") {
  # For inter_sensor, get one dataset per sensor
  psd_sensor1_data <- process_psd_data(data, channel, sensor1)
  psd_sensor2_data <- process_psd_data(data, channel, sensor2)
  
  # Change HG1 label to PPHG
  psd_sensor2_data$Sensor <- "PPHG"
  
  # Combine data for plotting
  psd_data <- rbind(psd_sensor1_data, psd_sensor2_data)
} else {
  # For intra_sensor, get raw and processed for each sensor
  psd_raw <- process_psd_data(data, channel, same_sensor, "raw")
  psd_proc <- process_psd_data(data, channel, same_sensor, "proc")
  
  # Combine data for plotting
  psd_data <- rbind(
    psd_raw, psd_proc
  )
}

# Extract frequency bands
get_bands <- function(data, channel_name) {
  bands <- data[[channel_name]]$bands
  band_list <- list()
  
  for (band_name in names(bands)) {
    range <- bands[[band_name]]$range
    band_list[[band_name]] <- list(
      min = as.numeric(range[1]),
      max = as.numeric(range[2]),
      name = band_name
    )
  }
  
  return(band_list)
}

bands <- get_bands(data, channel)

# Determine which bands to plot
if (bands_to_use == "all") {
  bands_to_plot <- names(bands)
} else {
  bands_to_plot <- intersect(bands_to_use, names(bands))
}

if(plot_type == "inter_sensor") {
  color_palette <- c(
    "AgCl" = "#d9020d", 
    "PPHG" = "#0066CC"
  )
  color_mapping <- "Sensor"
} else {
  if (same_sensor == "AgCl") {
    color_palette <- c(
      "Raw" = "#000000",  
      "Processed" = "#d9020d"  
    )
  } else {
    #Hydrogel
    color_palette <- c(
      "Raw" = "#000000", 
      "Processed" = "#0066CC"  
    )
  }
  color_mapping <- "Signal_Type"
}



create_main_plot <- function(psd_data) {
  if(plot_type == "inter_sensor") {
    plot_title <- "Inter-Sensor PSD Comparison"
    color_var <- "Sensor"
  } else {
    plot_title <- paste0(same_sensor, " Signal Processing Comparison")
    color_var <- "Signal_Type"
    # Update labels for raw and processed
    psd_data$Signal_Type <- factor(psd_data$Signal_Type, 
                                  levels = c("raw", "proc"),
                                  labels = c("Raw", "Processed"))
  }
  
  p <- ggplot(psd_data, aes_string(x = "Frequency", y = "Power", color = color_var)) +
    geom_line(size = 1, linetype = "solid") +
    labs(
      title = "",
      subtitle = "",
      x = "Frequency (Hz)",
      y = "PSD (10*log10(μV²/Hz))",
      color = ""
    ) +
    scale_x_continuous(limits = c(0, fs/2)) +
    scale_y_continuous(limits = c(-130, 10)) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 14),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 18),
      legend.position = "bottom",
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 18),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA)
    )
  
  
  p <- p + scale_color_manual(values = color_palette)
  
  return(p)
}

create_band_plot <- function(psd_data, band_name) {
  band <- bands[[band_name]]
  
  # Filter data for this frequency band
  band_data <- psd_data[psd_data$Frequency >= band$min & psd_data$Frequency <= band$max, ]
  
  # Get the color for the plot title only
  band_color <- viridis(length(names(bands)))[which(names(bands) == band_name)]
  
  # Set up color variable based on plot type
  if(plot_type == "inter_sensor") {
    color_var <- "Sensor"
  } else {
    color_var <- "Signal_Type"
    # Update labels for raw and processed
    band_data$Signal_Type <- factor(band_data$Signal_Type, 
                                   levels = c("raw", "proc"),
                                   labels = c("Raw", "Processed"))
  }
  
  p <- ggplot(band_data, aes_string(x = "Frequency", y = "Power", color = color_var)) +
    geom_line(size = 1, linetype = "solid") +
    labs(
      title = "",
      x = "",
      y = "",
      color = ""
    ) +
    # Simple x-axis limits without custom breaks
    scale_x_continuous(limits = c(band$min, band$max)) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold", color = band_color),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      legend.position = "none",
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 18),
      legend.key.size = unit(1, "cm"),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA)
    )
  
  p <- p + scale_color_manual(values = color_palette)
  
  return(p)
}

# Generate plots based on plot_type
if (plot_type == "inter_sensor") {
  # Create main plot for inter-sensor comparison
  p_main <- create_main_plot(psd_data)
  ggsave(file.path(output_dir, paste0("PSD_", participant, "_", sensor1, "_vs_", sensor2, "_inter_sensor.png")), 
         p_main, width = 10, height = 6, dpi = 300, bg = "white")
  
  # Create band plots
  band_plots <- list()
  
  for (band_name in bands_to_plot) {
    band_plots[[band_name]] <- create_band_plot(psd_data, band_name)
    
    # Save individual band plots
    ggsave(file.path(output_dir, paste0("PSD_", participant, "_", sensor1, "_vs_", sensor2, "_", band_name, ".png")), 
           band_plots[[band_name]], width = 6, height = 4, dpi = 300, bg = "white")
  }
  
  # Create multi-panel plots
  if (length(bands_to_plot) > 0) {
    multi_panel <- p_main + 
      wrap_plots(band_plots[bands_to_plot], ncol = min(3, length(bands_to_plot))) +
      plot_layout(ncol = 1, heights = c(2, length(bands_to_plot)/3)) +
      plot_annotation(
        title = "",
        subtitle = "",
        theme = theme(
          plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 14, hjust = 0.5),
          plot.background = element_rect(fill = "white", color = NA)
        )
      ) &
      theme(
        legend.position = "none",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )
    
    ggsave(file.path(output_dir, paste0("PSD_", participant, "_", sensor1, "_vs_", sensor2, "_multi_panel.png")), 
           multi_panel, width = 12, height = 8 + 2 * ceiling(length(bands_to_plot)/3), dpi = 300, bg = "white")
  }
  
  # Display the main plot
  print(p_main)
} else {
  # Create main plot for intra-sensor comparison
  p_main <- create_main_plot(psd_data)
  ggsave(file.path(output_dir, paste0("PSD_", participant, "_", same_sensor, "_raw_vs_proc.png")), 
         p_main, width = 10, height = 6, dpi = 300, bg = "white")
  
  # Create band plots
  band_plots <- list()
  
  for (band_name in bands_to_plot) {
    band_plots[[band_name]] <- create_band_plot(psd_data, band_name)
    
    # Save individual band plots
    ggsave(file.path(output_dir, paste0("PSD_", participant, "_", same_sensor, "_", band_name, ".png")), 
           band_plots[[band_name]], width = 6, height = 4, dpi = 300, bg = "white")
  }
  
  # Create multi-panel plots
  if (length(bands_to_plot) > 0) {
    multi_panel <- p_main + 
      wrap_plots(band_plots[bands_to_plot], ncol = min(3, length(bands_to_plot))) +
      plot_layout(ncol = 1, heights = c(2, length(bands_to_plot)/3)) +
      plot_annotation(
        title = "",
        subtitle = "",
        theme = theme(
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 14, hjust = 0.5),
          plot.background = element_rect(fill = "white", color = NA)
        )
      ) &
      theme(
        legend.position = "none",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )
    
    ggsave(file.path(output_dir, paste0("PSD_", participant, "_", same_sensor, "_multi_panel.png")), 
           multi_panel, width = 12, height = 8 + 2 * ceiling(length(bands_to_plot)/3), dpi = 300, bg = "white")
  }
  
  # Display the main plot
  print(p_main)
}