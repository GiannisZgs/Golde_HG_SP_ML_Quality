#
# ECG Biofouling Analysis Visualization
# Plots segmented ECG signals from multiple recordings with SNR annotations
# Maintains consistent styling and creates individual plots for each recording
# Examines possible signal quality degradation over time using SNR values

library(jsonlite)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(R.matlab)
library(scales)
library(grid)
library(cowplot)

lead_to_plot <- 3  # Choose lead 1, 2, or 3

data_dir <- file.path('.', 'data')
data <- fromJSON(txt = file.path(data_dir,'ECG_biofouling_data.json'), simplifyVector = FALSE)

# Create output directory
output_dir <- file.path(".", "imgs_figures", "figure_S13")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Function to create individual recording plots
create_recording_plot <- function(recording_data, recording_num, lead_to_plot) {
  # Extract the segmented ECG and time axis
  ecg_seg <- unlist(recording_data$ECG_seg)
  time_seg <- unlist(recording_data$t_axis_seg)
  snr_value <- recording_data$SNR[[lead_to_plot]]
  
  # Subtract time offset to start from 0
  time_seg <- time_seg - min(time_seg)
  
  # Create data frame
  plot_data <- data.frame(
    Time = time_seg,
    Amplitude = ecg_seg
  )
  
  # Calculate y-axis limits for consistent scaling
  y_min <- min(plot_data$Amplitude)
  y_max <- max(plot_data$Amplitude)
  y_range <- y_max - y_min
  
  # Create the plot with existing styling
  p <- ggplot(plot_data, aes(x = Time, y = Amplitude)) +
    geom_line(size = 0.3, color = "#0066CC") +
    labs(
      title = paste("Recording", recording_num, "- SNR =", round(snr_value, 2)),
      x = "Time (s)",
      y = "Amplitude (\u00B5V)"
    ) +
    theme_minimal() + 
    theme(
      plot.title = element_text(size = 35, hjust = 0.5, face = "bold"),
      axis.title = element_text(size = 30),
      axis.text = element_text(size = 30, color = "black"),
      axis.line = element_line(color = "black", size = 1), 
      axis.ticks = element_line(color = "black", size = 0.8),  
      axis.ticks.length = unit(0.3, "cm"),  
      panel.grid.major = element_line(color = "gray60"),  
      panel.grid.minor = element_line(color = "gray80"),  
      panel.background = element_rect(fill = "white", color = NA), 
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  return(p)
}

# Create plots for all 5 recordings
plot_list <- list()
for (i in 1:5) {
  if (!is.null(data[[i]])) {
    plot_list[[i]] <- create_recording_plot(data[[i]], i, lead_to_plot)
  }
}

# Create individual plots and save them
for (i in 1:length(plot_list)) {
  if (!is.null(plot_list[[i]])) {
    ggsave(file.path(output_dir, paste0("recording_", i, "_lead_", lead_to_plot, ".png")), 
           plot_list[[i]], width = 10, height = 6, dpi = 300, bg = "white")
  }
}

# Create a combined plot with all recordings
if (length(plot_list) > 0) {
  # Create a layout that centers the 5th plot in the bottom row
  # First create the first 4 plots in a 2x2 grid
  top_plots <- wrap_plots(plot_list[1:4], ncol = 2, nrow = 2)
  
  # Create a bottom row with the 5th plot centered
  if (length(plot_list) >= 5) {
    # Create a spacer and the 5th plot with specific widths
    empty_plot <- ggplot() + theme_void()
    
    # Use plot_layout with specific widths: 1 unit spacer, 2 units for plot, 1 unit spacer
    # This makes the 5th plot twice as wide as each spacer, matching the width of the plots above
    bottom_row <- wrap_plots(list(empty_plot, plot_list[[5]], empty_plot), ncol = 3, nrow = 1) +
      plot_layout(widths = c(1, 2, 1))
    
    # Combine top and bottom sections
    combined_plot <- top_plots / bottom_row + 
      plot_layout(heights = c(2, 1)) +
      plot_annotation(
        title = paste("Lead", lead_to_plot),
        theme = theme(
          plot.title = element_text(size = 35, hjust = 0.5, face = "bold")
        )
      )
  } else {
    # If less than 5 plots, use original layout
    combined_plot <- wrap_plots(plot_list, ncol = 2, nrow = 3) +
      plot_annotation(
        title = paste("Lead", lead_to_plot),
        theme = theme(
          plot.title = element_text(size = 35, hjust = 0.5, face = "bold")
        )
      )
  }
  
  # Save the combined plot
  ggsave(file.path(output_dir, paste0("combined_recordings_lead_", lead_to_plot, ".png")), 
         combined_plot, width = 20, height = 18, dpi = 300, bg = "white")
  
  # Display the combined plot
  print(combined_plot)
}

cat("Plots saved to:", output_dir, "\n")
cat("Lead analyzed:", lead_to_plot, "\n")